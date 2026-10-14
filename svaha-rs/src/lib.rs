use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::BufReader;
use memmap2::Mmap;
use noodles::fasta;
use anyhow::Result;
use std::cmp::min;

/// rGFA Node representation with Stable Coordinates
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Node {
    pub id: u64,
    pub sn: String,   // Stable Name (e.g., "chr1")
    pub so: u64,      // Stable Offset (0-based)
    pub sr: u32,      // Stable Rank (0=ref, 1+=alt)
    pub length: u32,
}

impl Node {
    pub fn to_gfa_s(&self, sequence: &str) -> String {
        format!(
            "S\t{}\t{}\tSN:Z:{}\tSO:i:{}\tSR:i:{}",
            self.id, sequence, self.sn, self.so, self.sr
        )
    }
}

/// Edge representation using IDs
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Edge {
    pub from: u64,
    pub to: u64,
    pub from_forward: bool,
    pub to_forward: bool,
}

impl Edge {
    pub fn to_gfa_l(&self) -> String {
        format!(
            "L\t{}\t{}\t{}\t{}\t0M",
            self.from,
            if self.from_forward { "+" } else { "-" },
            self.to,
            if self.to_forward { "+" } else { "-" }
        )
    }
}

/// Variant representation for Phase 0
#[derive(Debug, Clone)]
pub struct Variant {
    pub chrom: String,
    pub pos: u64,           // 0-based start
    pub end: u64,           // 0-based end
    pub ref_seq: String,
    pub packed_alt: Vec<u64>, // 2-bit packed alt sequence
    pub alt_len: u32,
    pub variant_type: String, // SNP, INS, DEL, INV, etc.
    pub chrom_2: Option<String>, // For translocations
    pub pos_2: Option<u64>,      // For translocations
    pub from_forward: bool,
    pub to_forward: bool,
}

/// Genomic Region representation
#[derive(Debug, Clone)]
pub struct Region {
    pub chrom: String,
    pub start: u64, // 0-based
    pub end: u64,   // 0-based
}

impl Region {
    pub fn parse(s: &str) -> anyhow::Result<Self> {
        let colon = s.find(':');
        match colon {
            None => Ok(Region {
                chrom: s.to_string(),
                start: 0,
                end: u64::MAX,
            }),
            Some(c) => {
                let chrom = s[..c].to_string();
                let rest = &s[c + 1..];
                let dash = rest.find('-');
                match dash {
                    None => {
                        let start = rest.parse::<u64>()?;
                        Ok(Region { chrom, start, end: u64::MAX })
                    }
                    Some(d) => {
                        let start = rest[..d].parse::<u64>()?;
                        let end = rest[d + 1..].parse::<u64>()?;
                        Ok(Region { chrom, start, end })
                    }
                }
            }
        }
    }

    pub fn contains(&self, c: &str, p: u64) -> bool {
        c == self.chrom && p >= self.start && p < self.end
    }

    pub fn overlaps(&self, c: &str, s: u64, e: u64) -> bool {
        c == self.chrom && s < self.end && e > self.start
    }
}

/// Simple Bit-Vector for breakpoints (Phase 0/1)
#[derive(Default)]
pub struct BreakpointVector {
    breakpoints: BTreeSet<u64>,
}

impl BreakpointVector {
    pub fn add(&mut self, pos: u64) {
        self.breakpoints.insert(pos);
    }

    pub fn get(&self) -> Vec<u64> {
        self.breakpoints.iter().cloned().collect()
    }

    pub fn size(&self) -> usize {
        self.breakpoints.len()
    }
}

/// 2-bit sequence packer
pub struct PackedSequence;

pub struct MafEntry {
    pub src: String,
    pub start: u64,
    pub size: u64,
    pub strand: char,
    pub src_size: u64,
    pub text: String,
}

pub mod parsers;

impl PackedSequence {
    pub fn pack(seq: &str) -> Vec<u64> {
        if seq == "-" || seq.is_empty() {
            return Vec::new();
        }
        let n_blocks = (seq.len() + 31) / 32;
        let mut packed = vec![0u64; n_blocks];
        for (i, c) in seq.chars().enumerate() {
            let val = match c.to_ascii_uppercase() {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' => 3,
                _ => 0, // Treat N as A for 2-bit
            };
            packed[i / 32] |= val << (2 * (i % 32));
        }
        packed
    }

    pub fn unpack(packed: &[u64], length: usize) -> String {
        if length == 0 {
            return String::new();
        }
        let mut seq = String::with_capacity(length);
        for i in 0..length {
            let val = (packed[i / 32] >> (2 * (i % 32))) & 0x3;
            let c = match val {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => unreachable!(),
            };
            seq.push(c);
        }
        seq
    }
}

pub fn phase_0_variant_parsing(
    var_path: &str,
    ref_path: &str,
    chrom_to_bps: &mut HashMap<String, BreakpointVector>,
    variants: &mut Vec<Variant>,
    region: &Option<Region>,
    explicit_format: &str,
    relax_chrom: bool,
) -> Result<()> {
    let fai_path = format!("{}.fai", ref_path);
    let mut fai_reader = File::open(&fai_path)
        .map(BufReader::new)
        .map(fasta::fai::io::Reader::new)?;
    let index = fai_reader.read_index()?;

    let normalize_chrom = |c: &str| -> String {
        if index.as_ref().iter().any(|r| r.name() == c.as_bytes()) {
            return c.to_string();
        }
        if relax_chrom {
            if c.starts_with("chr") {
                let no_chr = &c[3..];
                if index.as_ref().iter().any(|r| r.name() == no_chr.as_bytes()) {
                    return no_chr.to_string();
                }
            }
            let with_chr = format!("chr{}", c);
            if index.as_ref().iter().any(|r| r.name() == with_chr.as_bytes()) {
                return with_chr;
            }
        }
        c.to_string()
    };

    let format = if explicit_format.is_empty() {
        parsers::detect_format(var_path)?
    } else {
        explicit_format.to_string()
    };

    let mut process_variant = |mut v: Variant| {
        let orig_chrom = v.chrom.clone();
        v.chrom = normalize_chrom(&v.chrom);
        if let Some(ref mut c2) = v.chrom_2 {
            *c2 = normalize_chrom(c2);
        }

        if let Some(r) = region {
            if !r.overlaps(&v.chrom, v.pos, v.end) {
                return;
            }
        }

        if !index.as_ref().iter().any(|r| r.name() == v.chrom.as_bytes()) {
            eprintln!("Error: Chromosome '{}' not found in reference.", orig_chrom);
            if !relax_chrom {
                std::process::exit(1);
            }
            return;
        }

        chrom_to_bps.entry(v.chrom.clone()).or_default().add(v.pos);
        chrom_to_bps.entry(v.chrom.clone()).or_default().add(v.end);

        if let Some(ref c2) = v.chrom_2 {
            if &v.chrom != c2 {
                if index.as_ref().iter().any(|r| r.name() == c2.as_bytes()) {
                    let pos_2 = v.pos_2.unwrap_or(0);
                    chrom_to_bps.entry(c2.clone()).or_default().add(pos_2);
                    chrom_to_bps.entry(c2.clone()).or_default().add(pos_2 + 1);
                } else if !relax_chrom {
                    eprintln!("Error: Chromosome '{}' not found in reference.", c2);
                    std::process::exit(1);
                }
            }
        }
        variants.push(v);
    };

    match format.as_str() {
        "sv" | "bedpe" => parsers::SvParser::parse(var_path, &mut process_variant)?,
        "cbio_maf" => parsers::CbioMafParser::parse(var_path, &mut process_variant)?,
        "maf" => parsers::MafParser::parse(var_path, |ref_entry, alt_entry| {
            let clean_alt: String = alt_entry.text.chars().filter(|&c| c != '-').collect();
            let v = Variant {
                chrom: ref_entry.src.clone(),
                pos: ref_entry.start,
                end: ref_entry.start + ref_entry.size,
                ref_seq: String::new(),
                packed_alt: PackedSequence::pack(&clean_alt),
                alt_len: clean_alt.len() as u32,
                variant_type: "MAF_BLOCK".to_string(),
                chrom_2: Some(ref_entry.src),
                pos_2: Some(ref_entry.start + ref_entry.size),
                from_forward: true,
                to_forward: true,
            };
            process_variant(v);
        })?,
        _ => parsers::VcfParser::parse(var_path, &mut process_variant)?,
    }

    Ok(())
}

pub fn phase_1_topology(
    chrom_to_bps: &HashMap<String, BreakpointVector>,
    variants: &[Variant],
    nodes: &mut Vec<Node>,
    edges: &mut Vec<Edge>,
    pos_to_node_id: &mut HashMap<(String, u64), u64>,
    end_to_node_id: &mut HashMap<(String, u64), u64>,
) {
    let mut node_id = 1;

    // First, process backbone nodes
    let mut sorted_chroms: Vec<_> = chrom_to_bps.keys().collect();
    sorted_chroms.sort();

    for chrom in sorted_chroms {
        let bps = chrom_to_bps[chrom].get();
        for i in 0..bps.len().saturating_sub(1) {
            let n = Node {
                id: node_id,
                sn: chrom.clone(),
                so: bps[i],
                sr: 0,
                length: (bps[i + 1] - bps[i]) as u32,
            };
            node_id += 1;
            pos_to_node_id.insert((chrom.clone(), bps[i]), n.id);
            end_to_node_id.insert((chrom.clone(), bps[i + 1]), n.id);
            nodes.push(n);
            if i > 0 {
                edges.push(Edge {
                    from: nodes[nodes.len() - 2].id,
                    to: nodes.last().unwrap().id,
                    from_forward: true,
                    to_forward: true,
                });
            }
        }
    }

    let backbone_nodes = nodes.len();

    // Then, process variants
    for v in variants {
        match v.variant_type.as_str() {
            "SNP" | "INDEL" | "INS" | "MAF_BLOCK" | "DNP" | "TNP" | "ONP" => {
                let alt = Node {
                    id: node_id,
                    sn: v.chrom.clone(),
                    so: v.pos,
                    sr: 1,
                    length: v.alt_len,
                };
                node_id += 1;
                let pre_id = end_to_node_id.get(&(v.chrom.clone(), v.pos)).cloned();
                let post_id = pos_to_node_id.get(&(v.chrom.clone(), v.end)).cloned();
                nodes.push(alt.clone());
                if let Some(pid) = pre_id {
                    edges.push(Edge {
                        from: pid,
                        to: alt.id,
                        from_forward: true,
                        to_forward: true,
                    });
                }
                if let Some(pid) = post_id {
                    edges.push(Edge {
                        from: alt.id,
                        to: pid,
                        from_forward: true,
                        to_forward: true,
                    });
                }
            }
            "DEL" | "DELETION" => {
                let pre_id = end_to_node_id.get(&(v.chrom.clone(), v.pos)).cloned();
                let post_id = pos_to_node_id.get(&(v.chrom.clone(), v.end)).cloned();
                if let (Some(pre), Some(post)) = (pre_id, post_id) {
                    edges.push(Edge {
                        from: pre,
                        to: post,
                        from_forward: true,
                        to_forward: true,
                    });
                }
            }
            "INV" | "INVERSION" => {
                let pre_id = end_to_node_id.get(&(v.chrom.clone(), v.pos)).cloned();
                let post_id = pos_to_node_id.get(&(v.chrom.clone(), v.end)).cloned();
                let start_node_id = pos_to_node_id.get(&(v.chrom.clone(), v.pos)).cloned();
                let end_node_id = end_to_node_id.get(&(v.chrom.clone(), v.end)).cloned();

                if !v.from_forward || !v.to_forward {
                    if let (Some(pre), Some(post)) = (pre_id, post_id) {
                        edges.push(Edge {
                            from: pre,
                            to: post,
                            from_forward: v.from_forward,
                            to_forward: v.to_forward,
                        });
                    }
                } else {
                    if let (Some(pre), Some(end_id)) = (pre_id, end_node_id) {
                        edges.push(Edge {
                            from: pre,
                            to: end_id,
                            from_forward: true,
                            to_forward: false,
                        });
                    }
                    if let (Some(start_id), Some(post)) = (start_node_id, post_id) {
                        edges.push(Edge {
                            from: start_id,
                            to: post,
                            from_forward: false,
                            to_forward: true,
                        });
                    }
                }
            }
            "BND" | "TRANSLOCATION" | "FUSION" => {
                let pre_id = end_to_node_id.get(&(v.chrom.clone(), v.pos)).cloned();
                let c2 = v.chrom_2.as_ref().unwrap_or(&v.chrom);
                let p2 = v.pos_2.unwrap_or(0);
                let post_id = pos_to_node_id.get(&(c2.clone(), p2)).cloned();
                if let (Some(pre), Some(post)) = (pre_id, post_id) {
                    edges.push(Edge {
                        from: pre,
                        to: post,
                        from_forward: v.from_forward,
                        to_forward: v.to_forward,
                    });
                }
            }
            "DUPLICATION" => {
                let start_id = pos_to_node_id.get(&(v.chrom.clone(), v.pos)).cloned();
                let end_id = end_to_node_id.get(&(v.chrom.clone(), v.end)).cloned();
                if let (Some(start), Some(end)) = (start_id, end_id) {
                    edges.push(Edge {
                        from: end,
                        to: start,
                        from_forward: true,
                        to_forward: true,
                    });
                }
            }
            _ => {}
        }
    }

    eprintln!(
        "  Done. Created {} backbone nodes and {} alternate nodes.",
        backbone_nodes,
        nodes.len() - backbone_nodes
    );
    eprintln!("  Total edges: {}", edges.len());
}

pub fn fetch_ref_seq(
    mmap: &Mmap,
    index: &fasta::fai::Index,
    chrom: &str,
    start: u64,
    length: u32,
) -> String {
    let entry = match index.as_ref().iter().find(|r| r.name() == chrom.as_bytes()) {
        Some(e) => e,
        None => return "N".repeat(length as usize),
    };

    let mut seq = String::with_capacity(length as usize);
    let mut remaining = length as u64;
    let mut cur_pos = start;

    while remaining > 0 {
        let line_offset = entry.offset()
            + (cur_pos / entry.line_bases()) * entry.line_width()
            + (cur_pos % entry.line_bases());
        let bases_in_line = entry.line_bases() - (cur_pos % entry.line_bases());
        let to_read = min(remaining, bases_in_line);

        if line_offset + to_read > mmap.len() as u64 {
            break;
        }

        let chunk = &mmap[line_offset as usize..(line_offset + to_read) as usize];
        seq.push_str(std::str::from_utf8(chunk).unwrap_or(""));
        remaining -= to_read;
        cur_pos += to_read;
    }
    seq
}

use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

pub fn phase_2_materialize(
    ref_path: &str,
    nodes: &[Node],
    edges: &[Edge],
    variants: &[Variant],
) -> Result<()> {
    eprintln!("Phase 2: Materializing GFA...");
    let file = File::open(ref_path)?;
    let mmap = unsafe { Mmap::map(&file)? };

    let fai_path = format!("{}.fai", ref_path);
    let mut fai_reader = File::open(&fai_path)
        .map(BufReader::new)
        .map(fasta::fai::io::Reader::new)?;
    let index = fai_reader.read_index()?;

    let mut alt_map = HashMap::new();
    for v in variants {
        alt_map.insert((v.chrom.clone(), v.pos), v);
    }

    println!("H\tVN:Z:1.1");

    let nodes_processed = AtomicUsize::new(0);
    let node_strings: Vec<String> = nodes
        .par_iter()
        .map(|n| {
            let seq = if n.sr == 0 {
                fetch_ref_seq(&mmap, &index, &n.sn, n.so, n.length)
            } else {
                match alt_map.get(&(n.sn.clone(), n.so)) {
                    Some(v) => PackedSequence::unpack(&v.packed_alt, v.alt_len as usize),
                    None => "*".to_string(),
                }
            };

            let processed = nodes_processed.fetch_add(1, Ordering::Relaxed);
            if (processed + 1) % 1_000_000 == 0 {
                eprintln!("  Processed {} / {} nodes...", processed + 1, nodes.len());
            }

            n.to_gfa_s(&seq)
        })
        .collect();

    for s in node_strings {
        println!("{}", s);
    }

    eprintln!("  Writing edges...");
    for e in edges {
        println!("{}", e.to_gfa_l());
    }
    eprintln!("Phase 2 Complete.");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_region() {
        let r = Region::parse("chr1:100-200").unwrap();
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, 100);
        assert_eq!(r.end, 200);
        assert!(r.contains("chr1", 150));
        assert!(!r.contains("chr1", 250));
        assert!(r.overlaps("chr1", 50, 150));
        assert!(!r.overlaps("chr1", 250, 300));
    }

    #[test]
    fn test_packed_sequence() {
        let seq = "ACGTACGT";
        let packed = PackedSequence::pack(seq);
        let unpacked = PackedSequence::unpack(&packed, seq.len());
        assert_eq!(seq, unpacked);
    }

    #[test]
    fn test_breakpoint_vector() {
        let mut bv = BreakpointVector::default();
        bv.add(100);
        bv.add(50);
        bv.add(100);
        let bps = bv.get();
        assert_eq!(bps.len(), 2);
        assert_eq!(bps[0], 50);
        assert_eq!(bps[1], 100);
    }
}
