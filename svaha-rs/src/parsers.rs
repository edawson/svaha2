use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use crate::{Variant, PackedSequence, MafEntry};
use anyhow::Result;

pub struct VcfParser;

impl VcfParser {
    pub fn parse<P, F>(path: P, mut on_variant: F) -> Result<()>
    where
        P: AsRef<Path>,
        F: FnMut(Variant),
    {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with("##") && !line.starts_with("##CHROM") {
                continue;
            }
            if line.starts_with('#') {
                continue;
            }

            let tokens: Vec<&str> = line.split('\t').collect();
            if tokens.len() < 8 {
                continue;
            }

            let chrom = tokens[0].to_string();
            let pos = tokens[1].parse::<u64>()? - 1;
            let ref_seq = tokens[3].to_string();
            let alts = tokens[4];
            let info = tokens[7];

            for alt in alts.split(',') {
                let mut end = pos + ref_seq.len() as u64;
                let mut chrom_2 = chrom.clone();
                let mut pos_2 = end;

                let variant_type = if info.contains("SVTYPE=DEL") { "DEL".to_string() }
                else if info.contains("SVTYPE=INS") { "INS".to_string() }
                else if info.contains("SVTYPE=INV") { "INV".to_string() }
                else if info.contains("SVTYPE=BND") || info.contains("SVTYPE=TRA") { "BND".to_string() }
                else if ref_seq.len() == 1 && alt.len() == 1 { "SNP".to_string() }
                else { "INDEL".to_string() };

                if let Some(start_idx) = info.find("END=") {
                    let rest = &info[start_idx + 4..];
                    let end_idx = rest.find(';').unwrap_or(rest.len());
                    if let Ok(e) = rest[..end_idx].parse::<u64>() {
                        end = e;
                    }
                }

                if variant_type == "BND" {
                    if let Some(start_idx) = info.find("CHR2=") {
                        let rest = &info[start_idx + 5..];
                        let end_idx = rest.find(';').unwrap_or(rest.len());
                        chrom_2 = rest[..end_idx].to_string();
                    }
                    if let Some(start_idx) = info.find("POS2=") {
                        let rest = &info[start_idx + 5..];
                        let end_idx = rest.find(';').unwrap_or(rest.len());
                        if let Ok(p) = rest[..end_idx].parse::<u64>() {
                            pos_2 = p;
                        }
                    }
                }

                // Handle anchor bases for small variants
                let mut v_pos = pos;
                let mut v_end = end;
                let mut v_alt = alt.to_string();
                if (variant_type == "SNP" || variant_type == "INDEL") && !ref_seq.is_empty() && !alt.is_empty() {
                    if ref_seq.starts_with(alt.chars().next().unwrap()) && alt.starts_with(ref_seq.chars().next().unwrap()) {
                        v_pos += 1;
                        v_alt = alt[1..].to_string();
                        if ref_seq.len() == 1 {
                            v_end = v_pos; // Insertion
                        } else {
                            v_end = v_pos + ref_seq.len() as u64 - 1;
                        }
                    }
                }

                on_variant(Variant {
                    chrom: chrom.clone(),
                    pos: v_pos,
                    end: v_end,
                    ref_seq: ref_seq.clone(),
                    packed_alt: PackedSequence::pack(&v_alt),
                    alt_len: v_alt.len() as u32,
                    variant_type,
                    chrom_2: Some(chrom_2),
                    pos_2: Some(pos_2),
                    from_forward: true,
                    to_forward: true,
                });
            }
        }
        Ok(())
    }
}

pub struct CbioMafParser;

impl CbioMafParser {
    pub fn parse<P, F>(path: P, mut on_variant: F) -> Result<()>
    where
        P: AsRef<Path>,
        F: FnMut(Variant),
    {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut header: Option<HashMap<String, usize>> = None;

        for line in reader.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let tokens: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();

            if header.is_none() {
                let mut h = HashMap::new();
                for (i, token) in tokens.iter().enumerate() {
                    h.insert(token.clone(), i);
                }
                header = Some(h);
                continue;
            }

            let h = header.as_ref().unwrap();
            if tokens.len() < h.len() {
                continue;
            }

            let chrom = tokens[h["Chromosome"]].clone();
            let pos = tokens[h["Start_Position"]].parse::<u64>()? - 1;
            let mut end = tokens[h["End_Position"]].parse::<u64>()?;
            let ref_seq = tokens[h["Reference_Allele"]].clone();
            let mut alt = tokens[h["Tumor_Seq_Allele2"]].clone();
            let variant_type = tokens[h["Variant_Type"]].clone();

            if alt == "-" {
                alt = String::new();
            }
            let packed_alt = PackedSequence::pack(&alt);
            let alt_len = alt.len() as u32;

            if variant_type == "SNP" && end == pos {
                end = pos + 1;
            }

            on_variant(Variant {
                chrom,
                pos,
                end,
                ref_seq,
                packed_alt,
                alt_len,
                variant_type,
                chrom_2: None,
                pos_2: None,
                from_forward: true,
                to_forward: true,
            });
        }
        Ok(())
    }
}

pub struct SvParser;

impl SvParser {
    pub fn parse<P, F>(path: P, mut on_variant: F) -> Result<()>
    where
        P: AsRef<Path>,
        F: FnMut(Variant),
    {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut header: Option<HashMap<String, usize>> = None;

        for line in reader.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let tokens: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();

            if header.is_none() {
                let mut h = HashMap::new();
                for (i, token) in tokens.iter().enumerate() {
                    h.insert(token.clone(), i);
                }
                header = Some(h);
                continue;
            }

            let h = header.as_ref().unwrap();
            if tokens.len() < h.len() {
                continue;
            }

            let chrom = tokens[h["Site1_Chromosome"]].clone();
            let mut pos = tokens[h["Site1_Position"]].parse::<u64>()? - 1;
            let chrom_2 = tokens[h["Site2_Chromosome"]].clone();
            let pos_2 = tokens[h["Site2_Position"]].parse::<u64>()? - 1;
            let variant_type = tokens[h["Class"]].clone();

            let conn = &tokens[h["Connection_Type"]];
            let (from_forward, to_forward) = match conn.as_str() {
                "3to5" => (true, true),
                "5to5" => (false, true),
                "3to3" => (true, false),
                "5to3" => (false, false),
                _ => (true, true),
            };

            let mut end = pos + 1;
            if chrom == chrom_2 {
                end = pos_2;
                if end < pos {
                    std::mem::swap(&mut pos, &mut end);
                }
            }

            on_variant(Variant {
                chrom,
                pos,
                end,
                ref_seq: String::new(),
                packed_alt: Vec::new(),
                alt_len: 0,
                variant_type,
                chrom_2: Some(chrom_2),
                pos_2: Some(pos_2),
                from_forward,
                to_forward,
            });
        }
        Ok(())
    }
}

pub struct MafParser;

impl MafParser {
    pub fn parse<P, F>(path: P, mut on_block: F) -> Result<()>
    where
        P: AsRef<Path>,
        F: FnMut(MafEntry, MafEntry),
    {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        while let Some(line) = lines.next() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') || line.starts_with('a') {
                continue;
            }

            if line.starts_with('s') {
                let ref_entry = Self::parse_s_line(&line)?;
                if let Some(next_line) = lines.next() {
                    let next_line = next_line?;
                    if next_line.starts_with('s') {
                        let alt_entry = Self::parse_s_line(&next_line)?;
                        on_block(ref_entry, alt_entry);
                    }
                }
            }
        }
        Ok(())
    }

    fn parse_s_line(line: &str) -> Result<MafEntry> {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() < 7 {
            anyhow::bail!("Invalid s line: {}", line);
        }
        Ok(MafEntry {
            src: tokens[1].to_string(),
            start: tokens[2].parse()?,
            size: tokens[3].parse()?,
            strand: tokens[4].chars().next().unwrap_or('+'),
            src_size: tokens[5].parse()?,
            text: tokens[6].to_string(),
        })
    }
}

pub fn detect_format<P: AsRef<Path>>(path: P) -> Result<String> {
    let file = File::open(&path)?;
    let mut reader = BufReader::new(file);
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;

    if first_line.contains("Sample_Id") && (first_line.contains("Connection_Type") || first_line.contains("SV_Status")) {
        Ok("sv".to_string())
    } else if first_line.contains("Hugo_Symbol") || first_line.contains("#genome_nexus") {
        Ok("cbio_maf".to_string())
    } else if path.as_ref().extension().map_or(false, |ext| ext == "maf") && !first_line.is_empty() && first_line.starts_with('#') {
        Ok("maf".to_string())
    } else {
        Ok("vcf".to_string())
    }
}
