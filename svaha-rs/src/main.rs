use clap::{Parser, Subcommand};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use svaha_rs::*;
use anyhow::Result;
use noodles::fasta;

#[derive(Parser)]
#[command(name = "svaha")]
#[command(about = "High-performance variation graph construction", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Construct a GFA from reference and variants
    Build {
        /// Reference FASTA (required)
        #[arg(short, long)]
        ref_fasta: String,

        /// VCF variants
        #[arg(short, long)]
        vcf: Vec<String>,

        /// Mutation Annotation Format (cBioPortal)
        #[arg(long)]
        maf: Vec<String>,

        /// Structural Variants (cBioPortal data_sv.txt)
        #[arg(long)]
        sv: Vec<String>,

        /// BEDPE structural variants
        #[arg(long)]
        bedpe: Vec<String>,

        /// Max node size (default 32)
        #[arg(short, long, default_value_t = 32)]
        max_node_size: u32,

        /// Region (chrom:start-end)
        #[arg(short = 'R', long)]
        region: Option<String>,

        /// Allow 'chr1' to match '1' and vice-versa
        #[arg(long)]
        relax_chrom: bool,
    },
    /// Show statistics for a GFA file
    Stats {
        /// GFA file
        gfa: String,
    },
    /// Visualize a portion of the graph
    View {
        /// GFA file
        gfa: String,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Build {
            ref_fasta,
            vcf,
            maf,
            sv,
            bedpe,
            max_node_size,
            region,
            relax_chrom,
        } => {
            let mut var_inputs = Vec::new();
            for v in vcf { var_inputs.push((v, "vcf")); }
            for m in maf { var_inputs.push((m, "cbio_maf")); }
            for s in sv { var_inputs.push((s, "sv")); }
            for b in bedpe { var_inputs.push((b, "bedpe")); }

            if var_inputs.is_empty() || ref_fasta.is_empty() {
                anyhow::bail!("Error: Variant input and -r (Reference) are required for build.");
            }

            let region_parsed = region.as_ref().map(|s| Region::parse(s)).transpose()?;

            let mut chrom_to_bps = HashMap::new();
            let mut variants = Vec::new();
            for (path, format) in var_inputs {
                phase_0_variant_parsing(
                    &path,
                    &ref_fasta,
                    &mut chrom_to_bps,
                    &mut variants,
                    &region_parsed,
                    format,
                    relax_chrom,
                )?;
            }

            eprintln!("Finalizing breakpoints and discretizing reference...");
            let fai_path = format!("{}.fai", ref_fasta);
            let mut fai_reader = File::open(&fai_path)
                .map(BufReader::new)
                .map(fasta::fai::io::Reader::new)?;
            let index = fai_reader.read_index()?;

            let target_chroms: Vec<String> = if let Some(ref r) = region_parsed {
                vec![r.chrom.clone()]
            } else {
                chrom_to_bps.keys().cloned().collect()
            };

            for c in target_chroms {
                let entry = index.as_ref().iter().find(|r| r.name() == c.as_bytes());
                let c_start = region_parsed.as_ref().map(|r| r.start).unwrap_or(0);
                let c_end = region_parsed.as_ref().map(|r| r.end).unwrap_or(entry.map(|e| e.length()).unwrap_or(0));
                
                let bps = chrom_to_bps.entry(c.clone()).or_default();
                bps.add(c_start);
                bps.add(c_end);
                
                let mut b = c_start + max_node_size as u64;
                while b < c_end {
                    bps.add(b);
                    b += max_node_size as u64;
                }
            }

            // Note: BreakpointVector::get returns sorted and unique values
            // In C++ version we call finalize() which does sort and unique.
            // Our BreakpointVector uses BTreeSet, so it's already sorted and unique.

            let mut nodes = Vec::new();
            let mut edges = Vec::new();
            let mut pos_to_node_id = HashMap::new();
            let mut end_to_node_id = HashMap::new();
            
            phase_1_topology(
                &chrom_to_bps,
                &variants,
                &mut nodes,
                &mut edges,
                &mut pos_to_node_id,
                &mut end_to_node_id,
            );
            
            phase_2_materialize(&ref_fasta, &nodes, &edges, &variants)?;
        }
        Commands::Stats { gfa } => {
            let file = File::open(gfa)?;
            let reader = BufReader::new(file);
            let mut node_count = 0;
            let mut edge_count = 0;
            let mut total_len = 0;

            for line in reader.lines() {
                let line = line?;
                if line.is_empty() { continue; }
                if line.starts_with('S') {
                    node_count += 1;
                    let tokens: Vec<&str> = line.split('\t').collect();
                    if tokens.len() > 2 && tokens[2] != "*" {
                        total_len += tokens[2].len();
                    }
                } else if line.starts_with('L') {
                    edge_count += 1;
                }
            }
            println!("Nodes: {}\nEdges: {}\nTotal Sequence Length: {} bp", node_count, edge_count, total_len);
        }
        Commands::View { gfa } => {
            let file = File::open(gfa)?;
            let reader = BufReader::new(file);
            
            println!("digraph G {{");
            println!("  rankdir=LR;");
            println!("  node [shape=box, style=filled, fillcolor=white, fontname=\"Courier\"];");

            for line in reader.lines() {
                let line = line?;
                if line.is_empty() { continue; }
                let tokens: Vec<&str> = line.split('\t').collect();
                if tokens[0] == "S" {
                    let id = tokens[1];
                    let seq = tokens[2];
                    let mut color = "white";
                    
                    // Look for SR:i:0
                    for token in &tokens[3..] {
                        if token.starts_with("SR:i:") && *token != "SR:i:0" {
                            color = "lightblue";
                            break;
                        }
                    }
                    
                    let mut label = id.to_string();
                    let display_seq = if seq.len() > 20 {
                        format!("{}...", &seq[..17])
                    } else {
                        seq.to_string()
                    };
                    if seq != "*" {
                        label = format!("{}\n{}", label, display_seq);
                    }
                    
                    println!("  {} [label=\"{}\", fillcolor=\"{}\"];", id, label, color);
                } else if tokens[0] == "L" {
                    let from = tokens[1];
                    let from_ori = tokens[2];
                    let to = tokens[3];
                    let to_ori = tokens[4];
                    
                    let mut color = "black";
                    let mut style = "solid";
                    if from_ori == "-" || to_ori == "-" {
                        color = "red";
                        style = "dashed";
                    }
                    println!("  {} -> {} [label=\"{}{}\", color=\"{}\", style=\"{}\"];", from, to, from_ori, to_ori, color, style);
                }
            }
            println!("}}");
        }
    }
    Ok(())
}
