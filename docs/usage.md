# Svaha2 Extended Usage Guide

`svaha2` is designed for high-performance variation graph construction, specifically optimized for cancer genomics workflows where multiple variant types and formats (VCF, MAF, SV, BEDPE) must be integrated into a unified graph.

## Build Command Options

The `build` command is the primary entry point for graph construction.

```bash
./svaha build -r <reference.fa> [options]
```

### Required Options
*   `-r <file>`: The reference genome in FASTA format. A corresponding `.fai` index file must exist in the same directory.

### Variant Input Options
You can provide multiple variant files, even of different formats, in a single run:

*   `-v | --vcf <file>`: Standard VCF 4.0+ file. Handles SNPs, Indels, and SVs (DEL, INS, INV, BND).
*   `--maf <file>`: cBioPortal/TCGA Mutation Annotation Format. Optimized for somatic mutations (SNP, DNP, TNP, ONP, INS, DEL).
*   `--sv <file>`: cBioPortal structural variant format (`data_sv.txt`). Parses connection types (`3to5`, `5to5`, etc.) to set correct graph orientations.
*   `--bedpe <file>`: Standard BEDPE format for structural variants.

### Advanced Construction Options
*   `-m <int>`: **Maximum Node Size** (default: 32). This controls the discretization of the reference backbone. Smaller values increase graph resolution but result in more nodes/edges. Larger values reduce memory but might "hide" variants if they fall within the same node.
*   `-R <chrom:start-end>`: **Region Filtering**. Only variants and reference sequences within this coordinate range will be processed. Coordinates are 1-based.
    *   Example: `-R chr17:7500000-8000000`
*   `--relax-chrom`: **Chromosome Naming Normalization**. By default, `svaha2` requires exact matches between variant files and the reference genome. Use this flag to allow automated mapping between plain numbers (e.g., `17`) and `chr`-prefixed names (e.g., `chr17`).

---

## Detailed Format Specifications

### cBioPortal MAF (`--maf`)
`svaha2` parses the following columns from cBioPortal MAFs:
- `Chromosome`, `Start_Position`, `End_Position`
- `Reference_Allele`, `Tumor_Seq_Allele2`
- `Variant_Type` (SNP, DNP, TNP, ONP, INS, DEL)

Small variants are automatically normalized:
- Leading anchor bases are stripped to ensure clean graph bubbles.
- Coordinates are adjusted to represent only the novel sequence.

### cBioPortal Structural Variants (`--sv`)
Designed for `data_sv.txt` files, this parser interprets:
- `Site1_Chromosome`, `Site1_Position`
- `Site2_Chromosome`, `Site2_Position`
- `Class` (Inversion, Deletion, Translocation, etc.)
- `Connection_Type` (Used to determine head-to-head, tail-to-tail, or standard orientations).

### VCF (`-v`)
The VCF parser handles:
- **Small Variants**: Automatically detects SNPs and Indels.
- **Structural Variants**: Parses `SVTYPE` and `END` tags.
- **Breakends**: Supports `BND` and `TRA` types for translocations.
- **Multi-allelic Sites**: Correctly creates separate alternate paths for each allele at a single position.

---

## Progress and Monitoring

For large genomes (e.g., human hg38), construction can take several minutes. `svaha2` provides structured reporting to `stderr`:

1.  **Phase 0 (Parsing)**: Reports the file being read and the number of variants successfully loaded.
2.  **Discretization**: Informs you when the reference backbone is being calculated.
3.  **Phase 1 (Topology)**: Summarizes the number of backbone nodes, variant nodes, and edges created.
4.  **Phase 2 (Materialization)**: For full-genome builds, reports progress every 1,000,000 nodes.

## Example Workflow

```bash
# 1. Compile
make

# 2. Build graph for a specific cancer study
./svaha build -r /path/to/hg38.fa \
    --maf study_data/data_mutations.txt \
    --sv study_data/data_sv.txt \
    -v common_germline.vcf \
    --relax-chrom \
    -R chr17 > my_graph.gfa

# 3. Check graph stats
./svaha stats my_graph.gfa

# 4. Visualize the graph (regional only)
./svaha view my_graph.gfa > my_graph.dot
dot -Tpdf my_graph.dot -o my_graph.pdf
```

---

## Visualization Tools

1.  **Built-in `view` command**: Best for small regional graphs (e.g., < 100 nodes). Converts GFA to DOT format for use with Graphviz.
2.  **Bandage**: The gold standard for large-scale GFA visualization. It provides an interactive GUI for exploring complex variation graphs.
3.  **vg view**: If you have the `vg` toolkit, you can use `vg view -d` to generate more complex SVG/DOT visualizations.
