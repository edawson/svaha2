# svaha2: high-performance variation graph construction

`svaha2` is a high-performance tool for constructing variation graphs from a reference genome and multiple variant input formats. It uses a **breakpoint-based discretization** approach to divide the reference into non-overlapping nodes and overlays variation efficiently.

## Features
- **Comprehensive Variant Support**:
  - Small variants: SNPs, MNPs (DNP, TNP, ONP), and Indels.
  - Structural Variants (SVs): Deletions, Insertions, Inversions, Translocations, and Duplications.
- **Multi-format Input**:
  - Standard VCF files.
  - cBioPortal Mutation Annotation Format (MAF).
  - cBioPortal Structural Variant files (`data_sv.txt`).
  - BEDPE structural variants.
- **High Performance**:
  - **Zero-copy loading**: Uses memory-mapped files (`mmap`) for fast reference sequence access.
  - **Succinct storage**: 2-bit sequence packing for alternate alleles to minimize RAM usage.
  - **Parallelized Materialization**: Multi-threaded GFA generation for extreme throughput.
  - **rGFA Compatible**: Outputs GFA 1.1 with stable coordinate tags (`SN`, `SO`, `SR`).

## Quick Start

### 1. Compilation
Requires a C++17 compatible compiler (e.g., `clang++` or `g++`).
```bash
make
```

### 2. Basic Usage
Construct a GFA graph from a reference and variants:
```bash
./svaha build -r reference.fa -v variants.vcf > output.gfa
```

### 3. Cancer Genomics (cBioPortal/TCGA)
Load multiple files from different sources simultaneously:
```bash
./svaha build -r hg38.fa \
    --maf data_mutations.txt \
    --sv data_sv.txt \
    -v background_variants.vcf > cancer_graph.gfa
```

### 4. Docker
You can also run `svaha` using Docker:
```bash
# Build the image
docker build -t svaha .

# Run the image
docker run --rm -v $(pwd):/data svaha build -r /data/ref.fa -v /data/vars.vcf > output.gfa
```

## Command Line Interface

For detailed information on all input formats and construction options, see the [Extended Usage Guide](docs/usage.md).

### `build`
Construct a GFA from reference and variants.
- `-r <file>`: Reference FASTA (required).
- `-v <file>`: VCF variants.
- `--maf <file>`: Mutation Annotation Format (cBioPortal/TCGA).
- `--sv <file>`: Structural Variants (cBioPortal `data_sv.txt`).
- `--bedpe <file>`: BEDPE structural variants.
- `-m <int>`: Maximum node size (default: 32 bp).
- `-R <region>`: Restrict to a specific region (`chrom:start-end`).
- `--relax-chrom`: Allow matching between 'chr1' and '1' (useful for hg19/hg38 mismatches).

### `stats`
Show summary statistics for a GFA file.
```bash
./svaha stats output.gfa
```

## Algorithm
1. **Breakpoint identification**: Extract genomic positions where variants start or end.
2. **Backbone discretization**: Divide the reference into stable segments based on these breakpoints.
3. **Graph Materialization**: 
   - Emit reference nodes and edges.
   - Overlay variant nodes and oriented edges based on variant type and strand/orientation information.

## License
MIT

