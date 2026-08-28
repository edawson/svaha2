# Tutorial: Building Variation Graphs with cBioPortal Data

This tutorial walks you through using `svaha2` to construct a variation graph using cancer genomics data from cBioPortal. We will use real hg38 data from the Lung Adenocarcinoma (TCGA, PanCancer Atlas) study.

## Prerequisites

1.  **Compiled `svaha2`**: Run `make` in the repository root.
2.  **Reference Genome**: A human reference genome (GRCh38/hg38). You can download a standard version here:
    *   [Homo_sapiens_assembly38.fasta](https://public.betulalabs.com/References/Homo_sapiens_assembly38.fasta)
    *   *Note: Ensure you also have the `.fai` index file.*
3.  **Study Data**: We will use the `luad_tcga_gdc` dataset found in `ext-data/`.

---

## Step 1: Verify the Reference and Data

Ensure your reference genome is indexed:
```bash
# If you don't have the index, create it with samtools:
# samtools faidx Homo_sapiens_assembly38.fasta
ls Homo_sapiens_assembly38.fasta.fai
```

Check the LUAD mutations file:
```bash
head -n 5 ext-data/luad_tcga_gdc/data_mutations.txt | cut -f 1-6
```
You will see that `NCBI_Build` is `GRCh38` and the `Chromosome` column uses plain numbers (e.g., `1`).

---

## Step 2: Build a Regional Graph (hg38)

We will build a graph for a 50kb region on Chromosome 1 containing the `CHD5` gene.

```bash
./svaha build \
    -r Homo_sapiens_assembly38.fasta \
    --maf ext-data/luad_tcga_gdc/data_mutations.txt \
    --relax-chrom \
    -R chr1:6100000-6150000 > luad_chd5.gfa
```

### Key Parameters:
- `-r`: Points to your local copy of the hg38 reference.
- `--maf`: Loads the somatic mutations from the LUAD study.
- `--relax-chrom`: **Critical** for this dataset, as it maps the plain `1` in the MAF to `chr1` in the assembly38 reference.
- `-R`: Focuses the build on the specific `CHD5` genomic window.

---

## Step 3: Analyze the Graph

Once the build completes, verify the results:

```bash
./svaha stats luad_chd5.gfa
```

Expected output for this region:
```text
GFA Statistics for: luad_chd5.gfa
Nodes: 1669
Edges: 1707
Total Sequence Length: 50039 bp
```
This shows that 39 alternate alleles (variants) were successfully integrated into the reference backbone.

---

## Step 4: Full Genome Construction (Optional)

To build a comprehensive graph for the entire LUAD cohort across all chromosomes:

```bash
./svaha build \
    -r Homo_sapiens_assembly38.fasta \
    --maf ext-data/luad_tcga_gdc/data_mutations.txt \
    --relax-chrom > luad_full_genome.gfa
```
*Note: This will process approximately 280MB of mutation data and materialize the full hg38 backbone. Ensure you have sufficient disk space for the output GFA.*

---

## Troubleshooting

### Mismatched Builds
If your data is hg19 (GRCh37) but your reference is hg38, the variants will not map to the correct physical locations. Always check the `NCBI_Build` column in the MAF header.

### Empty Graphs
If `svaha stats` shows 0 alternate nodes, double-check that:
1.  The chromosome names match (or `--relax-chrom` is used).
2.  Your `-R` region uses the same coordinate system (hg38) as the reference genome.
