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

We will build a graph for a 100bp region on Chromosome 17 containing a `TP53` mutational hotspot. To show how SVs and SNVs integrate, we've provided a sample SV file in `docs/example_sv.txt`.

```bash
./svaha build \
    -r Homo_sapiens_assembly38.fasta \
    --maf ext-data/luad_tcga_gdc/data_mutations.txt \
    --sv docs/example_sv.txt \
    --relax-chrom \
    -R chr17:7675500-7676500 -m 32 > tp53_hotspot.gfa
```

### Key Parameters:
- `-r`: Points to your local copy of the hg38 reference.
- `--maf`: Loads somatic mutations from the LUAD study.
- `--sv`: Loads a sample inversion (both junctions) in the same region.
- `--relax-chrom`: Maps plain `17` in the data to `chr17` in the assembly38 reference.
- `-R`: Focuses on a 1kb window, providing enough **padding** to see the nodes connected by the structural variant.
- `-m 32`: Sets a small node size for higher resolution in the visualization.

---

## Step 3: Visualize the Graph

For small regions, you can generate a Graphviz DOT file to see the graph topology:

```bash
# 1. Convert GFA to DOT
./svaha view tp53_hotspot.gfa > tp53_hotspot.dot

# 2. Render to an image (requires Graphviz installed)
dot -Tpng tp53_hotspot.dot -o tp53_hotspot.png
```

### Interpreting the Visualization:
- **White Nodes**: Reference backbone segments.
- **Blue Nodes**: Alternate alleles (SNPs, insertions) from the MAF.
- **Black Solid Edges**: Standard linear reference flow.
- **Red Dashed Edges**: Non-standard connections (Inversions). You will see these skipping or reversing the flow of backbone nodes.

---

## Step 4: Prototyping a Visualization

If you have Graphviz installed, try this:
```bash
# Extract the backbone and variant nodes
grep "label" tp53_hotspot.dot | head -n 10
```

You will see the discretized reference nodes and the variant nodes (highlighted in blue in the rendered image). This 100bp region is dense enough to show the power of variation graphs in cancer genomics while being small enough to fit on a single screen.

---

## Step 5: Full Genome Construction (Optional)

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
