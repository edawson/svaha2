# Tutorial: Building Variation Graphs with cBioPortal Data

This tutorial walks you through using `svaha2` to construct a variation graph using cancer genomics data from cBioPortal. We will use a subset of the Pancreatic Adenocarcinoma (MSK, 2024) study as our example.

## Prerequisites

1.  **Compiled `svaha2`**: Run `make` in the repository root.
2.  **Reference Genome**: You need a human reference genome (e.g., hg19/GRCh37 or hg38/GRCh38) in FASTA format with a `.fai` index.
3.  **Study Data**: Download a study from [cBioPortal](https://www.cbioportal.org/datasets). You will specifically need:
    *   `data_mutations.txt` (MAF format for SNPs/Indels)
    *   `data_sv.txt` (Structural variants)

---

## Step 1: Identify your Genome Build

cBioPortal studies often use different reference builds. Before building the graph, check the `NCBI_Build` column in your data:

```bash
head -n 5 data_mutations.txt | cut -f 4
```

*   If the study uses **GRCh37** (like `pdac_msk_2024`), you should ideally use an hg19 reference.
*   If you use an hg38 reference with hg19 coordinates, you may need to liftover your data first, or use the `--relax-chrom` flag if the only difference is the `chr` prefix.

---

## Step 2: Build a Regional Graph

Building a full-genome graph takes time and memory. It is often best to start with a region of interest, such as the `TP53` locus.

Assuming you are using an hg38 reference but hg19-based data (where `TP53` is on chromosome 17):

```bash
./svaha build \
    -r path/to/hg38.fasta \
    --maf data_mutations.txt \
    --sv data_sv.txt \
    --relax-chrom \
    -R chr17:7500000-8000000 > tp53_cancer_graph.gfa
```

### What this command does:
- `-r`: Loads the reference genome.
- `--maf`: Parses somatic mutations (SNPs, small indels).
- `--sv`: Parses structural rearrangements (inversions, translocations).
- `--relax-chrom`: Automatically maps `17` in the data to `chr17` in the reference.
- `-R`: Restricts the graph to a 500kb window around `TP53`.

---

## Step 3: Verify the Construction

Check the summary statistics of your new graph:

```bash
./svaha stats tp53_cancer_graph.gfa
```

You should see:
- **Backbone Nodes**: The reference sequence divided by variant breakpoints.
- **Alternate Nodes**: Novel sequences (insertions, SNPs) introduced by the MAF.
- **Edges**: Connections representing deletions, inversions, or standard reference flow.

---

## Step 4: Understanding the GFA Output

Open the `.gfa` file to see how the variants are represented:

- **SNPs**: Look for `SR:i:1` tags on `S` (Segment) lines. These are variant nodes.
- **Deletions**: Look for `L` (Link) lines that skip over reference nodes.
- **Structural Variants**: `svaha2` uses the connection information in `data_sv.txt` to correctly orient edges (e.g., head-to-head links for certain inversions).

---

## Troubleshooting

### "Error: Chromosome '17' not found in reference"
This happens when your reference uses `chr17` but your data uses `17`. 
**Fix**: Add the `--relax-chrom` flag.

### "Done. Loaded 0 variants from this file."
This usually means your `-R` region coordinates do not overlap with any variants in the file, or there is a coordinate system mismatch (e.g., trying to use hg38 coordinates on an hg19 dataset).

---

## Next Steps
Once you are comfortable with regional graphs, you can build a study-wide graph by removing the `-R` flag. This is useful for analyzing the collective mutational landscape of a patient cohort.
