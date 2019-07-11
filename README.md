svaha2: linear time, low-memory construction of variation graphs
----------------------------------------------------------------
Eric T Dawson  
June 2019  


[![Build Status](https://dev.azure.com/ericco92/ericco92/_apis/build/status/edawson.svaha2?branchName=master)](https://dev.azure.com/ericco92/ericco92/_build/latest?definitionId=1&branchName=master)

## Intro
svaha2 is a reimplementation of the [svaha](https://github.com/edawson/svaha) algorithm for
variation graph construction. It cleans up the code significantly and makes use of better parsing
libraries. Performance is slightly improved.

## Usage
```
svaha2: linear-time, low-memory construction of variation graphs.
Usage: svaha2 [options] -r <ref.fa> -v <variants.vcf>
options:
-m / --max-node-size : maximum length (in basepairs) of a node in the graph.
-r / --reference     : The reference genome to use for construction.
-v / --vcf           : A VCF file to use for construction.
-b / --bed           : a bed file to restrict chromosomes to.
-T / --no-translocations  : ignore interchromosomal variants.
-f / --flat          : use flat alternates (every allele is represented by at least one node).
-p / --paths         : output path information for variants.
-I / --insertions         : FASTA file of insertion variant sequences.
```
## Algorithm

## Testing

## Getting Help

