svaha2: linear time, low-memory construction of variation graphs
----------------------------------------------------------------
Eric T Dawson  
June 2019  


[![Build Status](https://dev.azure.com/ericco92/ericco92/_apis/build/status/edawson.svaha2?branchName=master)](https://dev.azure.com/ericco92/ericco92/_build/latest?definitionId=1&branchName=master)

## Intro
svaha2 is a reimplementation of the [svaha](https://github.com/edawson/svaha) algorithm for
variation graph construction. It cleans up the code significantly and makes use of better parsing
libraries. Performance is slightly improved.


With svaha2, you can make variation graphs from structural variants:  
- [x] Deletions  
- [x] Inversions  
- [x] Insertions  
- [ ] SNPs  
- [ ] Duplications  
- [x] Translocations   
- [ ] Breakpoints  

In addition, there are currently some limitations on variants:

- While variants may be nested, they cannot share a given breakpoint.  
- Variants may only possess a single alternate allele  
- Variants must possess the following fields:  
  - SVTYPE (DEL, INS, INV, TRA)
  - END (DEL, INV, TRA)  
  - CHR2 (TRA)  
  - SEQ (INS)  

Given a well-formed VCF and a FASTA reference, svaha2 can produce a GFA graph in a matter of minutes for as many as several hundred thousand variants in under 16GB of RAM.

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

Roughly, the svaha algorithm works as follows:
```
1. Read in all variants. Compute the breakpoints for each variant, and store a map(basepair -> variant allele).
2. Add any breakpoints needed for maxmimum node size.
3. For each contig, sort its breakpoints. Remove any duplicates and any invalid ones (first/last base of sequence).
4. For each contig:
    5. For each breakpoint b in the contig's breakpoints  
        - Create a node
        - Fill it with the reference subseqence from the previous breakpoint to the current one.
        - If there is a variant at the node's position, cache the node in a map from position to node.
        - If the node is cached, also cache the preceding node.
        - If the variant is flat, or an insertion / substitution, create a node. Store it in a map from position to inserted nodes.  
        - Output the node in GFA format, add it the ref path if needed, and add an edge to the previous ref node if this node is also on the reference
        path.
6. For each contig  
    7. For each cached position-variant allele in the position->variant allele map. 
        - Collect the relevant nodes before / in / after the variant (depending on if it's flat or not).  
        - Create edges to represent the variant based on its SVTYPE.  
        - Ouput nodes and edges in GFA.  

```
## Feature Roadmap

- [ ] Lightweight storage of GFA paths (probably file backed). These are super expensive at the moment.
- [ ] Output maps to memory-mapped files, rather than keeping them in memory.  
- [ ] Big performance boost by reducing number of loops (all variant types except interchromosomals). This requires file-backed maps.
- [ ] Multiple variants at a site (will increase memory usage).  
- [ ] We build everything by contig right now. That's probably fine, but it could be too much for SNVs and indels.  
- [ ] Duplication handling. This is easy, but it's a bit moot given that aligners don't like cyclic structures.  
- [ ] Tabix / VCF-index handling. This will greatly improve memory usage by no longer requiring alleles to stay in memory.  
- [ ] GFA indexing using [tinyGFA](https://github.com/edawson/tinyGFA). This means we no longer need to cache nodes, just their IDs (128+bit struct -> 64bit Int).  
- [ ] Optional GFA1 / GFA2 output. Right now, we just output GFA1 (but this is easily converted with GFAKluge). All Links are valid Edges so this just requires adding the args and some hhandling logic.  
- [ ] Threading (zoom zoom). But we'll have to be careful about buffering output so the threads don't just thrash the disk and making sure to keep RAM usage down by getting rid of the reference sequences ASAP (although this is still a max of <4GB if we keep the whole Human genome).

