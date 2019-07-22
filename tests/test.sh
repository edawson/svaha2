#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap


plan tests 3

svaha=../svaha2
data=..//data


## Deletion tests
is "$(${svaha} -r ${data}/tgraph.fa -v ${data}/tdel.vcf  | grep -v "^P\|^H" | sort | md5sum | cut -f 1 -d ' ')" "239b9ae8a090f7e79c8f9c65f88a9ceb" "Graphical deletions are properly produced"
is "$(${svaha} -r ${data}/tgraph.fa -v ${data}/tdel.vcf -f | grep -v "^P\|^H" | sort | md5sum | cut -f 1 -d ' ')" "717cedc93001238cfbc64dd937ad74a7" "Flat deletions are properly produced"

## Inversion Tests

## Insertion Tests

## Translocation tests
is "$(${svaha} -r ${data}/translocation_graph.fa -v ${data}/ttra.vcf -T | sort | md5sum | cut -f 1 -d ' ' )" "7782740039ce19394c4ccd92d34da8e0" "Translocations are skipped with -T"
