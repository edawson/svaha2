#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap


plan tests 2

svaha=../svaha2
data=..//data


is "$(${svaha} -r ${data}/tgraph.fa -v ${data}/tdel.vcf  | grep -v "^P\|^H" | sort | md5sum | cut -f 1 -d ' ')" "239b9ae8a090f7e79c8f9c65f88a9ceb"
is "$(${svaha} -r ${data}/tgraph.fa -v ${data}/tdel.vcf -f | grep -v "^P\|^H" | sort | md5sum | cut -f 1 -d ' ')" "717cedc93001238cfbc64dd937ad74a7"
