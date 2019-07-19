#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap


plan tests 2

svaha=../svaha2
data=..//data

is "$(${svaha} -r ${data}/tgraph.fa -v ${data}/tdel.vcf  | sort | md5sum | cut -f 1 -d ' ')" "6bf0ee7e84500dc135a91b8db14289cd"
is "$(${svaha} -r ${data}/tgraph.fa -v ${data}/tdel.vcf -f | sort | md5sum | cut -f 1 -d ' ')" "0eb09ff47c49a4332ce79d187121ab68"
