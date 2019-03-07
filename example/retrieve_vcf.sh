#!/usr/bin/env bash

# downloads a file containing VCF data from UVA box
# reference allele coded as 1, alternate allele coded as 2
# missing data coded as 0
wget -O haplotypes.vcf.gz -v -L https://virginia.box.com/shared/static/8b59b50bf29k03om2sjynkmd0s7pl2ws.gz


# downloadsa  file containing haplotype paths from UVA box
wget -O haplotypes.txt.gz -v -L https://virginia.box.com/shared/static/iyf5303jtn0r8bqh5a9xqqfdjq48lbzt.gz
