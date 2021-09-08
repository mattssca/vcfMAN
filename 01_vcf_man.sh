#!/bin/bash

read -p "Is input VCF including small_variants, SVs or both? " VCF
if [ "$VCF" = "small_variants" ]; then
  Rscript scripts/read_vcf_smallvariants.R
elif [ "$VCF" = "SVs" ]; then
  Rscript scripts/read_vcf_structuralvariants.R
elif [ "$VCF" = "both" ]; then
  Rscript scripts/read_vcf_structuralvariants.R 
  Rscript scripts/read_vcf_smallvariants.R
else
  echo "Please enter a valid input"
  sh 01_vcf_man.sh 
fi