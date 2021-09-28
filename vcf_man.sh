#!/bin/bash

read -p "Is input VCF including small_variants, SVs or both? " VCF
if [ "$VCF" = "small_variants" ]; then
  sh scripts/unpacksmallvariants.sh
  Rscript scripts/read_vcf_smallvariants.R
  Rscript scripts/plot_smallvariants.R
  sh scripts/img_man_smallvariants.sh
elif [ "$VCF" = "SVs" ]; then
  sh scripts/unpackSVs.sh
  Rscript scripts/read_vcf_structuralvariants.R
  Rscript scripts/plot_structuralvariants.R
  sh scripts/img_man_structuralvariants.sh
elif [ "$VCF" = "both" ]; then
  sh scripts/unpacksmallvariants.sh
  sh scripts/unpackSVs.sh
  Rscript scripts/read_vcf_smallvariants.R
  Rscript scripts/read_vcf_structuralvariants.R
  Rscript scripts/plot_smallvariants.R
  Rscript scripts/plot_structuralvariants.R
  sh scripts/img_man_smallvariants.sh
  sh scripts/img_man_structuralvariants.sh
  sh scripts/img_man_combine.sh
else
  echo "Please enter a valid input"
  sh 01_vcf_man.sh 
fi