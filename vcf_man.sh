#!/bin/bash

vcf="$1"

if [ "$vcf" = "small_variants" ]; then
  sh scripts/unpacksmallvariants.sh
  Rscript scripts/read_vcf_smallvariants.R
  Rscript scripts/plot_smallvariants.R
  sh scripts/img_man_smallvariants.sh
elif [ "$vcf" = "SVs" ]; then
  sh scripts/unpackSVs.sh
  Rscript scripts/read_vcf_structuralvariants.R
  Rscript scripts/plot_structuralvariants.R
  sh scripts/bedtools.sh
  Rscript scripts/Rcircos.R
  sh scripts/img_man_structuralvariants.sh
elif [ "$vcf" = "both" ]; then
  sh scripts/unpacksmallvariants.sh
  sh scripts/unpackSVs.sh
  Rscript scripts/read_vcf_smallvariants.R
  Rscript scripts/read_vcf_structuralvariants.R
  Rscript scripts/plot_smallvariants.R
  Rscript scripts/plot_structuralvariants.R
  sh scripts/bedtools.sh
  Rscript scripts/Rcircos.R
  sh scripts/img_man_smallvariants.sh
  sh scripts/img_man_structuralvariants.sh
  sh scripts/img_man_combine.sh
else
  echo "Please enter a valid input\n"
  echo "Usage: vcf_man.sh {small_variants|SVs|both}\n"
  echo "Example for small variants: vcf_man.sh small_variants\n"
fi
