#!/bin/bash

BED="$1"

if [ "$BED" = "small_variants" ]; then
Rscript scripts/read_bed_dip_variants.R
convert -quiet out/small_variants/report.png out/BED/bed_chrdist_frac.png -append -quiet out/small_variants/report.png

elif [ "$BED" = "SVs" ]; then
Rscript scripts/read_bed_dip_variants.R
convert -quiet out/SVs/report.png out/BED/bed_chrdist_frac.png -append -quiet out/SVs/report.png

elif [ "$BED" = "both" ]; then
Rscript scripts/read_bed_dip_variants.R
convert -quiet out/report.png out/BED/bed_chrdist_frac.png -append -quiet out/report.png

else
  echo "Please enter a valid input! Specify the vcf input from vcf_man.sh (acceptable commands are; small_variants, SVs and both)\n"
fi

echo "BED coverage plot added to report...\n"