#!/bin/bash

#combine png for report format
convert out/small_variants/figs/*_smallvariants_summary.png out/small_variants/figs/*_smallvariants_genotypes.png +append -quiet out/small_variants/out1.png
convert out/small_variants/figs/*_smallvariants_large_variants.png out/small_variants/figs/*_smallvariants_nongen_variants.png +append -quiet out/small_variants/out2.png
convert out/small_variants/figs/*_01_smallvariants_size_violin.png out/small_variants/figs/*_03_snv_chr_dist.png +append -quiet out/small_variants/out3.png
convert out/small_variants/figs/*_header.png out/small_variants/out1.png out/small_variants/out2.png out/small_variants/out3.png out/small_variants/figs/*_02_snv_distance.png -append out/small_variants/report.png

#delete unused files
rm out/small_variants/out1.png
rm out/small_variants/out2.png
rm out/small_variants/out3.png
rm out/small_variants/figs/*_header.png

echo "Summary report compiled...\n"
echo "DONE!\n"