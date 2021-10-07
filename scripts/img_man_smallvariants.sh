#!/bin/bash

#combine png for report format
convert out/small_variants/figs/*_box1.png out/small_variants/figs/*_smallvariants_summary.png out/small_variants/figs/*_sum_genotypes.png out/small_variants/figs/*_box1.png +append -quiet out/small_variants/out1.png
convert out/small_variants/figs/*_header.png out/small_variants/figs/*_box2.png out/small_variants/out1.png out/small_variants/figs/*_box2.png -append -quiet out/small_variants/out1.png
convert out/small_variants/figs/*_01_smallvariants_size_violin.png out/small_variants/figs/*_03_snv_chr_dist.png +append -quiet out/small_variants/out2.png

# crop images
convert out/small_variants/figs/ideograms/chr1.png -crop 4424x193+132+415 out/small_variants/figs/ideograms/chr1.png
convert out/small_variants/figs/ideograms/chr2.png -crop 4424x193+132+415 out/small_variants/figs/ideograms/chr2.png
convert out/small_variants/figs/ideograms/chr3.png -crop 4424x193+132+415 out/small_variants/figs/ideograms/chr3.png
convert out/small_variants/figs/ideograms/chr4.png -crop 4424x193+132+415 out/small_variants/figs/ideograms/chr4.png
convert out/small_variants/figs/ideograms/chr5.png -crop 4424x193+132+415 out/small_variants/figs/ideograms/chr5.png
convert out/small_variants/figs/ideograms/chr6.png -crop 4424x193+132+415 out/small_variants/figs/ideograms/chr6.png
convert out/small_variants/figs/ideograms/chr7.png -crop 4424x193+132+415 out/small_variants/figs/ideograms/chr7.png
convert out/small_variants/figs/ideograms/chr8.png -crop 4424x193+132+415 out/small_variants/figs/ideograms/chr8.png
convert out/small_variants/figs/ideograms/chr9.png -crop 4424x193+132+415 out/small_variants/figs/ideograms/chr9.png
convert out/small_variants/figs/ideograms/chr10.png -crop 4424x193+132+415 out/small_variants/figs/ideograms/chr10.png
convert out/small_variants/figs/ideograms/chr11.png -crop 4424x193+132+415 out/small_variants/figs/ideograms/chr11.png
convert out/small_variants/figs/ideograms/chr12.png -crop 4424x193+132+415 out/small_variants/figs/ideograms/chr12.png
convert out/small_variants/figs/ideograms/chr13.png -crop 3036x193+132+415 out/small_variants/figs/ideograms/chr13.png
convert out/small_variants/figs/ideograms/chr14.png -crop 3036x193+132+415 out/small_variants/figs/ideograms/chr14.png
convert out/small_variants/figs/ideograms/chr15.png -crop 3036x193+132+415 out/small_variants/figs/ideograms/chr15.png
convert out/small_variants/figs/ideograms/chr16.png -crop 3036x193+132+415 out/small_variants/figs/ideograms/chr16.png
convert out/small_variants/figs/ideograms/chr17.png -crop 3036x193+132+415 out/small_variants/figs/ideograms/chr17.png
convert out/small_variants/figs/ideograms/chr18.png -crop 3036x193+132+415 out/small_variants/figs/ideograms/chr18.png
convert out/small_variants/figs/ideograms/chr19.png -crop 3036x193+132+415 out/small_variants/figs/ideograms/chr19.png
convert out/small_variants/figs/ideograms/chr20.png -crop 3036x193+132+415 out/small_variants/figs/ideograms/chr20.png
convert out/small_variants/figs/ideograms/chr21.png -crop 3036x193+132+415 out/small_variants/figs/ideograms/chr21.png
convert out/small_variants/figs/ideograms/chr22.png -crop 3036x193+132+415 out/small_variants/figs/ideograms/chr22.png
convert out/small_variants/figs/ideograms/plot.title.png -crop 9448x700+0+221 -quiet out/small_variants/figs/ideograms/plot.title.png 

# compile ideograms fro chr 1-12
convert out/small_variants/figs/ideograms/chr1.png out/small_variants/figs/ideograms/chr2.png out/small_variants/figs/ideograms/chr3.png out/small_variants/figs/ideograms/chr4.png out/small_variants/figs/ideograms/chr5.png out/small_variants/figs/ideograms/chr6.png out/small_variants/figs/ideograms/chr7.png out/small_variants/figs/ideograms/chr8.png out/small_variants/figs/ideograms/chr9.png out/small_variants/figs/ideograms/chr10.png out/small_variants/figs/ideograms/chr11.png out/small_variants/figs/ideograms/chr12.png out/small_variants/figs/*_box2.png -append -quiet out/small_variants/figs/ideograms/1-12.ideogram.png

# compile ideograms for chr 13-x(y)
convert out/small_variants/figs/ideograms/chr13.png out/small_variants/figs/ideograms/chr14.png out/small_variants/figs/ideograms/chr15.png out/small_variants/figs/ideograms/chr16.png out/small_variants/figs/ideograms/chr17.png out/small_variants/figs/ideograms/chr18.png out/small_variants/figs/ideograms/chr19.png out/small_variants/figs/ideograms/chr20.png out/small_variants/figs/ideograms/chr21.png out/small_variants/figs/ideograms/chr22.png out/small_variants/figs/*_box2.png -append -quiet out/small_variants/figs/ideograms/13-22.ideogram.png

# add both ideograms horizonatally
convert +append -quiet out/small_variants/figs/ideograms/1-12.ideogram.png out/small_variants/figs/ideograms/13-22.ideogram.png out/small_variants/figs/ideograms/ideogram.png

# add plot title
convert out/small_variants/figs/ideograms/plot.title.png out/small_variants/figs/ideograms/ideogram.png	-append -quiet out/small_variants/figs/ideograms/ideogram.png

#resize ideogram to fit report format
convert out/small_variants/figs/ideograms/ideogram.png -crop 7000x5106+0+0 out/small_variants/figs/ideograms/ideogram.png
convert -resize 4200 -quality 10 out/small_variants/figs/ideograms/ideogram.png out/small_variants/figs/ideograms/ideogram.png

#compile report
convert out/small_variants/out1.png out/small_variants/figs/*_box2.png out/small_variants/out2.png out/small_variants/figs/*_box2.png out/small_variants/figs/*_box2.png out/small_variants/figs/*_02_snv_distance.png out/small_variants/figs/ideograms/ideogram.png out/small_variants/figs/*_box2.png -append out/small_variants/report.png

#delete unused files
rm out/small_variants/out1.png
rm out/small_variants/out2.png
rm out/small_variants/figs/*_header.png
rm out/small_variants/figs/*_box2.png
rm out/small_variants/figs/ideograms/plot.title.png
rm out/small_variants/figs/ideograms/1-12.ideogram.png
rm out/small_variants/figs/ideograms/13-22.ideogram.png
rm out/small_variants/figs/ideograms/chr1.png
rm out/small_variants/figs/ideograms/chr2.png
rm out/small_variants/figs/ideograms/chr3.png
rm out/small_variants/figs/ideograms/chr4.png
rm out/small_variants/figs/ideograms/chr5.png
rm out/small_variants/figs/ideograms/chr6.png
rm out/small_variants/figs/ideograms/chr7.png
rm out/small_variants/figs/ideograms/chr8.png
rm out/small_variants/figs/ideograms/chr9.png
rm out/small_variants/figs/ideograms/chr10.png
rm out/small_variants/figs/ideograms/chr11.png
rm out/small_variants/figs/ideograms/chr12.png
rm out/small_variants/figs/ideograms/chr13.png
rm out/small_variants/figs/ideograms/chr14.png
rm out/small_variants/figs/ideograms/chr15.png
rm out/small_variants/figs/ideograms/chr16.png
rm out/small_variants/figs/ideograms/chr17.png
rm out/small_variants/figs/ideograms/chr18.png
rm out/small_variants/figs/ideograms/chr19.png
rm out/small_variants/figs/ideograms/chr20.png
rm out/small_variants/figs/ideograms/chr21.png
rm out/small_variants/figs/ideograms/chr22.png

echo "Summary report for small variants compiled...\n"
