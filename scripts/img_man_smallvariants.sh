#!/bin/bash

#combine png for report format
convert out/SVs/figs/*_box1.png out/SVs/figs/*_SV_summary.png out/SVs/figs/*04_sum_genotypes.png out/SVs/figs/*_box1.png +append -quiet out/SVs/out1.png
convert out/SVs/figs/*_header.png out/SVs/out1.png out/SVs/figs/*_box2.png -append -quiet out/SVs/out1.png
convert -quiet out/SVs/figs/*02_sv_size_violin.png out/SVs/figs/*03_sv_binned.png +append out/SVs/out2.png
convert -density 300x300 -units pixelsperinch out/SVs/figs/*_circos.pdf -background white -alpha background -alpha off -compress zip +adjoin out/SVs/figs/circos.png
convert -quiet out/SVs/out1.png out/SVs/figs/*01_sv_chrdist.png out/SVs/out2.png out/SVs/figs/circos.png -append out/SVs/report.png

#delete unused files
rm out/SVs/out1.png
rm out/SVs/out2.png
rm out/SVs/figs/*_box1.png
rm out/SVs/figs/*_box2.png
rm out/SVs/figs/*_header.png
rm out/SVs/figs/circos.png

echo "Summary report compiled...\n"
echo "DONE!\n"