suppressMessages(library(RCircos))

#get sample name
vcf.list = list.files(path = "in/SVs/", recursive = TRUE, pattern = "\\.vcf$", full.names = TRUE)
sample_name = gsub(".{4}$", '', vcf.list)
sample_name = substring(sample_name, 9)
now = format(Sys.time(), "%d_%m_%Y")

#set parameters
txtFileName <- paste0(sample_name, "_vcf_man_", now)

#red merged bed files back into R
sv_del_bed_merged = read.table("out/SVs/tables/bed/sv_del_merged.bed", sep = "\t", header = F)
sv_dup_bed_merged = read.table("out/SVs/tables/bed/sv_dup_merged.bed", sep = "\t", header = F)
sv_ins_bed_merged = read.table("out/SVs/tables/bed/sv_ins_merged.bed", sep = "\t", header = F)

#define reference build
data(UCSC.HG38.Human.CytoBandIdeogram)

#set core components
RCircos.Set.Core.Components(cyto.info = UCSC.HG38.Human.CytoBandIdeogram, chr.exclude = c("chrX", "chrY"), tracks.inside = 3, tracks.outside = 0)

#set plot parameters
rcircos.params = RCircos.Get.Plot.Parameters()

#define plotting parameters
out.file = paste0("out/SVs/figs/", txtFileName, "_circos.pdf")
pdf(out.file, height = 14, width = 14)
par(mai=c(0.25, 0.25, 0.25, 0.25));
plot.new();
plot.window(c(-1.3, 1.3), c(-1.3, 1.3));

#add tracks to plot
RCircos.Chromosome.Ideogram.Plot()

#del
rcircos.params$track.background = "coral2"
rcircos.params$max.layers = 1
RCircos.Reset.Plot.Parameters(rcircos.params)
RCircos.Tile.Plot(sv_del_bed_merged, track.num = 1, side = "in")

#dup
rcircos.params$track.background = "sandybrown"
rcircos.params$max.layers = 1
RCircos.Reset.Plot.Parameters(rcircos.params)
RCircos.Tile.Plot(sv_dup_bed_merged, track.num = 2, side = "in")

#ins
rcircos.params$track.background = "wheat4"
rcircos.params$max.layers = 1
RCircos.Reset.Plot.Parameters(rcircos.params)
RCircos.Tile.Plot(sv_ins_bed_merged, track.num = 3, side = "in")

#finish
whatever <- dev.off()

#prompt message to terminal
cat ("Circos plot (SVs) generated and exported...\n")