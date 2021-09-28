#load packages
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggthemr))
suppressMessages(library(gridExtra))
suppressMessages(library(RCircos))

#get sample name
vcf.list = list.files(path = "in/SVs/", recursive = TRUE, pattern = "\\.vcf$", full.names = TRUE)
sample_name = gsub(".{4}$", '', vcf.list)
sample_name = substring(sample_name, 9)
now = format(Sys.time(), "%d_%m_%Y")

#set parameters
txtFileName <- paste0(sample_name, "_vcf_man_", now)

#read SV callset into R data frame (skips first Header line)
sv_calls = read.table(file = paste0("out/SVs/tables/", txtFileName, "_structuralvariants.txt"), header = TRUE, sep = "\t")

#set titles, subtitles, axis, scales etc.
#violin
plot.title.violin = "Structural Variant Size Distribution"
y.axis.name.violin = "Size (bp)"

#chrdist
plot.title.chrbox = "Structural Variant Chromosome Distribution"
y.axis.name.chrbox = "Count (n)"

#pie
plot.title.pie = "Structural Variant Class"

#binned box
plot.title.binned = "Binned Structural Variant Size"
y.axis.name.binned = "Count (n)"

#genotype box
plot.title.genotype = "Genotypes"
y.axis.name.genotype = "Count (n)"

#violin
#convert variable to factor
sv_calls$sv_fam = as.factor(sv_calls$sv_fam)

#chrdist
#set max limits for y scale
sv_tmp = sv_calls %>% add_count(chr)
ymax = max(sv_tmp$n)

#subset dataframe on sv subtype
sv_del = filter(sv_calls, sv_fam == "del")
sv_dup = filter(sv_calls, sv_fam == "dup")
sv_ins = filter(sv_calls, sv_fam == "ins")

#add chr count variable
sv_del = sv_del %>% add_count(chr)
sv_dup = sv_dup %>% add_count(chr)
sv_ins = sv_ins %>% add_count(chr)

#remove duplicate chr rows
sv_del = distinct(sv_del, chr, .keep_all = TRUE)
sv_dup = distinct(sv_dup, chr, .keep_all = TRUE)
sv_ins = distinct(sv_ins, chr, .keep_all = TRUE)

#melt dataframes into one
sv_calls.count = rbind(sv_del, sv_dup, sv_ins)

#pie
#subset data for piechart
sv_calls_pie = sv_calls %>%
  group_by(sv_fam) %>%
  summarize(count = n())

#count number of variants
sv_n = sum(sv_calls_pie$count)

#count percentage of called sv classes
sv_calls_pie$percentage = ((sv_calls_pie$count / sv_n) * 100)

#binned box
#set up cut-off values 
breaks <- c(1, 10, 50, 100, 500, 1000, 10000, 100000, 1000000)

#specify interval/bin labels
tags <- c("1-10bp","10bp-50bb", "50bp-100bp", "100bp-500bp", "500bp-1kb", "1kb-10kb", "10kb-100kb", "100kb-1mb")

#bucketing values into bins
sv_tags <- cut(sv_calls$sv_length, breaks=breaks, include.lowest=TRUE, right=FALSE, labels=tags)

#genotype box
#create summary table
sv_tab = sv_calls %>% select(sv_length, sv_type, sv_fam, genotype)

#turn into factor
sv_tab$genotype = as.factor(sv_tab$genotype)

#sub genotype info
levels(sv_tab$genotype)[levels(sv_tab$genotype)=="1|1"] = "hom"
levels(sv_tab$genotype)[levels(sv_tab$genotype)=="1|0"] = "het"
levels(sv_tab$genotype)[levels(sv_tab$genotype)=="0|1"] = "het"

#subset dataframe on sv subtype
sv_tab_del = filter(sv_tab, sv_fam == "del")
sv_tab_dup = filter(sv_tab, sv_fam == "dup")
sv_tab_ins = filter(sv_tab, sv_fam == "ins")

#add chr count variable
sv_tab_del = sv_tab_del %>% add_count(genotype)
sv_tab_dup = sv_tab_dup %>% add_count(genotype)
sv_tab_ins = sv_tab_ins %>% add_count(genotype)

#remove duplicate chr rows
sv_tab_del = distinct(sv_tab_del, genotype, .keep_all = TRUE)
sv_tab_dup = distinct(sv_tab_dup, genotype, .keep_all = TRUE)
sv_tab_ins = distinct(sv_tab_ins, genotype, .keep_all = TRUE)

#melt dataframes into one
sv_tab_genotype_count = rbind(sv_tab_del, sv_tab_dup, sv_tab_ins)

#subset genotype table
sv_tab_genotype_count = sv_tab_genotype_count %>% select(sv_fam, genotype, n)

##plotting
#set plot theme
ggthemr("dust")

#violin plot
sv_size_violine = ggplot(sv_calls, aes(x = sv_fam, y = sv_length, fill = sv_fam)) + 
  labs(title = plot.title.violin, x = "", y = y.axis.name.violin) +
  geom_violin(trim = FALSE, scale = "width", color = NA) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +
  theme(legend.position = "none") +
  scale_y_log10() 

#chromosome distribution box plot
sv_chrdist_box = ggplot(sv_calls.count, aes(x = chr, y = n, fill = sv_fam)) +
  labs(title = plot.title.chrbox, x = "", y = y.axis.name.chrbox, fill = "") +
  scale_x_discrete(limits=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) +
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, ymax, by = 2))

#binned sv size
sv_binned = ggplot(data = as_tibble(sv_tags), mapping = aes(x=value)) + 
  geom_bar(fill = "#48727D") + 
  theme(axis.title.x = element_blank()) +
  labs(title = plot.title.binned, y = y.axis.name.binned, fill = "")

#save sv genotypes as table instead of piechart
sum_genotypes = as.data.frame(sv_tab_genotype_count)
colnames(sum_genotypes) <- c("SV Type", "Genotype", "Count")
sum_genotypes_grob = tableGrob(sum_genotypes, rows = NULL)

#plot header for raport
text = paste0(sample_name, " | GRCh38 | Structural Variants | ", now)
plot.title = ggplot() + 
  annotate("text", x = 7, y = 3, size=8, label = text) + 
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = "white"))

#plot empty box for spacing
box = ggplot() +
  geom_blank(mapping = NULL, data = NULL, stat = "identity", position = "identity", show.legend = NA, inherit.aes = TRUE) + 
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(), panel.background = element_rect(fill = "white", colour = "white"))

#circos plot
#subset dataframe on sv subtype
sv_del_bed = filter(sv_calls, sv_fam == "del") %>% 
  select(chr, start, end)

sv_dup_bed = filter(sv_calls, sv_fam == "dup") %>% 
  select(chr, start, end)

sv_ins_bed = filter(sv_calls, sv_fam == "ins") %>% 
  select(chr, start, end)

#define reference build
data(UCSC.HG38.Human.CytoBandIdeogram)

#set core components
RCircos.Set.Core.Components(cyto.info = UCSC.HG38.Human.CytoBandIdeogram, chr.exclude = NULL, tracks.inside = 3, tracks.outside = 0)

#set plor parameters
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
RCircos.Reset.Plot.Parameters(rcircos.params)
RCircos.Tile.Plot(sv_del_bed, track.num = 1, side = "in")

#dup
rcircos.params$track.background = "sandybrown"
RCircos.Reset.Plot.Parameters(rcircos.params)
RCircos.Tile.Plot(sv_dup_bed, track.num = 2, side = "in")

#ins
rcircos.params$track.background = "wheat4"
RCircos.Reset.Plot.Parameters(rcircos.params)
RCircos.Tile.Plot(sv_ins_bed, track.num = 3, side = "in")

# finnish
whatever <- dev.off()

#export plots
ggsave(sv_size_violine, file = paste0("out/SVs/figs/", txtFileName, "_02_sv_size_violin.png"), limitsize = FALSE, width = 7, height = 7, units = c("in"), dpi = 300)
ggsave(sv_chrdist_box, file = paste0("out/SVs/figs/", txtFileName, "_01_sv_chrdist.png"), limitsize = FALSE, width = 14, height = 6, units = c("in"), dpi = 300)
ggsave(sv_binned, file = paste0("out/SVs/figs/", txtFileName, "_03_sv_binned.png"), limitsize = FALSE, width = 7, height = 7, units = c("in"), dpi = 300)
ggsave(plot.title, file = paste0("out/SVs/figs/", txtFileName, "_header.png"), limitsize = FALSE, width = 14, height = 1, units = c("in"), dpi = 300)
ggsave(sum_genotypes_grob, file = paste0("out/SVs/figs/", txtFileName, "04_sum_genotypes.png"), limitsize = FALSE, width = 3, height = 2, units = c("in"), dpi = 300)
ggsave(box, file = paste0("out/SVs/figs/", txtFileName, "_box1.png"), limitsize = FALSE, width = 2, height = 2, units = c("in"), dpi = 300)
ggsave(box, file = paste0("out/SVs/figs/", txtFileName, "_box2.png"), limitsize = FALSE, width = 14, height = 0.3, units = c("in"), dpi = 300)

#prompt message to terminal
cat ("Figures and summaries exported...\n")
