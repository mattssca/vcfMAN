#remove warning messages
options(warn=-1)

#load packages
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggthemr))
suppressMessages(library(gridExtra))
suppressMessages(library(RCircos))
suppressMessages(library(ggthemr))
suppressMessages(library(data.table))
suppressMessages(library(psych))

#get sample name
vcf.list = list.files(path = "in/small_variants/", recursive = TRUE, pattern = "\\.vcf$", full.names = TRUE)
sample_name = gsub(".{4}$", '', vcf.list)
sample_name = substring(sample_name, 20)
now = format(Sys.time(), "%d_%m_%Y")

#set parameters
txtFileName <- paste0(sample_name, "_vcf_man_", now)

#read SV callset into R data frame (skips first Header line)
sv_calls = read.table(file = paste0("out/small_variants/tables/", txtFileName, "_smallvariants.txt"), header = TRUE, sep = "\t")

#set ggplot theme
ggthemr("dust")

#SNV distance dist
#subset all snv to new df
snv_calls = filter(sv_calls, sv_type == "SNV") %>%
  select(chr, start)

#set chr as factor
snv_calls$chr = as.factor(snv_calls$chr)

#subset every chr, remove first line, duplicate last row, merge df, calculate snv distance (bp), remove lst row
snv_chr1 = filter(snv_calls, chr =="chr1")
snv_chr1_tmp = snv_chr1[-1,]
snv_chr1_tmp =  snv_chr1_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr1$next_snv = snv_chr1_tmp$start
snv_chr1$snv_distance = snv_chr1$next_snv - snv_chr1$start
snv_chr1 = snv_chr1 %>% filter(row_number() <= n()-1)

snv_chr2 = filter(snv_calls, chr =="chr2")
snv_chr2_tmp = snv_chr2[-1,]
snv_chr2_tmp =  snv_chr2_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr2$next_snv = snv_chr2_tmp$start
snv_chr2$snv_distance = snv_chr2$next_snv - snv_chr2$start
snv_chr2 = snv_chr2 %>% filter(row_number() <= n()-1)

snv_chr3 = filter(snv_calls, chr =="chr3")
snv_chr3_tmp = snv_chr3[-1,]
snv_chr3_tmp =  snv_chr3_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr3$next_snv = snv_chr3_tmp$start
snv_chr3$snv_distance = snv_chr3$next_snv - snv_chr3$start
snv_chr3 = snv_chr3 %>% filter(row_number() <= n()-1)

snv_chr4 = filter(snv_calls, chr =="chr4")
snv_chr4_tmp = snv_chr4[-1,]
snv_chr4_tmp =  snv_chr4_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr4$next_snv = snv_chr4_tmp$start
snv_chr4$snv_distance = snv_chr4$next_snv - snv_chr4$start
snv_chr4 = snv_chr4 %>% filter(row_number() <= n()-1)

snv_chr5 = filter(snv_calls, chr =="chr5")
snv_chr5_tmp = snv_chr5[-1,]
snv_chr5_tmp =  snv_chr5_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr5$next_snv = snv_chr5_tmp$start
snv_chr5$snv_distance = snv_chr5$next_snv - snv_chr5$start
snv_chr5 = snv_chr5 %>% filter(row_number() <= n()-1)

snv_chr6 = filter(snv_calls, chr =="chr6")
snv_chr6_tmp = snv_chr6[-1,]
snv_chr6_tmp =  snv_chr6_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr6$next_snv = snv_chr6_tmp$start
snv_chr6$snv_distance = snv_chr6$next_snv - snv_chr6$start
snv_chr6 = snv_chr6 %>% filter(row_number() <= n()-1)

snv_chr7 = filter(snv_calls, chr =="chr7")
snv_chr7_tmp = snv_chr7[-1,]
snv_chr7_tmp =  snv_chr7_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr7$next_snv = snv_chr7_tmp$start
snv_chr7$snv_distance = snv_chr7$next_snv - snv_chr7$start
snv_chr7 = snv_chr7 %>% filter(row_number() <= n()-1)

snv_chr8 = filter(snv_calls, chr =="chr8")
snv_chr8_tmp = snv_chr8[-1,]
snv_chr8_tmp =  snv_chr8_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr8$next_snv = snv_chr8_tmp$start
snv_chr8$snv_distance = snv_chr8$next_snv - snv_chr8$start
snv_chr8 = snv_chr8 %>% filter(row_number() <= n()-1)

snv_chr9 = filter(snv_calls, chr =="chr9")
snv_chr9_tmp = snv_chr9[-1,]
snv_chr9_tmp =  snv_chr9_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr9$next_snv = snv_chr9_tmp$start
snv_chr9$snv_distance = snv_chr9$next_snv - snv_chr9$start
snv_chr9 = snv_chr9 %>% filter(row_number() <= n()-1)

snv_chr10 = filter(snv_calls, chr =="chr10")
snv_chr10_tmp = snv_chr10[-1,]
snv_chr10_tmp =  snv_chr10_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr10$next_snv = snv_chr10_tmp$start
snv_chr10$snv_distance = snv_chr10$next_snv - snv_chr10$start
snv_chr10 = snv_chr10 %>% filter(row_number() <= n()-1)

snv_chr11 = filter(snv_calls, chr =="chr11")
snv_chr11_tmp = snv_chr11[-1,]
snv_chr11_tmp =  snv_chr11_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr11$next_snv = snv_chr11_tmp$start
snv_chr11$snv_distance = snv_chr11$next_snv - snv_chr11$start
snv_chr11 = snv_chr11 %>% filter(row_number() <= n()-1)

snv_chr12 = filter(snv_calls, chr =="chr12")
snv_chr12_tmp = snv_chr12[-1,]
snv_chr12_tmp =  snv_chr12_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr12$next_snv = snv_chr12_tmp$start
snv_chr12$snv_distance = snv_chr12$next_snv - snv_chr12$start
snv_chr12 = snv_chr12 %>% filter(row_number() <= n()-1)

snv_chr13 = filter(snv_calls, chr =="chr13")
snv_chr13_tmp = snv_chr13[-1,]
snv_chr13_tmp =  snv_chr13_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr13$next_snv = snv_chr13_tmp$start
snv_chr13$snv_distance = snv_chr13$next_snv - snv_chr13$start
snv_chr13 = snv_chr13 %>% filter(row_number() <= n()-1)

snv_chr14 = filter(snv_calls, chr =="chr14")
snv_chr14_tmp = snv_chr14[-1,]
snv_chr14_tmp =  snv_chr14_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr14$next_snv = snv_chr14_tmp$start
snv_chr14$snv_distance = snv_chr14$next_snv - snv_chr14$start
snv_chr14 = snv_chr14 %>% filter(row_number() <= n()-1)

snv_chr15 = filter(snv_calls, chr =="chr15")
snv_chr15_tmp = snv_chr15[-1,]
snv_chr15_tmp =  snv_chr15_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr15$next_snv = snv_chr15_tmp$start
snv_chr15$snv_distance = snv_chr15$next_snv - snv_chr15$start
snv_chr15 = snv_chr15 %>% filter(row_number() <= n()-1)

snv_chr16 = filter(snv_calls, chr =="chr16")
snv_chr16_tmp = snv_chr16[-1,]
snv_chr16_tmp =  snv_chr16_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr16$next_snv = snv_chr16_tmp$start
snv_chr16$snv_distance = snv_chr16$next_snv - snv_chr16$start
snv_chr16 = snv_chr16 %>% filter(row_number() <= n()-1)

snv_chr17 = filter(snv_calls, chr =="chr17")
snv_chr17_tmp = snv_chr17[-1,]
snv_chr17_tmp =  snv_chr17_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr17$next_snv = snv_chr17_tmp$start
snv_chr17$snv_distance = snv_chr17$next_snv - snv_chr17$start
snv_chr17 = snv_chr17 %>% filter(row_number() <= n()-1)

snv_chr18 = filter(snv_calls, chr =="chr18")
snv_chr18_tmp = snv_chr18[-1,]
snv_chr18_tmp =  snv_chr18_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr18$next_snv = snv_chr18_tmp$start
snv_chr18$snv_distance = snv_chr18$next_snv - snv_chr18$start
snv_chr18 = snv_chr18 %>% filter(row_number() <= n()-1)

snv_chr19 = filter(snv_calls, chr =="chr19")
snv_chr19_tmp = snv_chr19[-1,]
snv_chr19_tmp =  snv_chr19_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr19$next_snv = snv_chr19_tmp$start
snv_chr19$snv_distance = snv_chr19$next_snv - snv_chr19$start
snv_chr19 = snv_chr19 %>% filter(row_number() <= n()-1)

snv_chr20 = filter(snv_calls, chr =="chr20")
snv_chr20_tmp = snv_chr20[-1,]
snv_chr20_tmp =  snv_chr20_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr20$next_snv = snv_chr20_tmp$start
snv_chr20$snv_distance = snv_chr20$next_snv - snv_chr20$start
snv_chr20 = snv_chr20 %>% filter(row_number() <= n()-1)

snv_chr21 = filter(snv_calls, chr =="chr21")
snv_chr21_tmp = snv_chr21[-1,]
snv_chr21_tmp =  snv_chr21_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr21$next_snv = snv_chr21_tmp$start
snv_chr21$snv_distance = snv_chr21$next_snv - snv_chr21$start
snv_chr21 = snv_chr21 %>% filter(row_number() <= n()-1)

snv_chr22 = filter(snv_calls, chr =="chr22")
snv_chr22_tmp = snv_chr22[-1,]
snv_chr22_tmp =  snv_chr22_tmp %>% group_by(chr) %>% slice(c(1:n(),n()))
snv_chr22$next_snv = snv_chr22_tmp$start
snv_chr22$snv_distance = snv_chr22$next_snv - snv_chr22$start
snv_chr22 = snv_chr22 %>% filter(row_number() <= n()-1)

snv_dist_df = rbind(snv_chr1, snv_chr2, snv_chr3, snv_chr4, snv_chr5, snv_chr6, snv_chr7, snv_chr8, snv_chr9, snv_chr10, snv_chr11, snv_chr12, snv_chr13, snv_chr14, snv_chr15, snv_chr16, snv_chr17, snv_chr18, snv_chr19, snv_chr20, snv_chr21, snv_chr22)

#calculate Q3 (to omit outliers)
quant3 = quantile(snv_dist_df$snv_distance, probs = 0.75)
q3_out = as.integer(quant3)
snv_distance_long = snv_dist_df %>% 
  filter(snv_distance> q3_out)
n_filtered_snv = nrow(snv_distance_long)

cat ("\n        Cut-off value (outliers, >q3) for SNV distances:")
q3_out
cat ("\n        Number of SNVs removed (>q3):")
n_filtered_snv
cat ("\n")

snv_distance_long = snv_dist_df %>% 
  filter(snv_distance> q3_out)

#chrdist
#subset dataframe on sv subtype
snvs = filter(sv_calls, sv_type == "SNV")

#add chr count variable
snvs = snvs %>% add_count(chr)

#remove duplicate chr rows
snvs_count = distinct(snvs, chr, .keep_all = TRUE)

#set y max
ymax = max(snvs_count$n)

#small variant lenght distribution (del & dup)
#subset data to only include del and dup
sv_del = filter(sv_calls, sv_type == "del")
sv_dup = filter(sv_calls, sv_type == "dup")

#rbind two data frames
sv_deldup = rbind(sv_del, sv_dup)

#remove any potential variants larger than 50bp (=SV)
sv_deldup = sv_deldup[sv_deldup[,4]<50,]

#SNV ideogram
#subset all snv to new df
snv_calls = filter(sv_calls, sv_type == "SNV") %>%
  select(chr, start)

#set end coordinates (default to 30000 bp)
endcoordval = 30000

#duplicate snv df
snv_calls_ideo = snv_calls

snv_calls_ideo$end = snv_calls_ideo$start + endcoordval

#set chr as factor
snv_calls_ideo$chr = as.factor(snv_calls_ideo$chr)

#subset every chr, remove first line, duplicate last row, merge df, calculate snv distance (bp), remove lst row
snv_chr1_ideo = filter(snv_calls_ideo, chr =="chr1")
snv_chr2_ideo = filter(snv_calls_ideo, chr =="chr2")
snv_chr3_ideo = filter(snv_calls_ideo, chr =="chr3")
snv_chr4_ideo = filter(snv_calls_ideo, chr =="chr4")
snv_chr5_ideo = filter(snv_calls_ideo, chr =="chr5")
snv_chr6_ideo = filter(snv_calls_ideo, chr =="chr6")
snv_chr7_ideo = filter(snv_calls_ideo, chr =="chr7")
snv_chr8_ideo = filter(snv_calls_ideo, chr =="chr8")
snv_chr9_ideo = filter(snv_calls_ideo, chr =="chr9")
snv_chr10_ideo = filter(snv_calls_ideo, chr =="chr10")
snv_chr11_ideo = filter(snv_calls_ideo, chr =="chr11")
snv_chr12_ideo = filter(snv_calls_ideo, chr =="chr12")
snv_chr13_ideo = filter(snv_calls_ideo, chr =="chr13")
snv_chr14_ideo = filter(snv_calls_ideo, chr =="chr14")
snv_chr15_ideo = filter(snv_calls_ideo, chr =="chr15")
snv_chr16_ideo = filter(snv_calls_ideo, chr =="chr16")
snv_chr17_ideo = filter(snv_calls_ideo, chr =="chr17")
snv_chr18_ideo = filter(snv_calls_ideo, chr =="chr18")
snv_chr19_ideo = filter(snv_calls_ideo, chr =="chr19")
snv_chr20_ideo = filter(snv_calls_ideo, chr =="chr20")
snv_chr21_ideo = filter(snv_calls_ideo, chr =="chr21")
snv_chr22_ideo = filter(snv_calls_ideo, chr =="chr22")

#ideograms (SNVs)
#read dependencies for plotting
chr1.table = read.table("dep/ideograms/chrtable/GRCh38.chr1.table.txt", header = FALSE, sep = "\t")
chr2.table = read.table("dep/ideograms/chrtable/GRCh38.chr2.table.txt", header = FALSE, sep = "\t")
chr3.table = read.table("dep/ideograms/chrtable/GRCh38.chr3.table.txt", header = FALSE, sep = "\t")
chr4.table = read.table("dep/ideograms/chrtable/GRCh38.chr4.table.txt", header = FALSE, sep = "\t")
chr5.table = read.table("dep/ideograms/chrtable/GRCh38.chr5.table.txt", header = FALSE, sep = "\t")
chr6.table = read.table("dep/ideograms/chrtable/GRCh38.chr6.table.txt", header = FALSE, sep = "\t")
chr7.table = read.table("dep/ideograms/chrtable/GRCh38.chr7.table.txt", header = FALSE, sep = "\t")
chr8.table = read.table("dep/ideograms/chrtable/GRCh38.chr8.table.txt", header = FALSE, sep = "\t")
chr9.table = read.table("dep/ideograms/chrtable/GRCh38.chr9.table.txt", header = FALSE, sep = "\t")
chr10.table = read.table("dep/ideograms/chrtable/GRCh38.chr10.table.txt", header = FALSE, sep = "\t")
chr11.table = read.table("dep/ideograms/chrtable/GRCh38.chr11.table.txt", header = FALSE, sep = "\t")
chr12.table = read.table("dep/ideograms/chrtable/GRCh38.chr12.table.txt", header = FALSE, sep = "\t")
chr13.table = read.table("dep/ideograms/chrtable/GRCh38.chr13.table.txt", header = FALSE, sep = "\t")
chr14.table = read.table("dep/ideograms/chrtable/GRCh38.chr14.table.txt", header = FALSE, sep = "\t")
chr15.table = read.table("dep/ideograms/chrtable/GRCh38.chr15.table.txt", header = FALSE, sep = "\t")
chr16.table = read.table("dep/ideograms/chrtable/GRCh38.chr16.table.txt", header = FALSE, sep = "\t")
chr17.table = read.table("dep/ideograms/chrtable/GRCh38.chr17.table.txt", header = FALSE, sep = "\t")
chr18.table = read.table("dep/ideograms/chrtable/GRCh38.chr18.table.txt", header = FALSE, sep = "\t")
chr19.table = read.table("dep/ideograms/chrtable/GRCh38.chr19.table.txt", header = FALSE, sep = "\t")
chr20.table = read.table("dep/ideograms/chrtable/GRCh38.chr20.table.txt", header = FALSE, sep = "\t")
chr21.table = read.table("dep/ideograms/chrtable/GRCh38.chr21.table.txt", header = FALSE, sep = "\t")
chr22.table = read.table("dep/ideograms/chrtable/GRCh38.chr22.table.txt", header = FALSE, sep = "\t")
chr1.cent = read.table("dep/ideograms/centromeres/GRCh38.chr1.centromeres.txt", header = FALSE, sep = "\t")
chr2.cent = read.table("dep/ideograms/centromeres/GRCh38.chr2.centromeres.txt", header = FALSE, sep = "\t")
chr3.cent = read.table("dep/ideograms/centromeres/GRCh38.chr3.centromeres.txt", header = FALSE, sep = "\t")
chr4.cent = read.table("dep/ideograms/centromeres/GRCh38.chr4.centromeres.txt", header = FALSE, sep = "\t")
chr5.cent = read.table("dep/ideograms/centromeres/GRCh38.chr5.centromeres.txt", header = FALSE, sep = "\t")
chr6.cent = read.table("dep/ideograms/centromeres/GRCh38.chr6.centromeres.txt", header = FALSE, sep = "\t")
chr7.cent = read.table("dep/ideograms/centromeres/GRCh38.chr7.centromeres.txt", header = FALSE, sep = "\t")
chr8.cent = read.table("dep/ideograms/centromeres/GRCh38.chr8.centromeres.txt", header = FALSE, sep = "\t")
chr9.cent = read.table("dep/ideograms/centromeres/GRCh38.chr9.centromeres.txt", header = FALSE, sep = "\t")
chr10.cent = read.table("dep/ideograms/centromeres/GRCh38.chr10.centromeres.txt", header = FALSE, sep = "\t")
chr11.cent = read.table("dep/ideograms/centromeres/GRCh38.chr11.centromeres.txt", header = FALSE, sep = "\t")
chr12.cent = read.table("dep/ideograms/centromeres/GRCh38.chr12.centromeres.txt", header = FALSE, sep = "\t")
chr13.cent = read.table("dep/ideograms/centromeres/GRCh38.chr13.centromeres.txt", header = FALSE, sep = "\t")
chr14.cent = read.table("dep/ideograms/centromeres/GRCh38.chr14.centromeres.txt", header = FALSE, sep = "\t")
chr15.cent = read.table("dep/ideograms/centromeres/GRCh38.chr15.centromeres.txt", header = FALSE, sep = "\t")
chr16.cent = read.table("dep/ideograms/centromeres/GRCh38.chr16.centromeres.txt", header = FALSE, sep = "\t")
chr17.cent = read.table("dep/ideograms/centromeres/GRCh38.chr17.centromeres.txt", header = FALSE, sep = "\t")
chr18.cent = read.table("dep/ideograms/centromeres/GRCh38.chr18.centromeres.txt", header = FALSE, sep = "\t")
chr19.cent = read.table("dep/ideograms/centromeres/GRCh38.chr19.centromeres.txt", header = FALSE, sep = "\t")
chr20.cent = read.table("dep/ideograms/centromeres/GRCh38.chr20.centromeres.txt", header = FALSE, sep = "\t")
chr21.cent = read.table("dep/ideograms/centromeres/GRCh38.chr21.centromeres.txt", header = FALSE, sep = "\t")
chr22.cent = read.table("dep/ideograms/centromeres/GRCh38.chr22.centromeres.txt", header = FALSE, sep = "\t")

#plot
sv_size_violine = ggplot(sv_deldup, aes(x = sv_type, y = sv_length, fill = sv_type)) + 
  labs(title = "Small Variants Size Distribution", subtitle = "Violin plot visualizing size-distributions for small variants (i.g deletions & duplications\n≤ 50bp) with SNVs excluded. Variant size (1-50bp) are arranged on the y-axis and\nvariant sub-type on the x-axis. Black dot annotates mean variant-size.", x = "", y = "Size (bp)") +
  geom_violin(trim = FALSE, scale = "width", color = NA) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +
  scale_y_continuous(breaks = seq(0, 50, by = 10)) +
  theme(legend.position = "none", plot.margin=unit(c(0,1.3,0,1.3),"cm"))

#chromosome distribution box plot
snv_chrdist_box = ggplot(snvs_count, aes(x = chr, y = n)) +
  labs(title = "SNV per Chromosome", subtitle = "Flipped histogram, showing the number of SNVs per chromosome.\nChromosomes are plotted on the y-axis and number of SNVs on the x-axis.\nExpected SNV frequency is on average 1 variant for every 1200 bases.\nExpected variant frequency for each chromosome are highlighted with grey dots.\n(Chaisson, Mark J P et al. “Multi-platform discovery of haplotype-resolved\nstructural variation in human genomes.”).", x = "", y = "Count (n)", fill = "") +
  scale_x_discrete(limits=c("chr22", "chr21", "chr20", "chr19", "chr18", "chr17", "chr16", "chr15", "chr14", "chr13", "chr12", "chr11", "chr10", "chr9", "chr8", "chr7", "chr6", "chr5", "chr4", "chr3", "chr2", "chr1")) +
  geom_bar(position = "stack", stat = "identity", fill = "#6FB392", width = 0.5) +
  theme(plot.margin=unit(c(0,1.3,0,1.3),"cm")) +
  coord_flip()

#density plot of SNV distances
snv_distance_plot = ggplot(snv_dist_df, aes(x = chr, y = snv_distance)) +
  geom_boxplot(outlier.shape = NA, fill = "#6FB392") + 
  labs(title = "SNV Distance", subtitle = "Histogram showing the distribution of distances between neighboring SNVs. SNV-distances above the 3rd quantile are excluded, to compensate for neighbouring variants with exceptionally\nlong distance between each other (e.g variants on opposite sides of the same centromere). Mean SNV distance for each chromosome is shown with black line inside each box.\nMetric can be used to interrogate the breadth of called SNVs compared to an expected frequency (SNV occurs on average every 1200 nucleotide). Chromosomes with a strong deviation\nfrom the sample-specific average SNV-distance will be highlited, revealing if SNV-calling is skewed towards certain chromosomes.", x = "", y = "Distance (bp)", fill = "") +
  scale_x_discrete(limits=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.margin=unit(c(0,1.3,0,1.3),"cm")) +
  ylim(0, quant3)

#plot header for report
text = paste0(sample_name, " | GRCh38 | Small Variants | ", now)
plot.title = ggplot() + 
  annotate("text", x = 7, y = 3, size=8, label = text) + 
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = "white"))

#plot ideograms
chr1.ideo = ggplot() +
  #chr table
  geom_segment(data = chr1.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr1.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr1_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr1 ", size = 8)

chr2.ideo = ggplot() +
  #chr table
  geom_segment(data = chr2.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr2.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr2_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr2 ", size = 8)

chr3.ideo = ggplot() +
  #chr table
  geom_segment(data = chr3.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr3.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr3_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr3 ", size = 8)

chr4.ideo = ggplot() +
  #chr table
  geom_segment(data = chr4.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr4.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr4_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr4 ", size = 8)

chr5.ideo = ggplot() +
  #chr table
  geom_segment(data = chr5.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr5.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr5_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr5 ", size = 8)

chr6.ideo = ggplot() +
  #chr table
  geom_segment(data = chr6.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr6.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr6_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr6 ", size = 8)

chr7.ideo = ggplot() +
  #chr table
  geom_segment(data = chr7.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr7.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr7_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr7 ", size = 8)

chr8.ideo = ggplot() +
  #chr table
  geom_segment(data = chr8.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr8.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr8_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr8 ", size = 8)

chr9.ideo = ggplot() +
  #chr table
  geom_segment(data = chr9.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr9.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr9_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr9 ", size = 8)

chr10.ideo = ggplot() +
  #chr table
  geom_segment(data = chr10.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr10.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr10_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr10 ", size = 8)

chr11.ideo = ggplot() +
  #chr table
  geom_segment(data = chr11.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr11.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr11_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr11 ", size = 8)

chr12.ideo = ggplot() +
  #chr table
  geom_segment(data = chr12.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr12.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr12_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr12 ", size = 8)

chr13.ideo = ggplot() +
  #chr table
  geom_segment(data = chr13.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr13.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr13_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr13 ", size = 8)

chr14.ideo = ggplot() +
  #chr table
  geom_segment(data = chr14.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr14.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr14_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr14 ", size = 8)

chr15.ideo = ggplot() +
  #chr table
  geom_segment(data = chr15.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr15.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr15_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr15 ", size = 8)

chr16.ideo = ggplot() +
  #chr table
  geom_segment(data = chr16.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr16.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr16_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr16 ", size = 8)

chr17.ideo = ggplot() +
  #chr table
  geom_segment(data = chr17.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr17.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr17_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr17 ", size = 8)

chr18.ideo = ggplot() +
  #chr table
  geom_segment(data = chr18.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr18.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr18_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr18 ", size = 8)

chr19.ideo = ggplot() +
  #chr table
  geom_segment(data = chr19.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr19.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr19_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr19 ", size = 8)

chr20.ideo = ggplot() +
  #chr table
  geom_segment(data = chr20.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr20.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr20_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr20 ", size = 8)

chr21.ideo = ggplot() +
  #chr table
  geom_segment(data = chr21.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr21.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr21_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr21 ", size = 8)

chr22.ideo = ggplot() +
  #chr table
  geom_segment(data = chr22.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey29", lineend = "butt", size = 10, stat = "identity", position = position_dodge()) +
  #centromere
  geom_segment(data = chr22.cent, aes(x = V2, xend = V3, y = 0, yend = 0), color = "white", size = 11, stat = "identity", position = position_dodge()) +
  #snvs
  geom_segment(data = snv_chr22_ideo, aes(x = start, xend = end, y = 0, yend = 0), color = "#6FB392", size = 10, stat = "identity", position = position_dodge()) + 
  #theme
  theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
  ylab("") + 
  xlab("") +
  ylim(-70, 70) +
  xlim(-15000000, 248956422) +
  annotate(geom = "text", x = -15000000, y = 0, label="chr22 ", size = 8)

#plot title
text.ideo = "SNV ideogram"

#calculate snv coverage (subtitle 1)
snv_cov = (sum(snvs_count$n) / 3088269832) * 100
snv_cov = format(round(snv_cov, 5), nsmall = 5)

#caculate refN (subtitle 2)
#read vcf into r
vcf.list = list.files(path = "in/small_variants/", recursive = TRUE, pattern = "\\.vcf$", full.names = TRUE)
sample_name = gsub(".{4}$", '', vcf.list)
sample_name = substring(sample_name, 20)
vcf_refcount = read.table(file = paste0("in/small_variants/", sample_name, ".vcf"), sep = "\t", header = F, comment.char="#")

#count nucleotides in ref 
vcf_refcount$refcount = nchar(vcf_refcount$V4)
refN = sum(vcf_refcount$refcount)

#plot title with subtitles
plot.title.ideo = ggplot() + 
  annotate("text", x = 1, y = 10, size = 10, label = text.ideo, fontface = "bold", hjust = 0) +
  annotate("text", x = 1, y = 6, size = 8, label = paste0("Number of reference (REF) nucleotides in VCF: ", refN), fontface = "bold", hjust = 0) +
  annotate("text", x = 1, y = 3, size = 8, label = "Chromosome dependant ideogram, SNVs are plotted (superimposed) as green vertical lines along each autosome (grey), revealing chromosomal regions with decreased SNV density.", hjust = 0) +
  annotate("text", x = 1, y = 0, size = 8, label = "regions with decreased SNV density.", hjust = 0) +
  xlim(0, 100) +
  ylim(0,20) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank())

#plot empty box for spacing
box = ggplot() +
  geom_blank(mapping = NULL, data = NULL, stat = "identity", position = "identity", show.legend = NA, inherit.aes = TRUE) + 
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(), panel.background = element_rect(fill = "white", colour = "white"))

#export plots
ggsave(sv_size_violine, file = paste0("out/small_variants/figs/", txtFileName, "_01_smallvariants_size_violin.png"), limitsize = FALSE, width = 7, height = 7, units = c("in"), dpi = 300)
ggsave(snv_distance_plot, file = paste0("out/small_variants/figs/", txtFileName, "_02_snv_distance.png"), limitsize = FALSE, width = 14, height = 7, units = c("in"), dpi = 300)
ggsave(snv_chrdist_box, file = paste0("out/small_variants/figs/", txtFileName, "_03_snv_chr_dist.png"), limitsize = FALSE, width = 7, height = 7, units = c("in"), dpi = 300)
ggsave(plot.title, file = paste0("out/small_variants/figs/", txtFileName, "_header.png"), limitsize = FALSE, width = 14, height = 1, units = c("in"), dpi = 300)
ggsave("plot.title.png", plot.title.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 80, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr1.png", chr1.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr2.png", chr2.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr3.png", chr3.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr4.png", chr4.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr5.png", chr5.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr6.png", chr6.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr7.png", chr7.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr8.png", chr8.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr9.png", chr9.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr10.png", chr10.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr11.png", chr11.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr12.png", chr12.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr13.png", chr13.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr14.png", chr14.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr15.png", chr15.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr16.png", chr16.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr17.png", chr17.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr18.png", chr18.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr19.png", chr19.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr20.png", chr20.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr21.png", chr21.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave("chr22.png", chr22.ideo, path = "out/small_variants/figs/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
ggsave(box, file = paste0("out/small_variants/figs/", txtFileName, "_box2.png"), limitsize = FALSE, width = 14, height = 0.3, units = c("in"), dpi = 300)

#prompt message to terminal
cat ("Figures and summaries (small variants) exported...\n")
