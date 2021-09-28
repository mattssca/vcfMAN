#load packages
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggthemr))
suppressMessages(library(gridExtra))
suppressMessages(library(RCircos))
suppressMessages(library(ggthemr))

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

#plot
sv_size_violine = ggplot(sv_deldup, aes(x = sv_type, y = sv_length, fill = sv_type)) + 
  labs(title = "Small Variants Size Distribution", x = "", y = "Size (bp)") +
  geom_violin(trim = FALSE, scale = "width", color = NA) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +
  theme(legend.position = "none") +
  scale_y_log10()

#chromosome distribution box plot
snv_chrdist_box = ggplot(snvs_count, aes(x = chr, y = n)) +
  labs(title = "SNV per Chromosome", subtitle = "Ref: 1 SNV Every 1000bp",x = "", y = "Count (n)", fill = "") +
  scale_x_discrete(limits=c("chr22", "chr21", "chr20", "chr19", "chr18", "chr17", "chr16", "chr15", "chr14", "chr13", "chr12", "chr11", "chr10", "chr9", "chr8", "chr7", "chr6", "chr5", "chr4", "chr3", "chr2", "chr1")) +
  geom_bar(position = "stack", stat = "identity", fill = "#3BA87B", width = 0.5) +
  coord_flip()

#density plot of SNV distances
snv_distance_plot = ggplot(snv_dist_df, aes(x = chr, y = snv_distance)) +
  geom_boxplot(outlier.shape = NA, fill = "#6E9BF5") + 
  labs(title = "SNV Distance", x = "", y = "Distance (bp)", fill = "") +
  scale_x_discrete(limits=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  ylim(0, quant3)

#plot header for report
text = paste0(sample_name, " | GRCh38 | Small Variants | ", now)
plot.title = ggplot() + 
  annotate("text", x = 7, y = 3, size=8, label = text) + 
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = "white"))

#export plots
ggsave(sv_size_violine, file = paste0("out/small_variants/figs/", txtFileName, "_01_smallvariants_size_violin.png"), limitsize = FALSE, width = 7, height = 7, units = c("in"), dpi = 300)
ggsave(snv_distance_plot, file = paste0("out/small_variants/figs/", txtFileName, "_02_snv_distance.png"), limitsize = FALSE, width = 14, height = 7, units = c("in"), dpi = 300)
ggsave(snv_chrdist_box, file = paste0("out/small_variants/figs/", txtFileName, "_03_snv_chr_dist.png"), limitsize = FALSE, width = 7, height = 7, units = c("in"), dpi = 300)
ggsave(plot.title, file = paste0("out/small_variants/figs/", txtFileName, "_header.png"), limitsize = FALSE, width = 14, height = 1, units = c("in"), dpi = 300)

