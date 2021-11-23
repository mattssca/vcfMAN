#R script for processing BED files (dip variants).
#author, Carl-Adam Mattsson. omYXX informatics.

#remove warning messages
options(warn=-1)

#load packages
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggthemr))

#get sample name
bed.list = list.files(path = "in/BED/", recursive = TRUE, pattern = "\\.bed$", full.names = TRUE)
sample_name = gsub(".{4}$", '', bed.list)
sample_name = substring(sample_name, 9)

#set parameters
now = format(Sys.time(), "%d_%m_%Y")
txtFileName = paste0(sample_name, "_bedspread_", now)

#read vcf into R (skip header)
bed_dip = read.table(file = paste0("in/BED/", sample_name, ".bed"), sep = "\t", header = F)

#subset on autosomes
selected = paste0("chr", c(1:22))
#selected = paste0("chr", c(1:22, "x,y"))
bed_dip = bed_dip[bed_dip$V1 %in% selected,]

#calculate variant-size
bed_dip$V4 = bed_dip$V3 - bed_dip$V2

#create new variable for sv-type (duplicate size variable, for now)
bed_dip$V5 = bed_dip$V4

#define variable names
bed_dip = setNames(bed_dip, c("chr", "start", "end", "size", "type"))

#compute variant type (deletion/insertion)
bed_dip$type[which(bed_dip$size > 50)] = "structural_variant"
bed_dip$type[which(bed_dip$size < 50)] = "dip_variant"

#transform chr and variant-type to factor
bed_dip$chr = as.factor(bed_dip$chr)
bed_dip$type = as.factor(bed_dip$type)

#sort bed file
bed_dip = bed_dip[order(bed_dip$chr, bed_dip$start),]

#compute number of bases included in bed file
n_bases_dip = sum(bed_dip$size)

#chr tables
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

#read blacklsited regions
chr1.black = read.table("dep/ideograms/black/GRCh38.chr1.black.txt", header = FALSE, sep = "\t")
chr2.black = read.table("dep/ideograms/black/GRCh38.chr2.black.txt", header = FALSE, sep = "\t")
chr3.black = read.table("dep/ideograms/black/GRCh38.chr3.black.txt", header = FALSE, sep = "\t")
chr4.black = read.table("dep/ideograms/black/GRCh38.chr4.black.txt", header = FALSE, sep = "\t")
chr5.black = read.table("dep/ideograms/black/GRCh38.chr5.black.txt", header = FALSE, sep = "\t")
chr6.black = read.table("dep/ideograms/black/GRCh38.chr6.black.txt", header = FALSE, sep = "\t")
chr7.black = read.table("dep/ideograms/black/GRCh38.chr7.black.txt", header = FALSE, sep = "\t")
chr8.black = read.table("dep/ideograms/black/GRCh38.chr8.black.txt", header = FALSE, sep = "\t")
chr9.black = read.table("dep/ideograms/black/GRCh38.chr9.black.txt", header = FALSE, sep = "\t")
chr10.black = read.table("dep/ideograms/black/GRCh38.chr10.black.txt", header = FALSE, sep = "\t")
chr11.black = read.table("dep/ideograms/black/GRCh38.chr11.black.txt", header = FALSE, sep = "\t")
chr12.black = read.table("dep/ideograms/black/GRCh38.chr12.black.txt", header = FALSE, sep = "\t")
chr13.black = read.table("dep/ideograms/black/GRCh38.chr13.black.txt", header = FALSE, sep = "\t")
chr14.black = read.table("dep/ideograms/black/GRCh38.chr14.black.txt", header = FALSE, sep = "\t")
chr15.black = read.table("dep/ideograms/black/GRCh38.chr15.black.txt", header = FALSE, sep = "\t")
chr16.black = read.table("dep/ideograms/black/GRCh38.chr16.black.txt", header = FALSE, sep = "\t")
chr17.black = read.table("dep/ideograms/black/GRCh38.chr17.black.txt", header = FALSE, sep = "\t")
chr18.black = read.table("dep/ideograms/black/GRCh38.chr18.black.txt", header = FALSE, sep = "\t")
chr19.black = read.table("dep/ideograms/black/GRCh38.chr19.black.txt", header = FALSE, sep = "\t")
chr20.black = read.table("dep/ideograms/black/GRCh38.chr20.black.txt", header = FALSE, sep = "\t")
chr21.black = read.table("dep/ideograms/black/GRCh38.chr21.black.txt", header = FALSE, sep = "\t")
chr22.black = read.table("dep/ideograms/black/GRCh38.chr22.black.txt", header = FALSE, sep = "\t")

#add all blacklisted regions to single df
blacklist = rbind(chr1.black, chr2.black, chr3.black, chr4.black, chr5.black, chr6.black, chr7.black, chr8.black, chr9.black, chr10.black, chr11.black, chr12.black, chr13.black, chr14.black, chr15.black, chr16.black, chr17.black, chr18.black, chr19.black, chr20.black, chr21.black, chr22.black)

#calculate size of blacklsited region
blacklist$V4 = blacklist$V3 - blacklist$V2

#add all size values, to get tot n bases blacklisted
black_sum = sum(blacklist$V4)

#calculate size of blacklisted regions for individual chromosomes
chr1.black$V4 = chr1.black$V3 - chr1.black$V2
chr2.black$V4 = chr2.black$V3 - chr2.black$V2
chr3.black$V4 = chr3.black$V3 - chr3.black$V2
chr4.black$V4 = chr4.black$V3 - chr4.black$V2
chr5.black$V4 = chr5.black$V3 - chr5.black$V2
chr6.black$V4 = chr6.black$V3 - chr6.black$V2
chr7.black$V4 = chr7.black$V3 - chr7.black$V2
chr8.black$V4 = chr8.black$V3 - chr8.black$V2
chr9.black$V4 = chr9.black$V3 - chr9.black$V2
chr10.black$V4 = chr10.black$V3 - chr10.black$V2
chr11.black$V4 = chr11.black$V3 - chr11.black$V2
chr12.black$V4 = chr12.black$V3 - chr12.black$V2
chr13.black$V4 = chr13.black$V3 - chr13.black$V2
chr14.black$V4 = chr14.black$V3 - chr14.black$V2
chr15.black$V4 = chr15.black$V3 - chr15.black$V2
chr16.black$V4 = chr16.black$V3 - chr16.black$V2
chr17.black$V4 = chr17.black$V3 - chr17.black$V2
chr18.black$V4 = chr18.black$V3 - chr18.black$V2
chr19.black$V4 = chr19.black$V3 - chr19.black$V2
chr20.black$V4 = chr20.black$V3 - chr20.black$V2
chr21.black$V4 = chr21.black$V3 - chr21.black$V2
chr22.black$V4 = chr22.black$V3 - chr22.black$V2

#summarize all regions for each chromosome
chr1.black.sum = sum(chr1.black$V4)
chr2.black.sum = sum(chr2.black$V4)
chr3.black.sum = sum(chr3.black$V4)
chr4.black.sum = sum(chr4.black$V4)
chr5.black.sum = sum(chr5.black$V4)
chr6.black.sum = sum(chr6.black$V4)
chr7.black.sum = sum(chr7.black$V4)
chr8.black.sum = sum(chr8.black$V4)
chr9.black.sum = sum(chr9.black$V4)
chr10.black.sum = sum(chr10.black$V4)
chr11.black.sum = sum(chr11.black$V4)
chr12.black.sum = sum(chr12.black$V4)
chr13.black.sum = sum(chr13.black$V4)
chr14.black.sum = sum(chr14.black$V4)
chr15.black.sum = sum(chr15.black$V4)
chr16.black.sum = sum(chr16.black$V4)
chr17.black.sum = sum(chr17.black$V4)
chr18.black.sum = sum(chr18.black$V4)
chr19.black.sum = sum(chr19.black$V4)
chr20.black.sum = sum(chr20.black$V4)
chr21.black.sum = sum(chr21.black$V4)
chr22.black.sum = sum(chr22.black$V4)

#total number of bases in grch38
n_bases_grch38 = 3058848640 - black_sum

#calculate fraction of bases include in input BED file
dip_var_frac = n_bases_dip / n_bases_grch38 * 100
dip_var_frac = round(dip_var_frac, digits = 0)

#breakdown by chromosome
chr1_dip_bed = filter(bed_dip, chr =="chr1")
chr1_frac_dip = sum(chr1_dip_bed$size) / (chr1.table$V3 - chr1.black.sum) * 100
chr1_frac_dip = round(chr1_frac_dip, digits = 0)

chr2_dip_bed = filter(bed_dip, chr =="chr2")
chr2_frac_dip = sum(chr2_dip_bed$size) / (chr2.table$V3 - chr2.black.sum) * 100
chr2_frac_dip = round(chr2_frac_dip, digits = 0)

chr3_dip_bed = filter(bed_dip, chr =="chr3")
chr3_frac_dip = sum(chr3_dip_bed$size) / (chr3.table$V3 - chr3.black.sum) * 100
chr3_frac_dip = round(chr3_frac_dip, digits = 0)

chr4_dip_bed = filter(bed_dip, chr =="chr4")
chr4_frac_dip = sum(chr4_dip_bed$size) / (chr4.table$V3 - chr4.black.sum) * 100
chr4_frac_dip = round(chr4_frac_dip, digits = 0)

chr5_dip_bed = filter(bed_dip, chr =="chr5")
chr5_frac_dip = sum(chr5_dip_bed$size) / (chr5.table$V3 - chr5.black.sum) * 100
chr5_frac_dip = round(chr5_frac_dip, digits = 0)

chr6_dip_bed = filter(bed_dip, chr =="chr6")
chr6_frac_dip = sum(chr6_dip_bed$size) / (chr6.table$V3 - chr6.black.sum) * 100
chr6_frac_dip = round(chr6_frac_dip, digits = 0)

chr7_dip_bed = filter(bed_dip, chr =="chr7")
chr7_frac_dip = sum(chr7_dip_bed$size) / (chr7.table$V3 - chr7.black.sum) * 100
chr7_frac_dip = round(chr7_frac_dip, digits = 0)

chr8_dip_bed = filter(bed_dip, chr =="chr8")
chr8_frac_dip = sum(chr8_dip_bed$size) / (chr8.table$V3 - chr8.black.sum) * 100
chr8_frac_dip = round(chr8_frac_dip, digits = 0)

chr9_dip_bed = filter(bed_dip, chr =="chr9")
chr9_frac_dip = sum(chr9_dip_bed$size) / (chr9.table$V3 - chr9.black.sum) * 100
chr9_frac_dip = round(chr9_frac_dip, digits = 0)

chr10_dip_bed = filter(bed_dip, chr =="chr10")
chr10_frac_dip = sum(chr10_dip_bed$size) / (chr10.table$V3 - chr10.black.sum) * 100
chr10_frac_dip = round(chr10_frac_dip, digits = 0)

chr11_dip_bed = filter(bed_dip, chr =="chr11")
chr11_frac_dip = sum(chr11_dip_bed$size) / (chr11.table$V3 - chr11.black.sum) * 100
chr11_frac_dip = round(chr11_frac_dip, digits = 0)

chr12_dip_bed = filter(bed_dip, chr =="chr12")
chr12_frac_dip = sum(chr12_dip_bed$size) / (chr12.table$V3 - chr12.black.sum) * 100
chr12_frac_dip = round(chr12_frac_dip, digits = 0)

chr13_dip_bed = filter(bed_dip, chr =="chr13")
chr13_frac_dip = sum(chr13_dip_bed$size) / (chr13.table$V3 - chr13.black.sum) * 100
chr13_frac_dip = round(chr13_frac_dip, digits = 0)

chr14_dip_bed = filter(bed_dip, chr =="chr14")
chr14_frac_dip = sum(chr14_dip_bed$size) / (chr14.table$V3 - chr14.black.sum) * 100
chr14_frac_dip = round(chr14_frac_dip, digits = 0)

chr15_dip_bed = filter(bed_dip, chr =="chr15")
chr15_frac_dip = sum(chr15_dip_bed$size) / (chr15.table$V3 - chr15.black.sum) * 100
chr15_frac_dip = round(chr15_frac_dip, digits = 0)

chr16_dip_bed = filter(bed_dip, chr =="chr16")
chr16_frac_dip = sum(chr16_dip_bed$size) / (chr16.table$V3 - chr16.black.sum) * 100
chr16_frac_dip = round(chr16_frac_dip, digits = 0)

chr17_dip_bed = filter(bed_dip, chr =="chr17")
chr17_frac_dip = sum(chr17_dip_bed$size) / (chr17.table$V3 - chr17.black.sum) * 100
chr17_frac_dip = round(chr17_frac_dip, digits = 0)

chr18_dip_bed = filter(bed_dip, chr =="chr18")
chr18_frac_dip = sum(chr18_dip_bed$size) / (chr18.table$V3 - chr18.black.sum) * 100
chr18_frac_dip = round(chr18_frac_dip, digits = 0)

chr19_dip_bed = filter(bed_dip, chr =="chr19")
chr19_frac_dip = sum(chr19_dip_bed$size) / (chr19.table$V3 - chr19.black.sum) * 100
chr19_frac_dip = round(chr19_frac_dip, digits = 0)

chr20_dip_bed = filter(bed_dip, chr =="chr20")
chr20_frac_dip = sum(chr20_dip_bed$size) / (chr20.table$V3 - chr20.black.sum) * 100
chr20_frac_dip = round(chr20_frac_dip, digits = 0)

chr21_dip_bed = filter(bed_dip, chr =="chr21")
chr21_frac_dip = sum(chr21_dip_bed$size) / (chr21.table$V3 - chr21.black.sum) * 100
chr21_frac_dip = round(chr21_frac_dip, digits = 0)

chr22_dip_bed = filter(bed_dip, chr =="chr22")
chr22_frac_dip = sum(chr22_dip_bed$size) / (chr22.table$V3 - chr22.black.sum) * 100
chr22_frac_dip = round(chr22_frac_dip, digits = 0)

#construct data frame with chromosome coverage
chr_fractions = rbind(chr1_frac_dip, chr2_frac_dip, chr3_frac_dip, chr4_frac_dip, chr5_frac_dip, chr6_frac_dip, chr7_frac_dip, chr8_frac_dip, chr9_frac_dip, chr10_frac_dip, chr11_frac_dip, chr12_frac_dip, chr13_frac_dip, chr14_frac_dip, chr15_frac_dip, chr16_frac_dip, chr17_frac_dip, chr18_frac_dip, chr19_frac_dip, chr20_frac_dip, chr21_frac_dip, chr22_frac_dip)
chrtable = rbind("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")

#cbind both tables
chr_fractions = cbind(chrtable, chr_fractions)

#remove row names
rownames(chr_fractions) = NULL

#convert to data frame
chr_fractions = as.data.frame(chr_fractions)

#convert variables to correct type
chr_fractions$V1 = as.factor(chr_fractions$V1)
chr_fractions$V2 = as.numeric(chr_fractions$V2)

#set plot subtitle
text1 = paste0("For each chromosome, number of bases included in BED file, divided by total number of bases for that chromosome, expressed in percentages.\nThe total n bases in BED file as a fraction of the complete genome (grch38) is: ", dip_var_frac, "% and total number of bases in input BED file is: ", n_bases_dip, ".")

#create a new row with bases in BED as fraction of complete genome
chr_fractions = chr_fractions %>% add_row(V1 = "All", V2 = dip_var_frac)

#define plot theme
ggthemr("dust")

#plot chr BED coverage as fractions
bed_chrdist_frac = ggplot(chr_fractions, aes(x = V1, y = V2)) +
  labs(title = "BED (dip variants) coverage per Chromosome", subtitle = text1, x = "", y = "Fraction", fill = "") +
  scale_x_discrete(limits=c("All", "chr22", "chr21", "chr20", "chr19", "chr18", "chr17", "chr16", "chr15", "chr14", "chr13", "chr12", "chr11", "chr10", "chr9", "chr8", "chr7", "chr6", "chr5", "chr4", "chr3", "chr2", "chr1")) +
  geom_bar(position = "stack", stat = "identity", fill = "#714E91", width = 0.5) +
  theme(plot.margin=unit(c(1.3,1.3,1.3,1.3),"cm")) +
  geom_text(aes(label = V2), vjust = 0.5, hjust = 1.2, size = 3, colour = "white") + 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  coord_flip()

#export plot
ggsave(bed_chrdist_frac, file = "out/BED/bed_chrdist_frac.png", limitsize = FALSE, width = 14, height = 7, units = c("in"), dpi = 300)
