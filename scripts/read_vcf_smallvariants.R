#R script for processing VCF for small variants.
#author, Carl-Adam Mattsson. omYXX informatics.

#remove warning messages
options(warn=-1)

#load packages
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(knitr))
suppressMessages(library(rvest))
suppressMessages(library(table1))
suppressMessages(library(openxlsx))
suppressMessages(library(gridExtra))
suppressMessages(library(ggplot2))

#get sample name
vcf.list = list.files(path = "in/small_variants/", recursive = TRUE, pattern = "\\.vcf$", full.names = TRUE)
sample_name = gsub(".{4}$", '', vcf.list)
sample_name = substring(sample_name, 20)

#set parameters
now = format(Sys.time(), "%d_%m_%Y")
txtFileName <- paste0(sample_name, "_vcf_man_", now)

#read vcf into R (skip header)
vcf = read.table(file = paste0("in/small_variants/", sample_name, ".vcf"), sep = "\t", header = F, comment.char="#")

#drop chrX and Y
#subset on autosomes
selected = paste0("chr", c(1:22))
#selected = paste0("chr", c(1:22, "x,y"))
vcf = vcf[vcf$V1 %in% selected,]

#extract vcf head
vcfHead = readLines(paste0("in/small_variants/", sample_name, ".vcf"))
vcfHead = vcfHead[grep('^##.*', vcfHead)]
vcfHead = as.data.frame(vcfHead)

#subset vcf file
vcf_sum = vcf %>% select(V1, V2, V4, V5, V6, V10)

#format genotype field
vcf_sum$V10 = sub(":.+", "", vcf_sum$V10)

#define variable names
vcf_sum = setNames(vcf_sum, c("chr", "start", "ref", "alt", "qual", "genotype"))

#count nucleotides in ref and alt sequence (to get variant length)
vcf_sum$ref_count = nchar(vcf_sum$ref)
vcf_sum$alt_count = nchar(vcf_sum$alt)
vcf_sum$sv_length = vcf_sum$alt_count - vcf_sum$ref_count

#calculate variant end coordinates
vcf_sum$end = vcf_sum$start + vcf_sum$sv_length

#compute variant type (deletion/insertion)
vcf_sum$sv_type = vcf_sum$sv_length
vcf_sum$sv_type[which(vcf_sum$sv_length == "0")] = "SNV"
vcf_sum$sv_type[which(vcf_sum$sv_length < 0)] = "del"
vcf_sum$sv_type[which(vcf_sum$sv_length > 0)] = "ins"

#create new data frame with selected variables (prepare for plotting and exporting BED)
vcf_sub = vcf_sum %>% select(chr, start, end, sv_length, sv_type, genotype)

#transform genotype to factor with set levels
vcf_sub$genotype = as.factor(vcf_sub$genotype)
vcf_sub$sv_type = as.factor(vcf_sub$sv_type)
vcf_sub$chr = as.factor(vcf_sub$chr)

#subset on autosomes
vcf_sub = vcf_sub[vcf_sub$chr %in% selected,]

#duplicate filtered vcf
vcf_up = vcf_sub

#transform haplotype information to het/hom
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="1|1"] = "hom"
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="1|0"] = "het"
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="0|1"] = "het"

#subset and remove non-hardcoded genotypes from main df
non_hardc_gen.1 = filter(vcf_sub, genotype != "het")
non_hardc_gen = filter(non_hardc_gen.1, genotype != "hom")

#create variable for non-hardcoded genotypes and print n variants
n_non_hardc_gen = nrow(non_hardc_gen)
cat ("\n        Number of variants with non hardcoded genotype")
n_non_hardc_gen

#select harcoded genotypes
vcf_sub = filter(vcf_sub, genotype == c("het", "hom"))
vcf_sub = droplevels(vcf_sub)

#replace end-coordinates for deletions with NA
vcf_sub$end[which(vcf_sub$sv_type == "del")] = "NA"
vcf_sub$end[which(vcf_sub$sv_type == "SNV")] = "NA"

#transform into absolute values
vcf_sub$sv_length = abs(vcf_sub$sv_length)

#transform end coordinates for deletions and SNVs (+1)
sv_tmp_del = filter(vcf_sub, sv_type == "del")
sv_tmp_del$end = sv_tmp_del$start + 1

sv_tmp_snv = filter(vcf_sub, sv_type == "SNV")
sv_tmp_snv$end = sv_tmp_snv$start + 1
sv_tmp_snv$sv_length = 1

#remove deletions from main data frame and add in transformed end coordinates for deletions
vcf_sub = droplevels(vcf_sub[!vcf_sub$sv_type == "del",])
vcf_sub = droplevels(vcf_sub[!vcf_sub$sv_type == "SNV",])
vcf_sub = rbind(vcf_sub, sv_tmp_del, sv_tmp_snv)
vcf_sub = vcf_sub[order(vcf_sub$chr, vcf_sub$start),]

#subset variants larger than 50 bp
vcf_sub_large = subset(vcf_sub, sv_length > 50)

#create variable for number of variants smaller or equal to 50 bp
n_largevariants = nrow(vcf_sub_large)
cat ("\n        Number of variants > 50bp:")
n_largevariants

#filter out variants large than 50 bp
vcf_sub = subset(vcf_sub, sv_length <= 50)

#revert vcf to standard formatted genotype
vcf_up$end[which(vcf_up$sv_type == "del")] = "NA"
vcf_up$end[which(vcf_up$sv_type == "SNV")] = "NA"

vcf_up$sv_length = abs(vcf_up$sv_length)

sv_tmp_del = filter(vcf_up, sv_type == "del")
sv_tmp_del$end = sv_tmp_del$start + 1

sv_tmp_snv = filter(vcf_up, sv_type == "SNV")
sv_tmp_snv$end = sv_tmp_snv$start + 1
sv_tmp_snv$sv_length = 1

vcf_up = droplevels(vcf_up[!vcf_up$sv_type == "del",])
vcf_up = droplevels(vcf_up[!vcf_up$sv_type == "SNV",])
vcf_up = rbind(vcf_up, sv_tmp_del, sv_tmp_snv)
vcf_up = vcf_up[order(vcf_up$chr, vcf_up$start),]

vcf_up_large = subset(vcf_up, sv_length > 50)

#select harcoded genotypes
vcf_up = filter(vcf_up, genotype == c("1|1", "1|0", "0|1"))
vcf_up = droplevels(vcf_up)

vcf_up = subset(vcf_up, sv_length <= 50)

#subset and 1,2 and 2,1 genotypes to put on report as table
non_hardc_gen.table = filter(non_hardc_gen, genotype == c("1|2", "2|1"))
non_hardc_gen.table = droplevels(non_hardc_gen.table)

#tables
#create data tables from input data
type_tab = table1(~vcf_sub$sv_length | vcf_sub$sv_type, data = vcf_sub)
gen_tab = table1(~vcf_sub$sv_length | vcf_sub$genotype, data = vcf_sub)
large_tab = table1(~vcf_up_large$sv_length | vcf_up_large$sv_type, data = vcf_up_large)
nongen_tab = table1(~non_hardc_gen.table$sv_type | non_hardc_gen.table$genotype, data = non_hardc_gen.table)

#transform data tables to data frames
type_tab = as.data.frame(read_html(type_tab) %>% html_table(fill=TRUE))
gen_tab = as.data.frame(read_html(gen_tab) %>% html_table(fill=TRUE))
large_tab = as.data.frame(read_html(large_tab) %>% html_table(fill=TRUE))
nongen_tab = as.data.frame(read_html(nongen_tab) %>% html_table(fill=TRUE))

#remove first row
type_tab = type_tab[-1,]
gen_tab = gen_tab[-1,]
large_tab = large_tab[-1,]
nongen_tab = nongen_tab[-1,]

#rename varaible
names(type_tab)[names(type_tab) == "Var.1"] <- "MCT"
names(gen_tab)[names(gen_tab) == "Var.1"] <- "MCT"
names(large_tab)[names(large_tab) == "Var.1"] <- "MCT"
names(nongen_tab)[names(nongen_tab) == "Var.1"] <- "MCT"

#rbind columns from individual data frames
sum_metric_tab = cbind(gen_tab[1], gen_tab[2], gen_tab[3], type_tab[2], type_tab[3], type_tab[4], type_tab[5])

#format new line in data frame
ldel = filter(vcf_sub, sv_type == "del") %>% nrow()
lins = filter(vcf_sub, sv_type == "ins") %>% nrow()
lSNV = filter(vcf_sub, sv_type == "SNV") %>% nrow()
lall = nrow(vcf_sub)

lhet = filter(vcf_sub, genotype == "het") %>% nrow()
lhom = filter(vcf_sub, genotype == "hom") %>% nrow()
lgen = lhet + lhom

lldel = filter(vcf_up_large, sv_type == "del") %>% nrow()
llins = filter(vcf_up_large, sv_type == "ins") %>% nrow()
llall = nrow(vcf_up_large)

l1_2 = filter(non_hardc_gen, genotype == "1|2") %>% nrow()
l2_1 = filter(non_hardc_gen, genotype == "2|1") %>% nrow()
lnongen = nrow(non_hardc_gen)

#combine all lists
sum_metric_tab[nrow(sum_metric_tab) + 1,] = c("n", lhet, lhom, lins, ldel, lSNV, lall)
sum_metric_tab[nrow(sum_metric_tab) + 1,] = c("Type", "Heterozygous", "Homozygous", "Insertions", "Deletions", "SNVs", "All Variants")

type_tab[nrow(type_tab) + 1,] = c("n", lins, ldel, lSNV, lall)
type_tab[nrow(type_tab) + 1,] = c("Small-Variants Summary", "Insertions", "Deletions", "SNVs", "All Variants")

gen_tab[nrow(gen_tab) + 1,] = c("n", lhet, lhom, lgen)
gen_tab[nrow(gen_tab) + 1,] = c("Genotypes", "Heterozygous", "Homozygous", "Total")

large_tab[nrow(large_tab) + 1,] = c("n", llins, lldel, llall)
large_tab[nrow(large_tab) + 1,] = c("Variants > 50 bp", "Insertions", "Deletions", "All Variants")

nongen_tab[nrow(nongen_tab) + 1,] = c("n", l1_2, l2_1, lnongen)
nongen_tab[nrow(nongen_tab) + 1,] = c("Non-hardcoded Genotypes", "1|2", "2|1", "Total")

#rename variable names
sum_metric_tab$MCT[1] = "Length bp - Mean (SD)"
sum_metric_tab$MCT[2] = "Length bp - Median [Min, Max]"

type_tab$MCT[1] = "Mean (SD)"
type_tab$MCT[2] = "Median [Min, Max]"

gen_tab$MCT[1] = "Mean (SD)"
gen_tab$MCT[2] = "Median [Min, Max]"

large_tab$MCT[1] = "Mean (SD)"
large_tab$MCT[2] = "Median [Min, Max]"

#select rows
sum_metric_tab_filt = sum_metric_tab[c(4,3,1,2), (1:7)]
type_tab_filt = type_tab[c(4,3,1,2), (1:5)]
gen_tab_filt = gen_tab[c(4,3,1,2), (1:4)]
large_tab_filt = large_tab[c(4,3,1,2), (1:4)]
nongen_tab_filt = nongen_tab[c(5,4,1,2,3), (1:4)]

#transform summary table
sum_metric_tab_filt = t(sum_metric_tab_filt)
sum_metric_tab_filt = as.data.frame(sum_metric_tab_filt)
colnames(sum_metric_tab_filt) = as.character(unlist(sum_metric_tab_filt[1,]))
sum_metric_tab_filt = sum_metric_tab_filt[-1, ]

type_tab_filt = t(type_tab_filt)
type_tab_filt = as.data.frame(type_tab_filt)
colnames(type_tab_filt) = as.character(unlist(type_tab_filt[1,]))
type_tab_filt = type_tab_filt[-1, ]

gen_tab_filt = t(gen_tab_filt)
gen_tab_filt = as.data.frame(gen_tab_filt)
colnames(gen_tab_filt) = as.character(unlist(gen_tab_filt[1,]))
gen_tab_filt = gen_tab_filt[-1, ]

large_tab_filt = t(large_tab_filt)
large_tab_filt = as.data.frame(large_tab_filt)
colnames(large_tab_filt) = as.character(unlist(large_tab_filt[1,]))
large_tab_filt = large_tab_filt[-1, ]

nongen_tab_filt = t(nongen_tab_filt)
nongen_tab_filt = as.data.frame(nongen_tab_filt)
colnames(nongen_tab_filt) = as.character(unlist(nongen_tab_filt[1,]))
nongen_tab_filt = nongen_tab_filt[-1, ]

#set table theme
theme_1 = ttheme_default(core = list(fg_params = list(hjust = 0, x = 0.1, fontsize = 9)), colhead = list(fg_params = list(fontsize = 12, fontface = "bold")))

#convert data frame into grob
sum_metric_grob = tableGrob(sum_metric_tab_filt, theme = theme_1, rows = NULL)
type_grob = tableGrob(type_tab_filt, theme = theme_1, rows = NULL)
gen_grob = tableGrob(gen_tab_filt, theme = theme_1, rows = NULL)
large_grob = tableGrob(large_tab_filt, theme = theme_1, rows = NULL) 
nongen_grob = tableGrob(nongen_tab_filt, theme = theme_1, rows = NULL)

#write each data frame as a separate sheet to xlsx
list_of_datasets = list("Small Variants - Genotype" = gen_tab, "Small Variant - SV Type" = type_tab)
write.xlsx(list_of_datasets, paste0("out/small_variants/tables/", txtFileName, "_smallvariants_summary.xlsx"))

#export summary table as png
ggsave(sum_metric_grob, file = paste0("out/small_variants/figs/", txtFileName, "_smallvariants_summary.png"), limitsize = FALSE, width = 7, height = 2, units = c("in"), dpi = 300)
ggsave(type_grob, file = paste0("out/small_variants/figs/", txtFileName, "_smallvariants_type.png"), limitsize = FALSE, width = 7, height = 2, units = c("in"), dpi = 300)
ggsave(gen_grob, file = paste0("out/small_variants/figs/", txtFileName, "_smallvariants_genotypes.png"), limitsize = FALSE, width = 7, height = 2, units = c("in"), dpi = 300)
ggsave(large_grob, file = paste0("out/small_variants/figs/", txtFileName, "_smallvariants_large_variants.png"), limitsize = FALSE, width = 7, height = 2, units = c("in"), dpi = 300)
ggsave(nongen_grob, file = paste0("out/small_variants/figs/", txtFileName, "_smallvariants_nongen_variants.png"), limitsize = FALSE, width = 7, height = 2, units = c("in"), dpi = 300)

#export updated vcf format and vcf head
write.table(vcf_up, file = paste0("out/small_variants/tables/", txtFileName, "_smallvariants.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(vcfHead, file = paste0("out/small_variants/tables/", txtFileName, "_HEADER.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(non_hardc_gen, file = paste0("out/small_variants/tables/", txtFileName, "_non-hardcoded-genotypes_smallvariants.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(vcf_up_large, file = paste0("out/small_variants/tables/", txtFileName, "_none_smallvariants.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#print SV summary
vcf_tmp = vcf_sub %>% select(sv_length, sv_type, genotype)
variant_summary = summary(vcf_tmp)
cat ("\n        Summary Small Variants Metrics:")
kable(variant_summary)
cat ("\nVCF (small variants) analysed, tsv and BED files exported...\n")
