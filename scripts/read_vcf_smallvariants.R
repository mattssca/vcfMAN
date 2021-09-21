#R script for processing VCF for small variants.
#author, Carl-Adam Mattsson. omYXX informatics.

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
vcf_sum$sv_type[which(vcf_sum$sv_length > 0)] = "dup"

#create new data frame with selected variables (prepare for plotting and exporting BED)
vcf_sub = vcf_sum %>% select(chr, start, end, sv_length, sv_type, genotype)

#transform genotype to factor with set levels
vcf_sub$genotype = as.factor(vcf_sub$genotype)
vcf_sub$sv_type = as.factor(vcf_sub$sv_type)
vcf_sub$chr = as.factor(vcf_sub$chr)

#subset on autosomes
selected = paste0("chr", c(1:22, "X"))
#selected = paste0("chr", c(1:22, "x,y"))
vcf_sub = vcf_sub[vcf_sub$chr %in% selected,]

#transform haplotype information to het/hom
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="1|1"] = "hom"
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="1|0"] = "het"
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="0|1"] = "het"
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="1|2"] = "het"
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="2|1"] = "het"

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

#create data tables from input data
type_tab = table1(~vcf_sub$sv_length | vcf_sub$sv_type, data = vcf_sub)
gen_tab = table1(~vcf_sub$sv_length | vcf_sub$genotype, data = vcf_sub)

#transform data tables to data frames
gen_tab = as.data.frame(read_html(gen_tab) %>% html_table(fill=TRUE))
type_tab = as.data.frame(read_html(type_tab) %>% html_table(fill=TRUE))

#remove first row
gen_tab = gen_tab[-1,]
type_tab = type_tab[-1,]

#rename varaible
names(gen_tab)[names(gen_tab) == "Var.1"] <- "MCT"
names(type_tab)[names(type_tab) == "Var.1"] <- "MCT"

#rbind columns from individual data frames
sum_metric_tab = cbind(gen_tab[1], gen_tab[2], gen_tab[3], type_tab[2], type_tab[3], type_tab[4], type_tab[5])

#format new line in data frame
lhet = filter(vcf_sub, genotype == "het") %>% nrow()
lhom = filter(vcf_sub, genotype == "hom") %>% nrow()
ldel = filter(vcf_sub, sv_type == "del") %>% nrow()
ldup = filter(vcf_sub, sv_type == "dup") %>% nrow()
lSNV = filter(vcf_sub, sv_type == "SNV") %>% nrow()
lall = nrow(vcf_sub)

#combine all lists
sum_metric_tab[nrow(sum_metric_tab) + 1,] = c("n", lhet, lhom, ldel, ldup, lSNV, lall)
sum_metric_tab[nrow(sum_metric_tab) + 1,] = c("Type", "Homozygous", "Heterozygous", "Duplications", "Deletions", "SNVs", "All Variants")

#rename variable names
sum_metric_tab$MCT[1] = "Length bp - Mean (SD)"
sum_metric_tab$MCT[2] = "Length bp - Median [Min, Max]"

#select rows
sum_metric_tab_filt = sum_metric_tab[c(4,3,1,2), (1:7)]

#transform summary table
sum_metric_tab_filt = t(sum_metric_tab_filt)
sum_metric_tab_filt = as.data.frame(sum_metric_tab_filt)
colnames(sum_metric_tab_filt) = as.character(unlist(sum_metric_tab_filt[1,]))
sum_metric_tab_filt = sum_metric_tab_filt[-1, ]

#convert data frame into grob
sum_grob = tableGrob(sum_metric_tab_filt, rows = NULL)

#export summary table as png
ggsave(sum_grob, file = paste0("out/small_variants/figs/", txtFileName, "_smallvariants_summary.png"), limitsize = FALSE, width = 7, height = 2, units = c("in"), dpi = 300)

#write each data frame as a separate sheet to xlsx
list_of_datasets <- list("Small Variants - Genotype" = gen_tab, "Small Variant - SV Type" = type_tab)
write.xlsx(list_of_datasets, paste0("out/small_variants/tables/", txtFileName, "_smallvariants_summary.xlsx"))

#export updated vcf format and vcf head
write.table(vcf_sub, file = paste0("out/small_variants/tables/", txtFileName, "_smallvariants.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(vcfHead, file = paste0("out/small_variants/tables/", txtFileName, "_HEADER.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#print SV summary
vcf_tmp = vcf_sub %>% select(sv_length, sv_type, genotype)
variant_summary = summary(vcf_tmp)
cat ("\n        Summary Small Variants Metrics:")
kable(variant_summary)
cat ("\n        :")
type_tab
cat ("\n        :")
gen_tab
cat ("\n        :")
cat ("\nVCF analysed, tsv and BED files exported...\n")
