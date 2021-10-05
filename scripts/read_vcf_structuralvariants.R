#R script for processing VCF for large and Structural Variants (SVs).
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
vcf.list = list.files(path = "in/SVs/", recursive = TRUE, pattern = "\\.vcf$", full.names = TRUE)
sample_name = gsub(".{4}$", '', vcf.list)
sample_name = substring(sample_name, 9)

#set parameters
now = format(Sys.time(), "%d_%m_%Y")
txtFileName <- paste0(sample_name, "_vcf_man_", now)

#read vcf into R (skip header)
vcf = read.table(file = paste0("in/SVs/", sample_name, ".vcf"), sep = "\t", header = F, comment.char="#")

#drop chrX and Y
vcf$V1 = as.factor(vcf$V1)
vcf = filter(vcf, V1 == c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) 
vcf = droplevels(vcf)

#extract vcf head
vcfHead = readLines(paste0("in/SVs/", sample_name, ".vcf"))
vcfHead = vcfHead[grep('^##.*', vcfHead)]
vcfHead = as.data.frame(vcfHead)

#subset vcf file
vcf_sum = vcf %>% select(V1, V2, V4, V5, V6, V8, V10)

#format genotype field
vcf_sum$V10 = sub(":.+", "", vcf_sum$V10)

#subset INFO column and split by field
vcf_sum$V8 = sub(";.+","", vcf_sum$V8)
vcf_sum$V8 = sub("^.*\\=","", vcf_sum$V8)

#define variable names
vcf_sum = setNames(vcf_sum, c("chr", "start", "ref", "alt", "qual", "sv_type", "genotype"))

#count nucleotides in ref and alt sequence (to get SV length)
vcf_sum$ref_count = nchar(vcf_sum$ref)
vcf_sum$alt_count = nchar(vcf_sum$alt)
vcf_sum$sv_length = vcf_sum$alt_count - vcf_sum$ref_count

#transform into absolute values
vcf_sum$sv_length = abs(vcf_sum$sv_length)

#calculate variant end coordinates
vcf_sum$end = vcf_sum$start + vcf_sum$sv_length

#create new data frame with selected variables (prepare for plotting and exporting BED)
vcf_sub = vcf_sum %>% select(chr, start, end, sv_length, sv_type, genotype)

#transform sv_type and genotype to factor with set levels
vcf_sub$sv_type = as.factor(vcf_sub$sv_type)
vcf_sub$genotype = as.factor(vcf_sub$genotype)

#transform haplotype information to het/hom
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="1|1"] = "hom"
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="1|0"] = "het"
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="0|1"] = "het"

#replace end-coordinates for deletions with NA
vcf_sub$end[which(vcf_sub$sv_type == "SIMPLEDEL")] = "NA"
vcf_sub$end[which(vcf_sub$sv_type == "SUBSDEL")] = "NA"
vcf_sub$end[which(vcf_sub$sv_type == "CONTRAC")] = "NA"

#add new variable with SV-type "family name" (duplicate sv_type)
vcf_sub$sv_fam = vcf_sub$sv_type

#rename factors for SV family name
levels(vcf_sub$sv_fam)[levels(vcf_sub$sv_fam)=="CONTRAC"] = "del"
levels(vcf_sub$sv_fam)[levels(vcf_sub$sv_fam)=="SIMPLEDEL"] = "del"
levels(vcf_sub$sv_fam)[levels(vcf_sub$sv_fam)=="SUBSDEL"] = "del"
levels(vcf_sub$sv_fam)[levels(vcf_sub$sv_fam)=="DUP"] = "dup"
levels(vcf_sub$sv_fam)[levels(vcf_sub$sv_fam)=="SIMPLEINS"] = "ins"
levels(vcf_sub$sv_fam)[levels(vcf_sub$sv_fam)=="SUBSINS"] = "ins"
levels(vcf_sub$sv_fam)[levels(vcf_sub$sv_fam)=="ARTIFACT"] = "art"

#drop artifacts
vcf_sub = droplevels(vcf_sub[!vcf_sub$sv_fam == "art",])

#transform end coordinates for deletions (+1)
sv_tmp_del = filter(vcf_sub, sv_fam == "del")
sv_tmp_del$end = sv_tmp_del$start + 1

#remove deletions from main data frame and add in transformed end coordinates for deletions
vcf_sub = droplevels(vcf_sub[!vcf_sub$sv_fam == "del",])
vcf_sub = rbind(vcf_sub, sv_tmp_del)
vcf_sub = vcf_sub[order(vcf_sub$chr, vcf_sub$start),]

#sort data frame
vcf_sub = vcf_sub %>% select(chr, start, end, sv_length, sv_type, sv_fam, genotype)
vcf_tmp = vcf_sub %>% select(sv_length, sv_type, sv_fam, genotype)

#create new subset with variants smaller than 50bp
vcf_sum_small = subset(vcf_sum, sv_length <= 50)

#create variable for number of variants smaller or equal to 50 bp
n_smallvariants = nrow(vcf_sum_small)
cat ("\n        Number of variants ≤ 50bp:")
n_smallvariants

#filter out variants smaller than 50 bp
vcf_sum = subset(vcf_sum, sv_length > 50)
vcf_tmp = subset(vcf_tmp, sv_length > 50)

#create data tables from input data
type_tab = table1(~vcf_tmp$sv_length | vcf_tmp$sv_type, data = vcf_tmp)
fam_tab = table1(~vcf_tmp$sv_length | vcf_tmp$sv_fam, data = vcf_tmp)
gen_tab = table1(~vcf_tmp$sv_length | vcf_tmp$genotype, data = vcf_tmp)

#transform data tables to data frames
fam_tab = as.data.frame(read_html(fam_tab) %>% html_table(fill=TRUE))
gen_tab = as.data.frame(read_html(gen_tab) %>% html_table(fill=TRUE))
type_tab = as.data.frame(read_html(type_tab) %>% html_table(fill=TRUE))

#remove first row
fam_tab = fam_tab[-1,]
gen_tab = gen_tab[-1,]
type_tab = type_tab[-1,]

#rename varaible
names(fam_tab)[names(fam_tab) == "Var.1"] <- "MCT"
names(gen_tab)[names(gen_tab) == "Var.1"] <- "MCT"
names(type_tab)[names(type_tab) == "Var.1"] <- "MCT"

#rbind columns from individual data frames
sum_metric_tab = cbind(gen_tab[1], gen_tab[2], gen_tab[3], fam_tab[2], fam_tab[3], fam_tab[4], type_tab[7])

#format new line in data frame
lhet = filter(vcf_tmp, genotype == "het") %>% nrow()
lhom = filter(vcf_tmp, genotype == "hom") %>% nrow()
ldel = filter(vcf_tmp, sv_fam == "del") %>% nrow()
ldup = filter(vcf_tmp, sv_fam == "dup") %>% nrow()
lins = filter(vcf_tmp, sv_fam == "ins") %>% nrow()
lall = nrow(vcf_tmp)

#combine all lists
sum_metric_tab[nrow(sum_metric_tab) + 1,] = c("n", lhet, lhom, ldel, ldup, lins, lall)
sum_metric_tab[nrow(sum_metric_tab) + 1,] = c("Type", "Homozygous", "Heterozygous", "Deletions", "Duplications", "Insertions", "All Variants")

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
ggsave(sum_grob, file = paste0("out/SVs/figs/", txtFileName, "_SV_summary.png"), limitsize = FALSE, width = 7, height = 2, units = c("in"), dpi = 300)

#write each data frame as a separate sheet to xlsx
list_of_datasets <- list("Structural Variants - Genotype" = gen_tab, "Structural Variants - SV-Type" = fam_tab, "Structural Variant - SV-SubType" = type_tab)
write.xlsx(list_of_datasets, paste0("out/SVs/tables/", txtFileName, "_structuralvariants_summary.xlsx"))

#revert vcf to standard formatted genotype
vcf_up = vcf_sum %>% select(chr, start, end, sv_length, sv_type, genotype)
vcf_up$sv_type = as.factor(vcf_up$sv_type)
vcf_up$genotype = as.factor(vcf_up$genotype)

vcf_up$end[which(vcf_up$sv_type == "SIMPLEDEL")] = "NA"
vcf_up$end[which(vcf_up$sv_type == "SUBSDEL")] = "NA"
vcf_up$end[which(vcf_up$sv_type == "CONTRAC")] = "NA"

vcf_up$sv_fam = vcf_up$sv_type

levels(vcf_up$sv_fam)[levels(vcf_up$sv_fam)=="CONTRAC"] = "del"
levels(vcf_up$sv_fam)[levels(vcf_up$sv_fam)=="SIMPLEDEL"] = "del"
levels(vcf_up$sv_fam)[levels(vcf_up$sv_fam)=="SUBSDEL"] = "del"
levels(vcf_up$sv_fam)[levels(vcf_up$sv_fam)=="DUP"] = "dup"
levels(vcf_up$sv_fam)[levels(vcf_up$sv_fam)=="SUBSINS"] = "ins"
levels(vcf_up$sv_fam)[levels(vcf_up$sv_fam)=="SIMPLEINS"] = "ins"
levels(vcf_up$sv_fam)[levels(vcf_up$sv_fam)=="ARTIFACT"] = "art"

#drop artifacts
vcf_up = droplevels(vcf_up[!vcf_up$sv_fam == "art",])

sv_tmp_del_up = filter(vcf_up, sv_fam == "del")
sv_tmp_del_up$end = sv_tmp_del_up$start + 1
vcf_up = droplevels(vcf_up[!vcf_up$sv_fam == "del",])
vcf_up = rbind(vcf_up, sv_tmp_del_up)

vcf_up = vcf_up[order(vcf_up$chr, vcf_up$start),]
vcf_up = vcf_up %>% select(chr, start, end, sv_length, sv_type, sv_fam, genotype)

#transform sv_type and genotype to factor with set levels
vcf_sum_small$sv_type = as.factor(vcf_sum_small$sv_type)
vcf_sum_small$genotype = as.factor(vcf_sum_small$genotype)

#replace end-coordinates for deletions with NA
vcf_sum_small$end[which(vcf_sum_small$sv_type == "SIMPLEDEL")] = "NA"
vcf_sum_small$end[which(vcf_sum_small$sv_type == "SUBSDEL")] = "NA"
vcf_sum_small$end[which(vcf_sum_small$sv_type == "CONTRAC")] = "NA"

#add new variable with SV-type "family name" (duplicate sv_type)
vcf_sum_small$sv_fam = vcf_sum_small$sv_type

#rename factors for SV family name
levels(vcf_sum_small$sv_fam)[levels(vcf_sum_small$sv_fam)=="CONTRAC"] = "del"
levels(vcf_sum_small$sv_fam)[levels(vcf_sum_small$sv_fam)=="SIMPLEDEL"] = "del"
levels(vcf_sum_small$sv_fam)[levels(vcf_sum_small$sv_fam)=="SUBSDEL"] = "del"
levels(vcf_sum_small$sv_fam)[levels(vcf_sum_small$sv_fam)=="DUP"] = "dup"
levels(vcf_sum_small$sv_fam)[levels(vcf_sum_small$sv_fam)=="SIMPLEINS"] = "ins"
levels(vcf_sum_small$sv_fam)[levels(vcf_sum_small$sv_fam)=="SUBSINS"] = "ins"
levels(vcf_sum_small$sv_fam)[levels(vcf_sum_small$sv_fam)=="ARTIFACT"] = "art"

#drop artifacts
vcf_sum_small = droplevels(vcf_sum_small[!vcf_sum_small$sv_fam == "art",])

#transform end coordinates for deletions (+1)
sv_tmp_del = filter(vcf_sum_small, sv_fam == "del")
sv_tmp_del$end = sv_tmp_del$start + 1

#remove deletions from main data frame and add in transformed end coordinates for deletions
vcf_sum_small = droplevels(vcf_sum_small[!vcf_sum_small$sv_fam == "del",])
vcf_sum_small = rbind(vcf_sum_small, sv_tmp_del)
vcf_sum_small = vcf_sum_small[order(vcf_sum_small$chr, vcf_sum_small$start),]

#sort data frame
vcf_sum_small = vcf_sum_small %>% select(chr, start, end, sv_length, sv_type, sv_fam, genotype)

#export updated vcf format and vcf head
write.table(vcf_up, file = paste0("out/SVs/tables/", txtFileName, "_structuralvariants.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(vcf_sum_small, file = paste0("out/SVs/tables/", txtFileName, "_none_structuralvariants.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(vcfHead, file = paste0("out/SVs/tables/", txtFileName, "_HEADER.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#print SV summary
variant_summary = summary(vcf_tmp)
cat ("\n        Summary Structural Variants Metrics:")
kable(variant_summary)
cat ("\nVCF (SVs) analysed, tsv and BED files exported...\n")