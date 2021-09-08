#R script for processing VCF for large and Structural Variants (SVs).
#author, Carl-Adam Mattsson. omYXX informatics.

#load packages
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(knitr))
suppressMessages(library(gtsummary))
suppressMessages(library(gt))

#get sample name
vcf.list = list.files(path = "in/SVs/", recursive = TRUE, pattern = "\\.vcf$", full.names = TRUE)
sample_name = gsub(".{4}$", '', vcf.list)
sample_name = substring(sample_name, 9)

#set parameters
now = format(Sys.time(), "%d_%m_%Y")
txtFileName <- paste0(sample_name, "_vcf_man_", now)

#read vcf into R (skip header)
vcf = read.table(file = paste0("in/SVs/", sample_name, ".vcf"), sep = "\t", header = F, comment.char="#")

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
vcf_sum$sv_lenght = vcf_sum$alt_count - vcf_sum$ref_count

#transform into absolute values
vcf_sum$sv_lenght = abs(vcf_sum$sv_lenght)

#calculate variant end coordinates
vcf_sum$end = vcf_sum$start + vcf_sum$sv_lenght

#create new data frame with selected variables (prepare for plotting and exporting BED)
vcf_sub = vcf_sum %>% select(chr, start, end, sv_lenght, sv_type, genotype)

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

#sort data frame
vcf_sub = vcf_sub %>% select(chr, start, end, sv_lenght, sv_type, sv_fam, genotype)
vcf_tmp = vcf_sub %>% select(sv_lenght, sv_type, sv_fam, genotype)

#export updated vcf format and BED
write.table(vcf_sub, file = paste0("out/SVs/", txtFileName, "_structuralvariants.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#print SV summary
variant_summary = summary(vcf_tmp)
cat ("\n        Summary Structural Variants Metrics:")
kable(variant_summary)
cat ("\nDone, VCF analysed and exported successfully! tsv and BED files are saved at /vcf_man/out/SVs/\n\n")

#print summary table to file
summary_figure = vcf_tmp %>% tbl_summary() %>%
  modify_caption("**Table 1. Structural Variants Summary Metrics**") %>%
  bold_labels() %>%
  as_gt() %>%
  gtsave(filename = paste0("out/SVs/", txtFileName, "_structuralvariants_summary.png"))