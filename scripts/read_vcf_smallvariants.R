#R script for processing VCF for small variants.
#author, Carl-Adam Mattsson. omYXX informatics.

#load packages
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(knitr))
suppressMessages(library(gtsummary))
suppressMessages(library(gt))

#get sample name
vcf.list = list.files(path = "in/small_variants/", recursive = TRUE, pattern = "\\.vcf$", full.names = TRUE)
sample_name = gsub(".{4}$", '', vcf.list)
sample_name = substring(sample_name, 20)

#set parameters
now = format(Sys.time(), "%d_%m_%Y")
txtFileName <- paste0(sample_name, "_vcf_man_", now)

#read vcf into R (skip header)
vcf = read.table(file = paste0("in/small_variants/", sample_name, ".vcf"), sep = "\t", header = F, comment.char="#")

#subset vcf file
vcf_sum = vcf %>% select(V1, V2, V4, V5, V6, V10)

#format genotype field
vcf_sum$V10 = sub(":.+", "", vcf_sum$V10)

#define variable names
vcf_sum = setNames(vcf_sum, c("chr", "start", "ref", "alt", "qual", "genotype"))

#count nucleotides in ref and alt sequence (to get SV length)
vcf_sum$ref_count = nchar(vcf_sum$ref)
vcf_sum$alt_count = nchar(vcf_sum$alt)
vcf_sum$sv_lenght = vcf_sum$alt_count - vcf_sum$ref_count

#calculate variant end coordinates
vcf_sum$end = vcf_sum$start + vcf_sum$sv_lenght

#compute variant type (deletion/insertion)
vcf_sum$sv_type = vcf_sum$sv_lenght
vcf_sum$sv_type[which(vcf_sum$sv_lenght == "0")] = "SNV"
vcf_sum$sv_type[which(vcf_sum$sv_lenght < 0)] = "del"
vcf_sum$sv_type[which(vcf_sum$sv_lenght > 0)] = "dup"

#create new data frame with selected variables (prepare for plotting and exporting BED)
vcf_sub = vcf_sum %>% select(chr, start, end, sv_lenght, sv_type, genotype)

#transform genotype to factor with set levels
vcf_sub$genotype = as.factor(vcf_sub$genotype)

#transform haplotype information to het/hom
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="1|1"] = "hom"
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="1|0"] = "het"
levels(vcf_sub$genotype)[levels(vcf_sub$genotype)=="0|1"] = "het"

#replace end-coordinates for deletions with NA
vcf_sub$end[which(vcf_sub$sv_type == "del")] = "NA"

#transform into absolute values
vcf_sub$sv_lenght = abs(vcf_sub$sv_lenght)

#create tmp data frame for summary statistics
vcf_tmp1 = vcf_sub %>% select(sv_lenght, sv_type, genotype)
vcf_tmp2 = vcf_sub %>% select(sv_type, genotype)

#transform sv_length to factor with set levels
vcf_tmp1$sv_type = as.factor(vcf_tmp1$sv_type)

#export updated vcf format and BED
write.table(vcf_sub, file = paste0("out/small_variants/", txtFileName, "_smallvariants.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#print SV summary
variant_summary = summary(vcf_tmp1)
cat ("\n        Summary Small Variants Metrics:")
kable(variant_summary)
cat ("\nDone, VCF analysed and exported successfully! tsv and BED files are saved at /vcf_man/out/small_variants/\n\n")

#print summary table to file
summary_figure = vcf_tmp2 %>% tbl_summary() %>%
  modify_caption("**Table 1. Small Variants Summary Metrics**") %>%
  bold_labels() %>%
  as_gt() %>%
  gtsave(filename = paste0("out/small_variants/", txtFileName, "_smallvariants_summary.png"))