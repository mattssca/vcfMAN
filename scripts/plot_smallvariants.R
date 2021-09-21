#load packages
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggthemr))
suppressMessages(library(gridExtra))
suppressMessages(library(RCircos))

HG002_GRCh38_CMRG_smallvar_v1.00_vcf_man_21_09_2021_smallvariants.txt

#get sample name
vcf.list = list.files(path = "in/small_variants/", recursive = TRUE, pattern = "\\.vcf$", full.names = TRUE)
sample_name = gsub(".{4}$", '', vcf.list)
sample_name = substring(sample_name, 20)
now = format(Sys.time(), "%d_%m_%Y")

#set parameters
txtFileName <- paste0(sample_name, "_vcf_man_", now)

#read SV callset into R data frame (skips first Header line)
sv_calls = read.table(file = paste0("out/small_variants/tables/", txtFileName, "_smallvariants.txt"), header = TRUE, sep = "\t")