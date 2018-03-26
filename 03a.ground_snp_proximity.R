#!/home/users/kjyi/bin/Rscript
#'
#'    Ground peak-gwas_snp proximity calculation
#'
#' Herein ground peak-gwas_snp proximity probability will be meaasured with
#' given randomized rules in each step of filtering. (e.g. Counts of
#' peak-proxmial DM.GWAS_SNP in differential peak vs. non-differential peak)
#'
#' The ground snp density will be calculated in following aspects
#' 1. differential ATAC peak vs.non-differential peak (MIN6 vs KO)
#' 2. irrelavent GWAS keywords or random loci in SNP data. (count matched)
#'
#' Additional filtering will be considered
#' 1. differentially expressed / commonly expressed (MIN6 vs KO, CON vs KO)
#' 2.

# load packages
library(tidyverse)

# Read data
peak_MIN6 <- read_tsv("processed/MIN6.bed")
peak_dif <- read_tsv("processed/depleted.bed")

gwas <- read_tsv("processed/gwas.position", col_types="ici")
#http://www.ebi.ac.uk/gwas/api/search/downloads/alternative
# rownumbers and coord only (awk '{printf "%d\t%s\n", NR, $0}' - > gwas.full ,  awk -F $'\t' 'BEGIN {OFS=FS} {print $1,$13,$14}' gwas.full > gwas.position)

gwas_search <- function(q){
	# full table also have same row number.
	i <- system(paste0("grep -E '", q, "' processed/gwas.full | cut -f1"), intern = T)
	filter(gwas, ROW %in% i) %>% select(-ROW)
}
gwas_search("diabetes")

gwas_search("iabete")
