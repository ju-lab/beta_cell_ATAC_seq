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

# load packages ----------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(tximport) # for RNA count matrix load
library(stringr)

# Read differential peak data --------------------------------------------------
peak_MIN6 <- read_tsv("processed/MIN6.bed")         # 0-based
peak_dif <- read_tsv("processed/depleted.bed")      # 0-based

# Biomart gene orthology data --------------------------------------------------
biomart <- read_tsv("biomart/mm10hg38_transcript.tsv") # 1-based

# Annotate nearest genes
# done by 03b.annotate_nearest_gene.sh
peak_MIN6 <-
peak_dif





gwas <- read_tsv("processed/gwas.position", col_types = "ici")   # maybe 1-based
#http://www.ebi.ac.uk/gwas/api/search/downloads/alternative
# rownumbers and coord only (awk '{printf "%d\t%s\n", NR, $0}' - > gwas.full ,  awk -F $'\t' 'BEGIN {OFS=FS} {print $1,$13,$14}' gwas.full > gwas.position)


# gwas search function ---------------------------------------------------------
gwas_search <- function(q){
	# full table also have same row number.
	i <- system(paste0("grep -E '", q, "' processed/gwas.full | cut -f1"),
				intern = T)
	filter(gwas, ROW %in% i) %>% select(-ROW)
}
# gwas_search("diabetes")
# gwas_search("iabete")

# RNA sequencing DEGs ----------------------------------------------------------
load_rsem <- function(type = "gene") {
	if (type == "gene") {
		pattern <- "ls rsem/*/*.genes.results"
		txIdCol <- FALSE
	} else {
		pattern <- "ls rsem/*/*.isoforms.results"
		txIdCol <- TRUE
	}
	files <- system(intern = T, pattern)
	names(files) <- str_extract(files, regex("(?<=/)[^/]*(?=.rsem)"))
	txi <- tximport(files,  type = "rsem", txIdCol = txIdCol)
	col_data <- data.frame(condition = str_replace(colnames(txi$counts), regex("-.$"), ""))
	rownames(col_data) <- colnames(txi$counts)
	txi$length[txi$length == 0] <- 1
	DESeqDataSetFromTximport(txi, col_data, ~ condition)
}

dds_g <- load_rsem()
dds_i <- load_rsem("isoform")
# ...defer...


