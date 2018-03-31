#!/bin/bash
# Anotate MIN6 peak ~ (1) nearest mouse gene, (2) ortholog human gene, (3) TAD


# Prepare annotation data
# basic data downloaded from biomart - mouse TSS, human TSS
# http://grch37.ensembl.org/biomart/martview/b241fe1c5a3f3ea94be0edfa70765082?VIRTUALSCHEMANAME=default&ATTRIBUTES=mmusculus_gene_ensembl.default.feature_page.chromosome_name|mmusculus_gene_ensembl.default.feature_page.transcription_start_site|mmusculus_gene_ensembl.default.feature_page.ensembl_gene_id|mmusculus_gene_ensembl.default.feature_page.ensembl_transcript_id|mmusculus_gene_ensembl.default.feature_page.external_gene_name|hsapiens_gene_ensembl.default.feature_page.chromosome_name|hsapiens_gene_ensembl.default.feature_page.transcription_start_site|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.external_gene_name&FILTERS=&VISIBLEPANEL=resultspanel
pushd processed
input=./biomart_mmhg.tsv
annot=./biomart_annotation_ready.tsv
min6=./MIN6_peak.bed
output=./MIN6_peak_annotated.bed
if false ; then #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# change order of column (to human first) and remove header
awk -F '\t' 'BEGIN {OFS=FS} NR>1 {print $6,$7-1,$7,$8,$9,$10,$1,$2-1,$2,$3,$4,$5}' $input |\
	sed '/^GL00/d; /^Y/d; /^MT/d;' > $input.tmp1

# sampling one rows per each mm gene
R --slave << EOF
library(tidyverse)
set.seed(42)
read_tsv("$input.tmp1", col_types = cols(.default = "c"), col_names = F) %>% 
	group_by(X4) %>%
	sample_n(1) %>% write_tsv("$input.tmp2", col_names = F)
EOF
wc -l $input.tmp* # 732579 ~> 18225
# Panc1.TAD downloaded, coord as hg19
#
tad=./original/Panc1.TAD.bed # no Y chr
cut -f1,2,3 $tad | sed 's/^chrM/chrMT/;s/^chr//' > $tad.cut
#annotate_nearest_feature.sh -t $input.tmp3 -f ${tad/bed/cut.bed} -o $input.tmp4 -N
sort-bed $input.tmp2 > $input.tmp3
sort-bed $tad.cut > $tad.sort
closest-features --closest --dist --delim $'\t' $input.tmp3 $tad.sort > $input.tmp4 
fi #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Replace TADs out of range of TSS to +- 500kb
# Reorder of column of bed to mouse-first order again
R --slave << EOF
library(tidyverse)
tmp <- read_tsv("$input.tmp4", col_types = "ciicccciicccciii", col_names = F)
tmp\$X13[tmp\$X16 > 0] <- tmp\$X1[tmp\$X16 > 0]
tmp\$X14[tmp\$X16 > 0] <- pmax(0,tmp\$X2[tmp\$X16 > 0] - 500000) %>% as.integer
tmp\$X15[tmp\$X16 > 0] <- tmp\$X14[tmp\$X16 > 0] + 1000000 %>% as.integer
tmp %>% select(-X16) %>% 
	select(X7:X12, X1:X6, everything()) %>% 
	write_tsv("$annot", col_names = F)
EOF
sed -i '/^GL/d; /^JH584/d;/^Y/d' $annot
# Annotate MIN6 peak
sed '/^1_GL4/d;/^4_/d;/^MT/d;/^Un_/d;/^X_/d;/^Y/d;' $min6 > $min6.tmp
annotate_nearest_feature.sh -t $min6.tmp -f $annot -o $output -N
exit 0
