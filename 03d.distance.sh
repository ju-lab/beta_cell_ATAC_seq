#!/bin/bash
cd processed3
min6=MIN6_peak_annotated_20180410.bed
min6_b=${min6/.bed/.bound.bed} # mapped region median +- 250kb boundary
min6_c=${min6/.bed/.cons.bed} # annotate phastcons overlap
min6_h=${min6/.bed/.bound_first.bed} # column start with boundary bed format
min6_hc=${min6/.bed/.count_annotated.bed} # count_annotated
min6_d=MIN6_peak_annotated_20180412.bed # differential peak annotated


cons_orig=phastConsElements60wayEuarchontoGlires.txt
conserved=conserved.mm.bed
conserved_ann_metrix=conserved_annotation_result_table.txt

gwas_orig=gwas_catalog_v1.0.1.tsv
gwas_dm=gwas_dm.tsv
gwas_dm_hg19=gwas_dm_hg19.tsv

diff=ATAC-seq_depleted_region.bed
difid=depleted_peak_id.txt
# trim chr name of phastcons file
cut -f 2,3,4 $cons_orig | sed 's/^chrM/chrMT/;s/^chr//' > $conserved

# trim chr name of mapped human coord in min6 peak file
# then, add range info +-250kb around mapped human region
R --slave << EOF 
suppressMessages(library(tidyverse))
library(stringr)
read_tsv("$min6", col_names = F, col_types = "ciicdddciicccciicccciiciiiddi",
		 na = c("", "NA", "x")) %>%
		 mutate(mapped_area_chr = str_replace(X23, "chr", "") %>% str_replace("M", "MT"), 
		 mapped_area_lower = pmax(0, (X24+X25) %/% 2 - 1) %>% as.integer(), #prevent negative coord
		 mapped_area_upper = (X24+X25) %/% 2 %>% as.integer())  %>%
		 write_tsv("$min6_b", col_names = F)
EOF
# calculalte intersection between min6 peak (mm10) and phastcon conserved region
bedtools intersect -a $min6_b -b $conserved -c > $min6_c
R --slave << EOF > $conserved_ann_metrix
suppressMessages(library(tidyverse))
library(stringr)
Data <- read_tsv("$min6_c", col_names = F, col_types = "ciicdddciicccciicccciiciiiddiccci")
table(mapped = ! is.na(Data\$X30), conserved = Data\$X33 > 0)
Data %>% select(X30:X32, everything()) %>% write_tsv("$min6_h", col_names=F)
EOF

# prepare DM associated loci
grep -e '[Dd]iabetes' $gwas_orig | awk -F$'\t' 'BEGIN {OFS=FS}{print $12,$13-1, $13}' > $gwas_dm

liftover=~/tools/liftover/liftOver
chain=~/tools/liftover/hg38ToHg19.over.chain

sed -i 's/^MT/M/;s/^/chr/;/-1/d;/;/d' $gwas_dm
$liftover $gwas_dm $chain $gwas_dm_hg19 /dev/null &> /dev/null
sed -i 's/^chrM/chrMT/;s/^chr//' $gwas_dm
sed -i 's/^chrM/chrMT/;s/^chr//' $gwas_dm_hg19
#wc -l $gwas_orig $gwas_dm $gwas_dm_hg19 > gwas_loci_count.txt;cat gwas_loci_count.txt


# count gwas loci in human ortholog region
sed 's/^NA\tNA\tNA/NA\t0\t0/' $min6_h | sort -k1,1 -k2,2n > ${min6_h/bed/sort.bed}
sort -k1,1 -k2,2n $gwas_dm_hg19 > ${gwas_dm_hg19/tsv/sort.tsv}
bedtools closest -a ${min6_h/bed/sort.bed} -b ${gwas_dm_hg19/tsv/sort.tsv} -d -t first | cut -f-33,37 > $min6_hc

# annotate is the peaks differential peak
tail -n+2 $diff | cut -f4 > $difid

R --slave << EOF > /dev/null
suppressMessages(library(tidyverse))
library(stringr)
diffID <- read_tsv("$difid", col_names = F, col_types = "c")\$X1
mydata <- read_tsv("$min6_hc", col_names = F, col_types = cols(.default="c")) %>%
	select(-X1,-X2,-X3) %>%
	mutate(diff = X7 %in% diffID)
colnames(mydata) <- c("mmdifp_chr", "mmdifp_start", "mmdifp_end", "peakid", "fc",
	"minuslogp", "minuslogq", "mmneartss_chr", "mmneartss_start", "mmneartss_end",
	"mmnear_gid", "mmnear_tid", "mmnear_gname", "hgtss_chr", "hgtss_start", "hgtss_end",
	"hg_gid", "hg_tid", "hg_gname", "hgtad_chr", "hgtad_start", "hgtad_end",
	"hgpeakalign_chr", "hgpeakalign_start", "hgpeakalign_end", "hgpeakalign_length",
	"hgpeakalign_eval", "hgpeakalign_match", "hgpeakalign_NmatchInTad",
	"phastcons_mm_overlap", "dist_closest_dmloci", "difpeak")
mydata %>%
	write_tsv("$min6_d")
EOF
rm $min6_b $min6_c $min6_h $min6_hc $conserved $gwas_dm $gwas_dm_hg19 $difid

# Example output importing R cmd
#> MIN6 <- readr::read_tsv("processed2/MIN6_peak_annotated_20180411.bed",
#+				 col_types = "ciicdddciicccciicccciiciiiddiiil")


