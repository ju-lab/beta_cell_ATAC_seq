# Statistical analysis on proximity of orthology regiond of differential peak
# and DM-associated loci

suppressMessages(library(tidyverse))
library(stringr)

MIN6 <- read_tsv("processed2/MIN6_peak_annotated_20180411.bed",
				 col_types = "ciicdddciicccciicccciiciiiddiiil")

# basic strategy and statistics
nrow(MIN6) # number of MIN6 peaks
table(MIN6$difpeak) # differential peak
(mis <- apply(MIN6, 2, function(x) sum(is.na(x)))) # missing values
mis["mmnear_gid"] # mm genes with no near genes within >20kb, no orthology also
MIN6 %>% filter(is.na(.$mmnear_gid))

MIN6 %>% filter(!is.na(.$mmnear_gid)) %>%
	{!is.na(.$hgpeakalign_chr)} %>% table # peak mapped to somewhere in tad habouring orgology gene
MIN6 %>% filter(!is.na(.$mmnear_gid)) $ # peak mapped to somewhere in tad habouring orgology gene

MIN6 %>% filter(!is.na(.$mmnear_gid)) %>%
 {.$phastcons_mm_overlap > 0} %>% table # phastcon_mm overlap: mm peak is conserved region


MIN6 %>% filter(!is.na(.$mmnear_gid)) %>%
	transmute(phastcons_overlap = phastcons_mm_overlap > 0,
			  eval_cut = (!is.na(hgpeakalign_chr)) & hgpeakalign_eval <= 0.0001) %>%
	table() %>% gmodels::CrossTable() # 92% of blast results were phastcons-conserved region

MIN6 %>% filter(!is.na(mmnear_gid),
				phastcons_mm_overlap > 0,
				(!is.na(hgpeakalign_chr)) & hgpeakalign_eval <= 0.0001) %>%
	{print(nrow(.)); .} %>% # final orthology matched peaks
	.$difpeak %>% table  #and differential peaks

conserved_peaks <- MIN6 %>% filter(!is.na(mmnear_gid),
								   phastcons_mm_overlap > 0,
								   (!is.na(hgpeakalign_chr)) & hgpeakalign_eval <= 0.0001)
conserved_peaks %>%
	ggplot(aes(x = dmloci_count_within_500kb, fill = difpeak)) +
	geom_histogram(binwidth = 1)

# with(conserved_peaks,
# 	 boxplot(dmloci_count_within_500kb[difpeak],
# 	 		dmloci_count_within_500kb[!difpeak])) #severely ugly plot

with(conserved_peaks,table(dmloci_count_within_500kb, difpeak))

conserved_peaks %>% select(dmloci_count_within_500kb, difpeak) %>%
	rename(ndmloci = dmloci_count_within_500kb) %>%
	group_by(difpeak) %>%
	summarize(median = median(ndmloci),
			  mean = mean(ndmloci),
			  max = max(ndmloci),
			  min = min(ndmloci))

with(conserved_peaks,
	 t.test(dmloci_count_within_500kb[difpeak],
	 	    dmloci_count_within_500kb[!difpeak]))
# distance calculation to closest

