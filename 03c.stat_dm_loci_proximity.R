# Statistical analysis on proximity of orthology regiond of differential peak
# and DM-associated loci

suppressMessages(library(tidyverse))
library(stringr)

MIN6 <- read_tsv("processed3/MIN6_peak_annotated_20180412.bed",
				 col_types = "ciicdddciicccciicccciiciiiddiiil")

conserved_peaks <- MIN6 %>% filter(!is.na(mmnear_gid),
								   phastcons_mm_overlap > 0,
								   (!is.na(hgpeakalign_chr)) & hgpeakalign_eval <= 0.0001)

conserved_peaks %>% select(difpeak, dist_closest_dmloci) %>%
	rename(dist = dist_closest_dmloci) %>%
	arrange(dist) %>%
	group_by(difpeak) %>%
	summarize(min = min(dist), median = median(dist), mean = mean(dist),
			  max = max(dist))

conserved_peaks %>% select(difpeak, dist_closest_dmloci) %>%
	rename(dist = dist_closest_dmloci) %>%
	group_by(difpeak) %>%
	arrange(dist) %>%
	mutate(n = row_number(),
		   r = n/max(n)) -> plotdata
plotdata %>%
	ggplot(aes(x = dist, y = r, color = difpeak)) + geom_point()
plotdata %>%
	ggplot(aes(x = log10(dist), y = r, color = difpeak)) + geom_point()
plotdata %>%
	ggplot(aes(x = dist, y = r, color = difpeak)) + geom_point() + coord_cartesian(xlim = c(0,500000), ylim = c(0, 0.45))
library(manipulate)
manipulate(ggplot(plotdata, aes(x = dist, y = r, color = difpeak)) + geom_point() + coord_cartesian(xlim = c(0,15000 * XLIM), ylim = c(0, YLIM / 1000)),
		   XLIM = slider(0, 100, initial = 6),
		   YLIM = slider(0, 1000, initial = 93))
plotdata
library(svglite)
svglite( width = 7, height = 6, file = "Rplot_pvalue.svg")
par(fig = c(0,1,0,1))
plot(plotdata$dist/1000, plotdata$r, pch = 20, col = c("black", "red")[plotdata$difpeak + 1], cex = 0.5,xlim = c(0, 15*200),
	 xlab = "Distances from ATAC-seq peaks to the closest diabetes-assotiated loci(kbps)",
	 ylab = "Cumulative proportion", las = 1, xaxt = "n")
X <- 0:3 * 1000
axis(1, X, labels = prettyNum(X, big.mark = ","))
legend("topleft", bty = "n", legend = c("", "", "p-value* = 0.0003739"))
legend("topleft", bty = "n", legend = c("Non-differential peaks","Differential peaks (Prmt1-KO vs. WT)"), pch = 20, col = c("black", "red"))
rect(0, 0, 90, 93/1000, density = NULL, angle = 45, lty = 2)
par(fig = c(0.4, 0.95, 0.08, 0.7), new = TRUE)
plot(plotdata$dist/1000, plotdata$r, pch = 20, col = c("black", "red")[plotdata$difpeak + 1], cex = 0.5, xlim = c(0, 15*6), ylim = c(0, 93/1000),
	 xlab = "",
	 ylab = "", las = 1, xaxt = "n")
X <- 0:3 * 30
axis(1, X, labels = prettyNum(X, big.mark = ","))
par(fig = c(0,1,0,1))
dev.off()
 # ,

plotdata

difpeak_dist <- plotdata$dist[plotdata$difpeak == TRUE]
nondifpeak_dist <- plotdata$dist[plotdata$difpeak == FALSE]
A <- length(difpeak_dist)
B <- length(nondifpeak_dist)
C <- seq(1, B, length.out = A) %>% round
#C <- (1:A) * B %/% A
nondifpeak_dist_sampled <- nondifpeak_dist[C]
t.test(difpeak_dist, nondifpeak_dist)
nondifpeak_dist_sampled - difpeak_dist -> dist_dif
table(dist_dif[1:(A %/% 2)] > 0)
nondifpeak_dist_sampled[1:(A %/% 2)]
t.test(difpeak_dist, nondifpeak_dist_sampled, paired=TRUE)
t.test(difpeak_dist[1:(A %/% 2)], nondifpeak_dist_sampled[1:(A %/% 2)], paired=TRUE)

