#!/bin/bash

for i in `find -L . | grep 1.fastq.gz | grep RNA`; do
	sample_name=`echo $i | sed 's!.*/\([^/]*\)_1.fastq.gz$!\1!'`
	outdir_star=`echo $i | sed 's!./fastq/RNA-seq!./star!;s!\(.*\)/[^/]*!\1!'` 
	outdir_rsem=`echo $i | sed 's!./fastq/RNA-seq!./rsem!;s!\(.*\)/[^/]*!\1!'` 
	run_star.sh \
		-s $sample_name \
		--outdir_star $outdir_star \
		--outdir_rsem $outdir_rsem \
		--reference mm10 \
		$i ${i/1.fastq.gz/2.fastq.gz} \
		--process star,rsem
done
