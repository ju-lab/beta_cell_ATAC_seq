#!/bin/bash
for i in `find -L fastq/ATAC-seq | grep "1.fastq.gz$" | grep trim_ `; do
	sample=`echo $i | sed 's!.*/\([^/]*\)_1.fastq.gz$!\1!;s!trim_!!'`
	dir=`echo $i | sed 's!fastq!bowtie2!;s!\(.*\)/[^/]*!\1!;s!trim_!!'`
	mkdir -p $dir
	~/src/atac/run_atac_pipe.sh \
		--output_bam $dir \
		--process spr \
		$sample \
		$i \
		${i/1.fastq.gz/2.fastq.gz} 
done
