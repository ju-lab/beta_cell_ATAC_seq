#!/bin/bash
~/src/atac/run_atac_pipe.sh \
	--thread 22 \
	--memory 40GB \
	--output_bam test1a \
	--output_qc test1a \
	--log test1a \
	--process bowtie2,postalign \
	test01a \
	./fastq/ATAC-seq/MIN6/ATAC_MIN6_1.fastq.gz \
	./fastq/ATAC-seq/MIN6/ATAC_MIN6_2.fastq.gz

~/src/atac/run_atac_pipe.sh \
	--thread 22 \
	--memory 40GB \
	--output_bam test1b \
	--output_qc test1b \
	--log test1b \
	--process bowtie2,postalign \
	test01b \
	./fastq/ATAC-seq/MIN6/ATAC_MIN6_1.fastq.gz 
