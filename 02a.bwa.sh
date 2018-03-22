#!/bin/bash
cd /home/users/kjyi/Projects/temp_20180306

for i in `find -L . | grep "1.fastq.gz$" | grep -v RNA`; do
	sample_name=`echo $i | sed 's!.*/\([^/]*\)_1.fastq.gz$!\1!'`
	dir=`echo $i | sed 's!./fastq!./bam!;s!\(.*\)/[^/]*!\1!'`
	mkdir -p $dir
	run_bwa.sh \
		-s $sample_name \
		-o $dir \
		$i ${i/1.fastq.gz/2.fastq.gz}
done

exit 0
fastq
├── ATAC-seq -> /home/users/kjyi/Downloads/20180306/ATAC-seq
│   ├── ATAC-YFP-beta-cell
│   │   ├── ATAC-YFP-biKO-1_1.fastq.gz
│   │   ├── ATAC-YFP-biKO-1_2.fastq.gz
│   │   ├── ATAC-YFP-biKO-2_1.fastq.gz
│   │   ├── ATAC-YFP-biKO-2_2.fastq.gz
│   │   ├── ATAC-YFP-biKO-3_1.fastq.gz
│   │   ├── ATAC-YFP-biKO-3_2.fastq.gz
│   │   ├── ATAC-YFP-CON-1_1.fastq.gz
│   │   ├── ATAC-YFP-CON-1_2.fastq.gz
│   │   ├── ATAC-YFP-CON-2_1.fastq.gz
│   │   ├── ATAC-YFP-CON-2_2.fastq.gz
│   │   ├── ATAC-YFP-CON-3_1.fastq.gz
│   │   └── ATAC-YFP-CON-3_2.fastq.gz
│   └── MIN6
│       ├── ATAC_MIN6_1.fastq.gz
│       └── ATAC_MIN6_2.fastq.gz
├── H4R3me2a_ChIP-seq -> /home/users/kjyi/Downloads/20180306/H4R3me2a_ChIP-seq
│   ├── H4R3me2a_ChIP-1_1.fastq.gz
│   ├── H4R3me2a_ChIP-1_2.fastq.gz
│   ├── H4R3me2a_ChIP-2_1.fastq.gz
│   ├── H4R3me2a_ChIP-2_2.fastq.gz
│   ├── H4R3me2a_ChIP-3_1.fastq.gz
│   ├── H4R3me2a_ChIP-3_2.fastq.gz
│   ├── H4R3me2a_Input-1_1.fastq.gz
│   ├── H4R3me2a_Input-1_2.fastq.gz
│   ├── H4R3me2a_Input-2_1.fastq.gz
│   └── H4R3me2a_Input-2_2.fastq.gz
├── RNA-seq -> /home/users/kjyi/Downloads/20180306/RNA-seq
│   ├── Prmt1_biKO_12w
│   │   ├── Prmt1_biKO_12w_biko-1_1.fastq.gz
│   │   ├── Prmt1_biKO_12w_biko-1_2.fastq.gz
│   │   ├── Prmt1_biKO_12w_biko-2_1.fastq.gz
│   │   ├── Prmt1_biKO_12w_biko-2_2.fastq.gz
│   │   ├── Prmt1_biKO_12w_con-1_1.fastq.gz
│   │   ├── Prmt1_biKO_12w_con-1_2.fastq.gz
│   │   ├── Prmt1_biKO_12w_con-2_1.fastq.gz
│   │   └── Prmt1_biKO_12w_con-2_2.fastq.gz
│   ├── Prmt1_biKO_8w
│   │   ├── Prmt1_biKO_8w_biKO-1_1.fastq.gz
│   │   ├── Prmt1_biKO_8w_biKO-1_2.fastq.gz
│   │   ├── Prmt1_biKO_8w_biKO-2_1.fastq.gz
│   │   ├── Prmt1_biKO_8w_biKO-2_2.fastq.gz
│   │   ├── Prmt1_biKO_8w_con-1_1.fastq.gz
│   │   ├── Prmt1_biKO_8w_con-1_2.fastq.gz
│   │   ├── Prmt1_biKO_8w_con-2_1.fastq.gz
│   │   └── Prmt1_biKO_8w_con-2_2.fastq.gz
│   └── Prmt1_bKO_12w
│       ├── Prmt1_bKO_12w_bKO-1_1.fastq.gz
│       ├── Prmt1_bKO_12w_bKO-1_2.fastq.gz
│       ├── Prmt1_bKO_12w_bKO-2_1.fastq.gz
│       ├── Prmt1_bKO_12w_bKO-2_2.fastq.gz
│       ├── Prmt1_bKO_12w_con-1_1.fastq.gz
│       ├── Prmt1_bKO_12w_con-1_2.fastq.gz
│       ├── Prmt1_bKO_12w_con-2_1.fastq.gz
│       └── Prmt1_bKO_12w_con-2_2.fastq.gz
└── tree.sh -> /home/users/kjyi/tools/kjscript/tree.sh

8 directories, 49 files
