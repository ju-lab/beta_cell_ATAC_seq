#!/bin/bash
#I prepared a annotation table(biomart) containing TSS of mouse genes, which has ortholog human genes.
tail -n+2 processed/MIN6.bed > tmp.tmp
tail -n+2 biomart/mm10.tss.orth.inner.bed > tmp.tmp2
~/src/annotate_nearest_genes.sh -t tmp.tmp -f tmp.tmp2 -o processed/MIN6.nearest.bed -n 2
rm tmp.tmp tmp.tmp2
