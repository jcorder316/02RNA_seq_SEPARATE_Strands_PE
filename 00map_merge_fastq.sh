#!/bin/bash

source  ~/.bash_profile

A=path/to/the/file_location/file1.fastq.gz
B=path/to/the/file_location/file2.fastq.gz


############## this line is to merge the fastq files from different runs
zcat $A $B > | gzip -c > $pa/${Sample}/merge.fq.gz


########## this line is to merge the files after mapping , bam files### you could use this one after you see that both runs give you similar results , at least in the 
##########  target genes

A1=path/to/the/file_location/file1.fastq.gz
B1=path/to/the/file_location/file2.fastq.gz

Eq=folder/to/save/output

bamtools merge -in $A 1-in $B1 -out $Eq/${Sample}/merge.bam

	samtools sort $Eq/${Sample}/merge.bam -o $Eq/${Sample}/merge_s.bam

		bamToBed -i $Eq/${Sample}/merge_s.bam > $Eq/${Sample}/merge_s.bed

			makeTagDirectory $Eq/${Sample}/merge_ucsc $Eq/${Sample}/merge_s.bed -single

				#makeUCSCfile $Eq/${Sample}/merge_ucsc -o $Eq/${Sample}/merge_ucsc_norm		

					samtools index $Eq/${Sample}/merge_s.bam

						bamCoverage -b $Eq/${Sample}/merge_s.bam -bs 20 --smoothLength 40 -p max  --normalizeUsing RPKM -e 150 -o $Eq/00dpi_${IP}_danRer11.bw
				
				
				
