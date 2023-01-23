#! /bin/bash

for file in `ls ../raw_data/bulk_RNAseq_KIR_C1_C2/*gz`
do
	cc=${file##*/}
	qq=${cc%%_R1_merged.fastq.gz}
	STAR --runMode alignReads --runThreadN 8 --genomeDir /home/ql312/rds/hpc-work/STAR_index/GRCh38_release37_PRI --readFilesIn /home/ql312/rds/hpc-work/raw_data/bulk_RNAseq_KIR_C1_C2/$cc --readFilesCommand zcat --outFilterMultimapNmax 5 --outFilterMismatchNmax 1 --outFileNamePrefix /home/ql312/rds/hpc-work/KIR_C1_C2/mapping_res/$qq --outSAMtype BAM SortedByCoordinate --outFilterMatchNminOverLread 0.95 --outFilterIntronMotifs RemoveNoncanonicalUnannotated
done
