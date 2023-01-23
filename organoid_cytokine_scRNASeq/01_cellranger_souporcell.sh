#! /bin/bash

cc=$1
echo $cc
dir=/home/ql312/rds/rds-turco-lab2-jxGoj1xLQV4/trophoblast_organoid_from_Sanger

if [[ $cc = 6044* ]]
then
	sdir=32734
else
	sdir=37295
fi

##run cellranger
mkdir $dir/$cc/cellranger_res/
cd $dir/$cc/cellranger_res/
cellranger count --id=$cc --transcriptome=/home/ql312/rds/hpc-work/fluffy/run_cellranger/refdata-gex-GRCh38-2020-A --fastqs=$dir/$cc/$sdir --sample=$cc --localcores=12 1> log_$cc.out 2> log_$cc.err

##run souporcell
gzip -cd $dir/$cc/cellranger_res/$cc/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $dir/$cc/cellranger_res/barcodes.tsv 
mkdir $dir/$cc/souporcell_res
singularity exec /home/ql312/software/souporcell_latest.sif souporcell_pipeline.py -i $dir/$cc/cellranger_res/$cc/outs/possorted_genome_bam.bam -b $dir/$cc/cellranger_res/barcodes.tsv -f ../fluffy/run_cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa -t 12 -o $dir/$cc/souporcell_res -k 3 --common_variants ../fluffy/analyses/common_variants_grch38_modify.vcf 
