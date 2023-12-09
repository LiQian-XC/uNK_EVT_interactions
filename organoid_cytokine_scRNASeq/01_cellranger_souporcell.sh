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

##run souporcell shared samples
n1=($(seq 1 1 6))
n2=($(seq 28 1 35))
cq=(${n1[@]/#/6044STDY864056} ${n2[@]/#/Pla_Camb101239})
for cc in 1 2
do
	if [[ $cc -eq 1 ]]
	then
		ss=0
		dd1=4
		dd2=5
	else
		ss=6
		dd1=12
		dd2=13
	fi
	for i in `seq $ss 1 $dd1`
	do
		si=$((i + 1))
		for j in `seq $si 1 $dd2`
		do
			echo $i $j
			singularity exec /home/ql312/software/souporcell_latest.sif shared_samples.py -1 $dir/${cq[$i]}/souporcell_res -2 $dir/${cq[$j]}/souporcell_res -n 3 > ${cq[$i]}_${cq[$j]}_correspondence.txt
		done
	done	
done
