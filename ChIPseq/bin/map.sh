#!/bin/bash

export PYTHONPATH=$HOME/lib/python3.7/site-packages/:$PYTHONPATH
export PATH=$HOME/bin/:$PATH

GENOMESIZE=/mnt/projects/rpd/genomes/mm10/mm10.fa.fai

F1=$1
F2=$2
SAMPLE=$3


TF1=$SAMPLE"_1_val_1.fq.gz"
TF2=$SAMPLE"_2_val_2.fq.gz"

INDEX=/mnt/projects/rpd/genomes/mm10/bwa_path/nucleotide/mm10.fa

mkdir $SAMPLE
cd $SAMPLE

trim_galore --fastqc --gzip --length 100 --paired $F1 $F2
bwa mem -t 4 $INDEX $TF1 $TF2 > alignment.sam
/mnt/software/bin/samtools view -bS alignment.sam | /mnt/software/bin/samtools sort - sorted_alignment
/mnt/software/bin/samtools flagstat sorted_alignment.bam > statistcis.txt
/mnt/software/bin/samtools rmdup sorted_alignment.bam rmdup_sorted_alignment.bam
/mnt/software/bin/samtools view -h -F 4 -b sorted_alignment.bam > sorted_alignment_only_mapped.bam
/mnt/software/bin/bedtools bamtobed -i sorted_alignment_only_mapped.bam > sorted_alignment_only_mapped.bed

qsub -pe OpenMP 1 -l h_rt=5:00:00,mem_free=30G -V -cwd -b y samtools sort -n sorted_alignment_only_mapped.bam sorted_alignment_byname

run_dfilter.sh -d=sorted_alignment_only_mapped.bam -o=ATACseq.bed -f=bam -ks=50 -bs=100 -lpval=2 -wig


/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/macs2/MACS2-2.2.6/bin/macs2  callpeak -t pcr_removed_sorted_alignment.bam -f BAM -g 1.89e9 -n Atac_PCR --outdir Atac -B
cd Atac_PCR
bedGraphToBigWig Atac_treat_pileup.bdg $GENOMESIZE $SAMPLE".atac.bw"
