#### RNA-seq / single-cell RNA-seq (except 10x) analysis ####
#### Author : Wilson Tan
#### Date   : 1 October 2020
#### Version: 1
#############################################################

#!/bin/bash

## Usage
## sh rnaseq_align.sh <SAMPLENAME> <READ1> <READ2> <GENOME>
## Details:
## SAMPLENAME = unique id/name that represents the sample. Eg. Sham1
## READ1      = FASTQ File pair 1 (can be in both compressed or uncompressed format)
## READ2      = FASTQ File pair 2 (can be in both compressed or uncompressed format)
## GENOME     = Genome build (hg19 / hg38 / mm9 / mm10)
##
## Note:
## If there are multiple read1/2 file, the files should be concanated first using cat in a single file
##
## Example:
## sh rnaseq_align.sh Sham1 read1.fq.gz read2.fq.gz mm9
##
##

##INPUT ARGUMENT (DO NOT EDIT)
##BASH DEFAULT USER ARGUMENT READING FROM COMMAND LINE
SAMPLENAME=$1
FASTQ_READ1=$2
FASTQ_READ2=$3
GENOME=$4

#AFTER TRIMMING (DO NOT EDIT)
TF1=$SAMPLENAME"_1_val_1.fq.gz"
TF2=$SAMPLENAME"_2_val_2.fq.gz"

##EDIT DURING SETUP
###TOOLS:
TRIM=trim_galore
STAR=/mnt/projects/rpd/apps/star-2.5.3a/bin/STAR
HTSEQCOUNT=/mnt/software/bin/htseq-count
SAMTOOLS=/mnt/bin/software

###DATABASE/ANNOTATION:
ADAPTOR=/mnt/projects/wlwtan/cardiac_epigenetics/foolab/jenny/mar2020/rnaseq/analysis_adaptor/illumina.fa
STARIND="/mnt/projects/rpd/genomes/"$GENOME"/star"
STARFASTA="/mnt/projects/rpd/genomes/"$GENOME"/"$GENOME".fa"
GTF="/mnt/projects/rpd/genomes/"$GENOME"/gtf/"$GENOME"_annotation.gtf"

###PARAMETERS:
THREAD=4


##ACTUAL COMMANDS:

mkdir $SAMPLENAME
cd $SAMPLENAME

## Trimming of adaptors and base quality
$TRIM --fastqc --gzip --length 100 --paired $FASTQ_READ1 $FASTQ_READ2


## STAR alignment
$STAR --runThreadN $THREAD --outFileNamePrefix rnaseqtrimmed --outSAMtype BAM Unsorted --genomeDir $STARIND --readFilesCommand zcat --readFilesIn $TF1 $TF2

## sort the bam files by name and count by htseq-count for EdgeR/DESeq analysis
$SAMTOOLS sort -n rnaseqtrimmedAligned.out.bam name_rnaseqtrimmedAligned.out
$HTSEQCOUNT -f bam -r name -s no -m union name_rnaseqtrimmedAligned.out.bam $GTF > count.txt



