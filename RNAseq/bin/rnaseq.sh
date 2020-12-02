####Pipeline created for single-cell RNA-seq analysis


F1=$1
F2=$2
SAMPLE=$3


##STAR alignment: align fastq files to genome
###Parameters:
STAR=/mnt/projects/rpd/apps/star-2.5.3a/bin/STAR
RSEM=/mnt/projects/rpd/apps/rsem-1.3.0/bin/rsem-calculate-expression
THREAD=4
STARIND=/mnt/projects/rpd/genomes/hg38/star
STARFASTA=/mnt/projects/rpd/genomes/hg38/hg38.fa
GTF=/mnt/projects/rpd/genomes/hg38/gtf/hg38_annotation.gtf
RSEMGENOME=/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/RSEM/index/hg38/hg38

TF1=$(echo $(basename $F1) | sed 's/.fastq.gz//')"_val_1.fq.gz"
TF2=$(echo $(basename $F2) | sed 's/.fastq.gz//')"_val_2.fq.gz"


mkdir $SAMPLE
cd $SAMPLE


trim_galore --fastqc --gzip --length 50 --paired $F1 $F2

$STAR --runThreadN $THREAD --genomeDir $STARIND --readFilesCommand zcat --outFileNamePrefix RNASEQ --outSAMtype BAM Unsorted --readFilesIn $TF1 $TF2
$STAR --runThreadN $THREAD --genomeDir $STARIND --sjdbFileChrStartEnd RNASEQSJ.out.tab --readFilesCommand zcat --outFileNamePrefix RNASEQ.2Pass --outSAMtype BAM Unsorted --readFilesIn $TF1 $TF2


## RUN RSEM
$RSEM --bowtie2 --bowtie2-path /mnt/software/bin/ --num-threads $THREAD --estimate-rspd --append-names --output-genome-bam --paired-end $TF1 $TF2 $RSEMGENOME RSEM

## Run Kallisto

## HTSEQCOUNT FROM STAR alignment
samtools sort -n RNASEQ.2PassAligned.out.bam name_RNASEQ.2PassAligned.out
/mnt/software/bin/htseq-count -f bam -r name -s no -m union name_RNASEQ.2PassAligned.out.bam $GTF > count.txt



