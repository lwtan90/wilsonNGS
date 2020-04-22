## RNA-seq pipelines  
  
  
This pipeline will take in paired-end FASTQ files, and produce gene count as output.  
This is a skeleton for basic analysis.  
There are various added analysis described at the bottom of the page.  
  

```  

#!/bin/bash

#FASTQ FILE R1 $F1; R2 $F2;  
F1=$1
F2=$2
SAMPLE=$3

##STAR alignment: align fastq files to genome
###Parameters:
STAR=/mnt/projects/rpd/apps/star-2.5.3a/bin/STAR  
THREAD=1  
STARIND=/mnt/projects/rpd/genomes/mm9/star  
STARFASTA=/mnt/projects/rpd/genomes/mm9/mm9.fa  
GTF=/mnt/projects/rpd/genomes/mm9/gtf/mm9_annotation.gtf  

mkdir $SAMPLE
cd $SAMPLE

java -jar /mnt/projects/wlwtan/cardiac_epigenetics/pipeline/trimmomatics/Trimmomatic-0.39/trimmomatic-0.39.jar PE $F1 $F2 read1.fastq.gz read1.un.fastq.gz read2.fastq.gz read2.un.fastq.gz ILLUMINACLIP:illumina.fa:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:100

$STAR --runThreadN $THREAD --outFileNamePrefix rnaseqtrimmed --outSAMtype BAM Unsorted --genomeDir $STARIND --readFilesCommand zcat --readFilesIn read1.fastq.gz read2.fastq.gz
samtools sort -n rnaseqtrimmedAligned.out.bam name_rnaseqtrimmedAligned.out
/mnt/software/bin/htseq-count -f bam -r name -s no -m union name_rnaseqtrimmedAligned.out.bam $GTF > count.txt

```

  
  
## Detailed Description of the basic pipeline  
#### 1. Trimming of the reads before mapping  

```
java -jar trimmomatic-0.39.jar PE $F1 $F2 read1.fastq.gz read1.un.fastq.gz read2.fastq.gz read2.un.fastq.gz ILLUMINACLIP:illumina.fa:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:100
```  

Input files:  
```  
a. F1 = FASTQ File Read 1  
b. F2 = FASTQ File Read 2  
```  

The details on various parameters can be found [here](http://www.usadellab.org/cms/?page=trimmomatic).    

  
This command trim the reads:  
a. base quality  
b. adaptors / index  

Running LOG:  
```
Using Medium Clipping Sequence: 'AGATCGGAAGAGCACACGTCT'  
Using Medium Clipping Sequence: 'CTGTCTCTTATACACATCT'  
Using Long Clipping Sequence: 'CTGTCTCTTATACACATCT+AGATGTGTATAAGAGACAG'  
Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'  
Using Short Clipping Sequence: 'AGATCGGAAGAG'  
Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'  
Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAAC'  
Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGA'  
Using Medium Clipping Sequence: 'TGGAATTCTCGGGTGCCAAGG'  
ILLUMINACLIP: Using 0 prefix pairs, 9 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences  
Quality encoding detected as phred33  
Input Read Pairs: 57690695 Both Surviving: 40872736 (70.85%) Forward Only Surviving: 9407255 (16.31%) Reverse Only Surviving: 1621021 (2.81%) Dropped: 5789683 (10.04%)  
TrimmomaticPE: Completed successfully  

``` 
   
#### 2. Mapping of sequencing reads  



