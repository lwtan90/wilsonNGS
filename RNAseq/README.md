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
```
star --runThreadN $THREAD --outFileNamePrefix rnaseqtrimmed --outSAMtype BAM Unsorted --genomeDir $STARIND --readFilesCommand zcat --readFilesIn read1.fastq.gz read2.fastq.gz  
```  

Input files:  
```
a. read1.fastq.gz (trimmed fastq file R1)  
b. read2.fastq.gz (trimmed fastq file R2)  
c. $STARIND = indices required for the mapping of sequencing reads.  

```

The details on various parameters can be found [here](https://github.com/alexdobin/STAR).
  
Running LOG:
```
Mar 13 23:10:09 ..... started STAR run  
Mar 13 23:10:10 ..... loading genome  
Mar 13 23:11:59 ..... started mapping  
Mar 14 02:11:40 ..... finished successfully  
```

  
#### 3. Generation of RAW gene count using htseq-count
```
/mnt/software/bin/htseq-count -f bam -r name -s no -m union name_rnaseqtrimmedAligned.out.bam $GTF > count.txt  
```   

Input files:
```
a. BAM file (sorted by read name / position )  
b. GTF files  
```

Output file (you need this for most analysis)
```
Filename: count.txt (gene count needed for EdgeR / DESeq)  
Format:  
<Gene ID> <unnormalized count>  
Example:
ENSMUSG00000093369.1	0  
ENSMUSG00000093370.1	3  
ENSMUSG00000093371.1	15  
ENSMUSG00000093373.1	0  
ENSMUSG00000093374.1	0  
__no_feature	9269724  
__ambiguous	97362  
__too_low_aQual	0  
__not_aligned	0  
__alignment_not_unique	3113139  

```
  
Description of the last 5 lines of each files:  
__no_feature: reads (or read pairs) which could not be assigned to any feature   
__ambiguous:  reads (or read pairs) which could have been assigned to more than one feature  
__too_low_aQual: reads (or read pairs) which were skipped due to the -a option  
__not_aligned: reads (or read pairs) in the SAM file without alignment  
__alignment_not_unique: reads (or read pairs) with more than one reported alignment.  
  
  
Running LOG:  
```
...  
31700000 SAM alignment record pairs processed.  
31800000 SAM alignment record pairs processed.  
31900000 SAM alignment record pairs processed.  
32000000 SAM alignment record pairs processed.  
32100000 SAM alignment record pairs processed.  
32200000 SAM alignment record pairs processed.  
32300000 SAM alignment record pairs processed.  
32400000 SAM alignment record pairs processed.  
32500000 SAM alignment record pairs processed.  
32600000 SAM alignment record pairs processed.  
32700000 SAM alignment record pairs processed.  
32800000 SAM alignment record pairs processed.  
32889607 SAM alignment pairs processed.  
```
