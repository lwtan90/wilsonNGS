## RNA-seq pipelines  
  
  
This pipeline will take in paired-end FASTQ files, and produce gene count as output.  



```  
F1=$1
F2=$2
SAMPLE=$3

mkdir $SAMPLE
cd $SAMPLE

##STAR alignment: align fastq files to genome
###Parameters:
STAR=/mnt/projects/rpd/apps/star-2.5.3a/bin/STAR
THREAD=1
STARIND=/mnt/projects/rpd/genomes/mm9/star
STARFASTA=/mnt/projects/rpd/genomes/mm9/mm9.fa
GTF=/mnt/projects/rpd/genomes/mm9/gtf/mm9_annotation.gtf

java -jar /mnt/projects/wlwtan/cardiac_epigenetics/pipeline/trimmomatics/Trimmomatic-0.39/trimmomatic-0.39.jar PE $F1 $F2 read1.fastq.gz read1.un.fastq.gz read2.fastq.gz read2.un.fastq.gz ILLUMINACLIP:/mnt/projects/wlwtan/cardiac_epigenetics/foolab/jenny/mar2020/rnaseq/analysis_adaptor/illumina.fa:2:30:10 LEADING:28 TRAILING:28 SLIDINGWINDOW:4:28 MINLEN:100

$STAR --runThreadN $THREAD --outFileNamePrefix rnaseqtrimmed --outSAMtype BAM Unsorted --genomeDir $STARIND --readFilesCommand zcat --readFilesIn read1.fastq.gz read2.fastq.gz
samtools sort -n rnaseqtrimmedAligned.out.bam name_rnaseqtrimmedAligned.out
/mnt/software/bin/htseq-count -f bam -r name -s no -m union name_rnaseqtrimmedAligned.out.bam $GTF > count.txt
```

