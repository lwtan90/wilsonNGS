#!/bin/bash

export PYTHONPATH=$HOME/lib/python3.7/site-packages/:$PYTHONPATH
export PATH=$HOME/bin/:$PATH

GENOMESIZE=/mnt/projects/rpd/genomes/mm10/mm10.fa.fai

SAMPLE=$1
INPUT=$2



INDEX=/mnt/projects/rpd/genomes/mm10/bwa_path/nucleotide/mm10.fa

cd $SAMPLE



run_dfilter.sh -d=sorted_alignment_only_mapped.bed -c=$INPUT -o=H3K27ac.bed -f=bed -ks=60 -bs=100 -lpval=6 -wig


