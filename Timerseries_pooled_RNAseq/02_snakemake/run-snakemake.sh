#! /bin/bash
#
#$ -S /bin/bash
#$ -N snakejob
#$ -cwd
#$ -j n
#$ -m e
#$ -M m.brasseur@leibniz-zfmk.de

module load anaconda3/2020.02
conda activate snakemake
module load samtools/1.10
module load opt-python

a=1; b=$NSLOTS
THREADS=$((b-a))

snakemake --cores ${THREADS} 
