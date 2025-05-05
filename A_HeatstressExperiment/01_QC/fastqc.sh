#! /bin/bash
#
# -S /bin/bash
#$ -N fastqc
#$ -cwd
#$ -j n
#$ -m bea
#$ -M m.brasseur@leibniz-zfmk.de

module load fastqc/0.11.9

fastqc --noextract -o QC/ *.fq.gz 
