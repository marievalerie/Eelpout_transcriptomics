#! /bin/bash
#
#$ -S /bin/bash
#$ -N concat
#$ -cwd 
#$ -j n
#$ -m a
#$ -M m.brasseur@leibniz-zfmk.de


fwd='_1_polyx_trimmed_val_1.fq'
rev='_2_polyx_trimmed_val_2.fq' 

f='_trimmed_cat_1.fq'
r='_trimmed_cat_2.fq'


#sample='RNA_VaM_21'
#cat $sample*$fwd > $sample$f
#cat $sample*$rev > $sample$r

while read sample; do
  cat $sample*$fwd > $sample$f
  cat $sample*$rev > $sample$r
done < samples
