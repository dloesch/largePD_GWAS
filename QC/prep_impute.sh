#!/usr/bin/bash
#$ -q all.q  
#$ -P toconnor-lab
#$ -l mem_free=10G
#$ -e PREP2.log
#$ -o PREP2.log
#$ -cwd
#$ -N prep2
#$ -hold_jid prep

#chromosome
chr=$1

#reference genotype. Used a subset of TOPMed to align large-PD with TOPMed
ref=$2

#path to plink binary
plink=/usr/local/packages/plink-1.90.beta-3.6/bin/plink


#input genotype file
geno=large-PD.hg38.vcf.gz
 
#path to CONFORM-GT program
conform=/usr/local/packages/conform-gt/conform-gt.jar


#run conform-gt to "align" large-PD with TOPMed for imputation
java -jar $conform ref=$ref gt=$geno chrom=chr$chr out=large-PD.chr$chr.HG38.CONFORM match=POS


