#!/usr/bin/bash
#$ -q all.q  
#$ -P toconnor-lab
#$ -l mem_free=25G
#$ -e PREP.log
#$ -o PREP.log
#$ -N prep
#$ -cwd

#path to plink
plink=/usr/local/packages/plink-1.90.beta-3.6/bin/plink

#large-PD genotype file
geno=large.autosome.final

#memory parameter
M=5000

#create temporary vcf file with plink
$plink --bfile $geno --recode vcf-iid bgz --out temp.relabeled --output-chr 26 --memory $M --threads 2

new=temp.relabled.vcf.gz

#picardtools inputs
chain=/home/dloesch/WORKSPACE/lift/hg19ToHg38.over.chain.gz
ref=/home/dloesch/WORKSPACE/lift/hg38.fa.gz
reject=rejected_vars.vcf
input=temp.relabled.vcf.gz
output=large-PD.hg38.vcf

#run picard to liftover vcf
java -Xmx20g -jar /usr/local/packages/picard-2.18.7/picard.jar LiftoverVcf \
    I=$input \
    O=$output \
    CHAIN=$chain \
    REJECT=$reject \
    R=$ref \
    RECOVER_SWAPPED_REF_ALT=true \
    MAX_RECORDS_IN_RAM=8000	

#bgzip and index output
bgzip -f $reject 
bgzip -f $output
tabix -f $output.gz

#clean-up temp files
rm temp.*
