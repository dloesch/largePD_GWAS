#!/usr/bin/bash

###Can use wrapper script to submit as batch jobs###
prefix=large

###Step1: convert to GDS format#######

for chr in {1..22};
do
	vcffile=../vcf/chr$chr.dose.vcf.gz
	gdsfile=../gds/$prefix.chr$chr.gds

	Rscript ../scripts/GENESIS/vcf2gds.R $vcffile $gdsfile
done

###step2: prune######### 
#using genotype snps only for GRM and PCA
#ld prning of genotyped file using PLINK
binary=../plink/large.autosome.final
pruned_vcf=../vcf/$prefix.pruned.vcf.gz
dupes=../plink/dupevars.after_flipping.txt #remove variants flagged as duplicates

plink --bfile $binary --exclude $dupes --make-bed --out temp
 
plink --bfile temp --indep-pairwise 50 5 0.2 --threads 2 --out ld

plink --bfile temp --recode vcf-iid bgz --exclude ld.prune.out --out ../vcf/$prefix.pruned


##convert pruned vcf file into gds file
pruned_vcf=../vcf/$prefix.pruned.vcf.gz
pruned_gds=../gds/$prefix.pruned.gds

Rscript ../scripts/GENESIS/vcf2gds.genotyped.R $pruned_vcf $pruned_gds

#make GCTA GRM
Rscript ../scripts/GENESIS/grm_allchrs.R $prefix $pruned_gds

####step4: kinship via KING robust algorithim
snpgds=$prefix.pruned.snp.gds
Rscript ../scripts/GENESIS/ibd.R $prefix $snpgds

########step4: PCA##########

#step4A:PC-AIR
king_file=$prefix.ibd.RData
snpgds=$prefix.pruned.snp.gds
Rscript ../scripts/GENESIS/pcair.R $snpgds $king_file $prefix king


###step4B: PC-Relate
gdsfile=$pruned_gds
Rscript ../scripts/GENESIS/pcrelate.R $prefix $gdsfile

####step4B: PC-Air, take 2 #################
king_file=$prefix.ibd.RData
snpgds=$prefix.pruned.snp.gds
Rscript ../scripts/GENESIS/pcair.R $snpgds $king_file $prefix pcrelate


##########Step5: Create phenotype file#############
#phenotype file. needs to be csv or tab separated
#first column should be subject ids
pheno=large-PD.pheno.05_2019.csv
covars="AGE,SEX" #separated by comma
n_pcs=5 #specify desired number of pcs
trait=PD_STATUS #trait

Rscript ../scripts/GENESIS/pheno.R $pheno $prefix $covars $n_pcs $trait

#####step6: Fit null model #########
grm=$prefix.ALL.grm.RData
grm_type=GCTA #can also use pcrelate 

Rscript ../scripts/GENESIS/null_model.R $prefix $grm $grm_type


#####Step 7: score test ###########
for chr in {1..22}; do 
	gdsfile=../gds/$prefix.chr$chr.gds
	Rscript ../scripts/GENESIS/assoc.R $chr $gdsfile $prefix 
done
