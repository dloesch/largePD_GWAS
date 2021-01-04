#!/bin/bash

#path to plink
plink=/usr/local/packages/plink-1.90.beta-3.6/bin/plink

#raw illumina data was parsed to lgen format using awk command
original=../original_data/large 

#convert to binary plink
M=5000
$plink --lfile $original --1 --snps-only 'just-acgt' --biallelic-only 'strict' --missing-genotype "-" \
    --output-missing-genotype 0 --make-bed --out large --memory $M

#basic qc
geno=large

#sex-check
$plink --bfile $geno --chr X, Y, XY --make-bed --out large.XY --memory $M
$plink --bfile large.XY --indep-pairwise 50 5 0.2 --out XY --memory $M
$plink --bfile large.XY --exclude XY.prune.out --check-sex --out XY --memory $M
#plot the output to see if there are any outliers

#sample missingness and remove X chromosome ()
$plink --bfile $geno --autosome --geno 0.2 --mind 0.1 --make-bed --out large.temp1 --memory $M

$plink --bfile large.temp1 --geno 0.05 --make-bed --out large.autosome --memory $M 

#hw equilibrium
#exclude variant if pvalue < 1E-6 in controls or < 1E-10 in cases
$plink --bfile large.autosome --hardy --out qc_hardy --memory $M

awk '($3 == "AFF" && $9 < 1E-10) || ($3 == "UNAFF" && $9 <1E-06) {print $2}' qc_hardy.hwe > HWE.failed_snps.txt

#identify duplicate sites

$plink --bfile large.autosome --list-duplicate-vars 'ids-only' 'suppress-first' --out qc_dupes --memory $M 

#identify ambiguous sites (CG/AT)
cat large.autosome.bim |awk -v OFS='\t' '($5 == "C" && $6 == "G") || ($5 == "G" && $6 == "C") || ($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") {print $2}' > CG_AT_sites.txt


#final filesets
cat qc_dupes.dupvar HWE.failed_snps.txt | sort -u > drop_snps.txt

$plink --bfile large.autosome --mac 1 --exclude drop_snps.txt --make-bed --out large.autosome.final --memory $M 

cat drop_snps.txt CG_AT_sites.txt | sort -u > drop_snps.CGAT.txt

$plink --bfile large.autosome --mac 1 --exclude drop_snps.CGAT.txt --make-bed --out large.autosome.final.noCGAT --memory $M 

#pruned filesets

$plink --bfile large.autosome.final --indep-pairwise 50 5 0.2 --out ld1 --memory $M
$plink --bfile large.autosome.final.noCGAT --indep-pairwise 50 5 0.2 --out ld2 --memory $M

$plink --bfile large.autosome.final --exclude ld1.prune.out --make-bed --out large.autosome.final.pruned --memory $M
$plink --bfile large.autosome.final.noCGAT --exclude ld2.prune.out --make-bed --out large.autusome.final.noCGAT.pruned --memory $M

#run genome to id any duplicates

$plink --bfile large.autosome.final.pruned --genome --out qc_dupe_subjects --memory $M 
awk -v OFS='\t' '($9 > 0.99){print $2,$4}' qc_dupe_subjects.genome > dupe_subject_pairs.txt
awk 'NR>1{print $1}' dupe_subject_pairs.txt >> dupe_subjects.txt && awk 'NR>1{print $2}' dupe_subject_pairs.txt >> dupe_subjects.txt

#use king to generate set of unrelated individuals
/usr/local/bin/king -b large.autosome.final.bed --unrelated --degree 2 --prefix large

mv largeunrelated_toberemoved.txt large.unrelated_toberemoved.txt

