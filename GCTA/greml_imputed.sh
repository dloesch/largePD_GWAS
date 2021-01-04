#!/usr/bin/bash
#$ -q threaded.q
#$ -pe thread 4
#$ -P toconnor-lab
#$ -l mem_free=10G
#$ -e GCTA.log
#$ -o GCTA.log
#$ -cwd
#$ -N GCTA
#$ -hold_jid LD.GCTA


#set up files for merging
echo -n > multi_GRMs.test_group1.txt
echo -n	> multi_GRMs.test_group2.txt
echo -n	> multi_GRMs.test_group3.txt
echo -n	> multi_GRMs.test_group4.txt

for chr in {1..22};do
	echo chr$chr.test_group1 >> multi_GRMs.test_group1.txt
	echo chr$chr.test_group2 >> multi_GRMs.test_group2.txt
	echo chr$chr.test_group3 >> multi_GRMs.test_group3.txt
	echo chr$chr.test_group4 >> multi_GRMs.test_group4.txt
done


#merge grms by bin across chromosomes. Model does not converge if use all 88 stratified GRMs (4 bins *22 chromosomes) 
gcta64  --mgrm multi_GRMs.test_group1.txt  --make-grm  --out group1 --thread-num 4
gcta64  --mgrm multi_GRMs.test_group2.txt  --make-grm  --out group2 --thread-num 4
gcta64  --mgrm multi_GRMs.test_group3.txt  --make-grm  --out group3 --thread-num 4
gcta64  --mgrm multi_GRMs.test_group4.txt  --make-grm  --out group4 --thread-num 4


#inputs for GCTA analysis
pheno=large-PD.pheno
covar=large-PD.covar
qcovar=large-PD.qcovar
rels=large.unrelated_toberemoved.txt

echo "group1" > multi_GRMs.txt
echo "group2" >> multi_GRMs.txt
echo "group3" >> multi_GRMs.txt
echo "group4" >> multi_GRMs.txt

#run GCTA
gcta64 --reml --mgrm multi_GRMs.txt --out large.imputed.site --pheno $pheno --covar $covar --qcovar $qcovar --prevalence 0.005 --remove $rels --thread-num 4
