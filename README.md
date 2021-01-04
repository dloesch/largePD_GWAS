# largePD_GWAS
 Scripts for running GWAS and related analyses for the LARGE-PD GWAS manuscript

## Scripts from GWAS manuscript 
1. QC: contains bash scripts for QC and prepping LARGE-PD genotype data for imputation (process_LARGEPD.sh, lift_LARGEPD.sh, prep_impute.sh)
2. GWAS: contains a bash script outlining GENESIS workflow (GWAS_workflow.sh) and R scripts for running GENESIS modfied from the TOPMed pipeline https://github.com/UW-GAC/analysis_pipeline
4. POST_GWAS: contains R scripts for performing conditional analysis (conditional_GMMAT.R) and replicated known PD loci (replicate_GENESIS.R and replicate_PD_GWAS.R)
5. ADMIX_MAPPING: contains R scripts overlaying adnixture mapping and GWAS results (overlay_AM_GWAS.R) and for plotting admixture mapping results (plot_AM.R)
6. GCTA: contains scripts for running GCTA to esimate additive heritability using imputed data (ldscore_GRM.sh, ld_score_bin.R, greml_imputed.sh)
7. ADMIXTURE: contains script for running ADMIXTURE (run_admix.R) and for plotting admixture results (plot_admixture.R)

## Additional analyses:
1. Local ancestry was conducted using RFMIX and scripts developed by Martin et al. 2017 (https://github.com/armartin/ancestry_pipeline/)
2. Admixture mapping analysis was conducted using a logistic mixed model implemented with the admixMap function in the GENESIS R package (https://www.bioconductor.org/packages/release/bioc/html/GENESIS.html). 
3. Statistical fine mapping of GWAS results was done with PAINTOR 3.0 using a custom pipeline (https://github.com/dloesch/postGWAS/Finemap_pipeline)
4. GWAS plot were genereated using previously generated code (https://github.com/dloesch/postGWAS/Plot)

## Relevant Citations
### LARGE-PD
1. Loesch D, Horimoto ARVR, Heilbron K, et al. Characterizing the genetic architecture of Parkinson’s disease in Latinos. medRxiv. Published online November 12, 2020:2020.11.09.20227124. 
2. Zabetian CP, Mata IF, Latin American Research Consortium on the Genetics of PD (LARGE-PD). LARGE-PD: Examining the genetics of Parkinson’s disease in Latin America. Mov Disord. 2017;32(9):1330-1331. 
### Software
1. CONFORM-GT: https://faculty.washington.edu/browning/conform-gt.html
2. GENESIS: Conomos MP, Gogarten SM, Brown L, et al. GENESIS: GENetic EStimation and Inference in Structured Samples (GENESIS): Statistical Methods for Analyzing Genetic Data from Samples with Population Structure and/or Relatedness. Bioconductor version: Release (3.10); 2020. 
3. GCTA: Yang J, Bakshi A, Zhu Z, et al. Genetic variance estimation with imputed variants finds negligible missing heritability for human height and body mass index. Nature Genetics. 2015;47(10):1114-1120. doi:10.1038/ng.3390
4. PLINK 1.9: Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience. 2015;4.
5. PicardTools: http://broadinstitute.github.io/picard/
6. TOPMed Imputation panel: Taliun D, Harris DN, Kessler MD, et al. Sequencing of 53,831 Diverse Genomes from the NHLBI TOPMed Program. Genomics; 2019.
7. GMMAT: Chen H, Wang C, Conomos MP, et al. Control for Population Structure and Relatedness for Binary Traits in Genetic Association Studies via Logistic Mixed Models. Am J Hum Genet. 2016;98(4):653-666.
8. PAINTOR: Kichaev G, Yang W-Y, Lindstrom S, et al. Integrating Functional Data to Prioritize Causal Variants in Statistical Fine-Mapping Studies. PLOS Genetics. 2014;10(10):e1004722.  
9. RFMix: Maples BK, Gravel S, Kenny EE, Bustamante CD. RFMix: A Discriminative Modeling Approach for Rapid and Robust Local-Ancestry Inference. Am J Hum Genet. 2013;93(2):278-288. doi:10.1016/j.ajhg.2013.06.020
### PD GWAS
1. Nalls MA, Blauwendraat C, Vallerga CL, et al. Identification of novel risk loci, causal insights, and heritable risk for Parkinson’s disease: a meta-analysis of genome-wide association studies. The Lancet Neurology. 2019;18(12):1091-1102. 
2. Foo JN, Chew EGY, Chung SJ, et al. Identification of Risk Loci for Parkinson Disease in Asians and Comparison of Risk Between Asians and Europeans: A Genome-Wide Association Study. JAMA Neurol. 