# LINKAGE_DRAG

[HA412v2_chromosome.txt](https://github.com/m-jahani/LINKAGE_DRAG/blob/main/HA412v2_chromosome.txt)

> Each chromosome length for Ha412v2.0 refrence genome

### pcadmix.xxx.regions.txt

> pcadmix result

### SAM_LIST

> list of samples in SAM population

### PCADMIX2TPED_EMMAX.R

> finds the overlaps between PCadmix result and 1kb genomewide windows and save the introgression variants in plink transposed format tped specificly for EMMAX, introgressed allele saved as allele 2

> runs as: PCADMIX2TPED_EMMAX.R HA412v2_chromosome.txt pcadmix.xxx.regions.txt SAM_LIST /save/directory

### GWAS_EMMAX_INROG.sh

> Runs EMMAX GWAS on introgression variants

### manhattan_beta.R

> Draws Manhattan and barcode plot for introgression GWAS, in the manhattan plots the position of introgressions in Y axis is base on P value and beta sign

### GWAS_SUMMARY.R

> Extracting ouliers from EMMAX GWAS result

### PCADMIX2TPED_rrBLUP.R

> finds the overlaps between PCadmix result and 1kb genomewide windows and save the introgression variants in plink transposed format tped specificly for rrBLUP, introgressed allele saved as allele 1

> runs as: PCADMIX2TPED_rrBLUP.R HA412v2_chromosome.txt pcadmix.xxx.regions.txt SAM_LIST /save/directory

### sample_list_for_phenotype

> list of samples with phenotype

### TPED2RRBLUP.sh

> convert tped to vcf and rrblup fomat and filters for maf

> runs as: PCADMIX2TPED_rrBLUP.R tpedprefix maf sample_list_for_phenotype

### GP_Permutation.R

> Rscript to calculate average effect size of all introgressions on each trait. To check the hypothesis of:
> difference between observed effect size and random event : 0
> It shuffles introgressions among samples for 10000 times and built a null distribution to calculate a P-value.
> P<0.05 was considered a significant difference

> Runs as: GP_Permutation.R MARKER PHENOTYPE SAVE_DIR DONOR

### CC_GP_Permutation.R

> same function as GP_Permutation.R, but it can take arrays on slurm Computecanada

> Runs as: GP_Permutation.R MARKER PHENOTYPE SAVE_DIR SLURM_ARRAY_TASK_ID DONOR

### slurm_GP_Permutation.sh

> script to submit CC_GP_Permutation.R on slurm system

### Heterosis_introgression_test.R

> compares average effect size and prioportion of positive effect between low hetrozygousity and full data set for introgression variants

> runs as: Heterosis_introgression_test.R MARKER MARKER_HOMOZ PHENOTYPE PHENOTYPE_HOMOZ DONOR SAVE_DIR

### VCF2FRQ.sh

> calculate frequency from VCF

> runs as: VCF2FRQ.sh VCF

### TEST_FRQUENCY_INTROGRESSION_EFFECT.R

> calculates if introgression frequency has significant impact on introgression effect

> runs as: TEST_FRQUENCY_INTROGRESSION_EFFECT.R MARKER PHENOTYPE IDs FRQs DONOR SAVE_DIR
