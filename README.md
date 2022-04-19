# LINKAGE_DRAG

### HA412v2_chromosome.txt

> Each chromosome length for Ha412v2.0 refrence genome

### pcadmix.xxx.regions.txt

> pcadmix result

### SAM_LIST

> list of samples in SAM population

### PCADMIX2TPED_EMMAX.R

> finds the overlaps between PCadmix result and 1kb genomewide windows and save the introgression variants in plink transposed format tped specificly for EMMAX, introgressed allele saved as allele 2

> runs as: PCADMIX2TPED_EMMAX.R HA412v2_chromosome.txt pcadmix.xxx.regions.txt SAM_LIST /save/directory

### PCADMIX2TPED_rrBLUP.R

> finds the overlaps between PCadmix result and 1kb genomewide windows and save the introgression variants in plink transposed format tped specificly for rrBLUP, introgressed allele saved as allele 1

> runs as: PCADMIX2TPED_rrBLUP.R HA412v2_chromosome.txt pcadmix.xxx.regions.txt SAM_LIST /save/directory
