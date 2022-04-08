#!/bin/bash

PATH=/DATA/home/mjahani/LINKADE_DRAG/new_method/GWAS/genotype
KINSHIP=${PATH}/Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.aIBS.kinf
COVARIATE=${PATH}/Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.cov
MAP=${PATH}/SAM_introgression_donor_ANNUUS.map
OUT_DIR=/DATA/home/mjahani/LINKADE_DRAG/new_method/GWAS/result

for PHENO in /DATA/home/mjahani/LINKADE_DRAG/new_method/GWAS/phenotype/*; do
    /DATA/home/mjahani/bin/emmax-intel64 -v -d 10 -t ${PATH}/SAM_introgression_donor_ANNUUS_maf0.03 -p $PHENO -k $KINSHIP -c $COVARIATE -o ${OUT_DIR}/SAM_introgression_donor_ANNUUS_maf0.03_${PHENO##*/}
    /DATA/home/mjahani/bin/emmax-intel64 -v -d 10 -t ${PATH}/SAM_introgression_donor_2nd_GERMPLASM_maf0.03 -p $PHENO -k $KINSHIP -c $COVARIATE -o ${OUT_DIR}/SAM_introgression_donor_2nd_GERMPLASM_maf0.03_${PHENO##*/}
    Rscript ${PATH}/manhattan_beta.R ${OUT_DIR}/SAM_introgression_donor_ANNUUS_maf0.03_${PHENO##*/}.ps ${OUT_DIR}/SAM_introgression_donor_2nd_GERMPLASM_maf0.03${PHENOTYPE##*/}.ps $MAP $OUT_DIR
done
