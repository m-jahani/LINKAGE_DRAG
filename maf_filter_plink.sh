#!/bin/bash
plink --tfile SAM_introgression_donor_ANNUUS --keep sample_list_for_phenotype --maf 0.03 --make-bed --out SAM_introgression_donor_ANNUUS_maf0.03        #505038
plink --tfile SAM_introgression_donor_2nd_GERMPLASM --keep sample_list_for_phenotype --maf 0.03 --make-bed --out SAM_introgression_donor_ANNUUS_maf0.03 #5243
