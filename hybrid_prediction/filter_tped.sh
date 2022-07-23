#!/bin/bash

awk '{print $1,"\t",$1}' /DATA/home/mjahani/LINKADE_DRAG/new_method/GP1/genotype/SAM_introgression_donor_ANNUUS_sample_list > SAM_introgression_donor_ANNUUS_sample_list_2column
plink --tfile /DATA/home/mjahani/LINKADE_DRAG/new_method/GP1/genotype/SAM_introgression_donor_ANNUUS --maf 0.03 --keep SAM_introgression_donor_ANNUUS_sample_list_2column --recode --transpose --out SAM_introgression_donor_ANNUUS.286_sample_MAF0.03

awk '{print $1,"\t",$1}' /DATA/home/mjahani/LINKADE_DRAG/new_method/GP1/genotype/SAM_introgression_donor_2nd_GERMPLASM_sample_list > SAM_introgression_donor_2nd_GERMPLASM_sample_list_2column
plink --tfile /DATA/home/mjahani/LINKADE_DRAG/new_method/GP1/genotype/SAM_introgression_donor_2nd_GERMPLASM --maf 0.03 --keep SAM_introgression_donor_ANNUUS_sample_list_2column --recode --transpose --out SAM_introgression_donor_2nd_GERMPLASM.286_sample_MAF0.03
