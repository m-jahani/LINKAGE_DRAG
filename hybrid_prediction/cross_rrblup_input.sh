#!/bin/bash

conda activate VCFTOOLS

TEPD=$1     #/DATA/home/mjahani/LINKADE_DRAG/new_method/GP1/hybrid_prediction/SAM_introgression_donor_ANNUUS.286_sample_MAF0.03.tped
TFAM=$2     #/DATA/home/mjahani/LINKADE_DRAG/new_method/GP1/hybrid_prediction/SAM_introgression_donor_ANNUUS.286_sample_MAF0.03.tfam
SAVE_DIR=$3 #/DATA/home/mjahani/LINKADE_DRAG/new_method/GP1/hybrid_prediction/crosses/ANNUUS

Rscript /DATA/home/mjahani/LINKADE_DRAG/new_method/GP1/hybrid_prediction/Cross_simulation.R $TEPD $TFAM ${SAVE_DIR}/plink

cd ${SAVE_DIR}/plink

for FILE in $(ls -l S* | awk '{print $9}' | sed 's/.tped//g' | sed 's/.tfam//g' | sed 's/.map//g' | sort -u); do
    plink --tfile $FILE --recode vcf-iid --out ${SAVE_DIR}/vcf/${FILE}
    java -jar /DATA/home/mjahani/bin/NGSEPcore_4.2.0/NGSEPcore.jar VCFConverter -rrBLUP -i ${SAVE_DIR}/vcf/${FILE}.vcf -o ${SAVE_DIR}/rrblup/${FILE}
done
