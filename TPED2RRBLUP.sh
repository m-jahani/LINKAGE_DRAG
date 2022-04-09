#!/bin/bash
TFILE=$1
plink --tfile $TFILE --recode vcf --out $TFILE
java -jar /DATA/home/mjahani/bin/NGSEPcore_4.2.0/NGSEPcore.jar VCFConverter -rrBLUP -i ${TFILE}.vcf -o $TFILE
