#!/bin/bash

for GWA in *ps; do
    awk '{print FILENAME , $1,"\t",$2,"\t",$4}' $GWA | sed 's/^Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99_//g' | sed 's/.ps//g' | sed 's/__/\t/g' >>all_GWAS
done
