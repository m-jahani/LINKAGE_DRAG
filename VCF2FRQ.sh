#!/bin/bash
VCF=$1

vcftools --vcf $VCF --freq --out ${VCF%%.vcf}

awk -F "\t|:" '{print $1":"$2,$8}' ${VCF%%.vcf}.frq >${VCF%%.vcf}.frequency
rm ${VCF%%.vcf}.frq
