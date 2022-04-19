#!/bin/bash
TFILE=$1
MAF=$2

plink --tfile $TFILE --maf $MAF --recode vcf-iid --out $TFILE

bcftools query -l ${TFILE}.vcf >${TFILE}_sample_list
grep -v "^#" ${TFILE}.vcf | awk '{print $1":"$2}' >${TFILE}_variant_list
while read SAMPLE; do
    echo "$SAMPLE"
    vcftools --vcf ${TFILE}.vcf --indv $SAMPLE --recode --recode-INFO-all --out $SAMPLE
    java -jar /DATA/home/mjahani/bin/NGSEPcore_4.2.0/NGSEPcore.jar VCFConverter -rrBLUP -i ${SAMPLE}.recode.vcf -o ${TFILE}_${SAMPLE}
    cat ${TFILE}_${SAMPLE}_rrBLUP.in >>${TFILE}_rrBLUP.in
    rm ${TFILE}_${SAMPLE}_rrBLUP.in ${TFILE}_${SAMPLE}_rrBLUP_samples.txt ${SAMPLE}.log ${SAMPLE}.recode.vcf
done <${TFILE}_variant_list
