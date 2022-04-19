#!/bin/bash
TFILE=$1
MAF=$2
SAMPLE=$3

plink --tfile $TFILE --maf $MAF --keep $SAMPLE --recode vcf-iid --out ${TFILE}_maf${MAF}

bcftools query -l ${TFILE}_maf${MAF}.vcf >${TFILE}_sample_list
grep -v "^#" ${TFILE}.vcf | awk '{print $1":"$2}' >${TFILE}_variant_list
while read SAMPLE; do
    echo "$SAMPLE"
    vcftools --vcf ${TFILE}_maf${MAF}.vcf --indv $SAMPLE --recode --recode-INFO-all --out $SAMPLE
    java -jar /DATA/home/mjahani/bin/NGSEPcore_4.2.0/NGSEPcore.jar VCFConverter -rrBLUP -i ${SAMPLE}.recode.vcf -o ${TFILE}_${SAMPLE}
    cat ${TFILE}_${SAMPLE}_rrBLUP.in | sed 's/-1/2/g' | sed 's/1/-1/g' | sed 's/2/1/g' >>${TFILE}_maf${MAF}_rrBLUP.in
    rm ${TFILE}_${SAMPLE}_rrBLUP.in ${TFILE}_${SAMPLE}_rrBLUP_samples.txt ${SAMPLE}.log ${SAMPLE}.recode.vcf
done <${TFILE}_sample_list

rm {TFILE}_maf${MAF}.vcf
