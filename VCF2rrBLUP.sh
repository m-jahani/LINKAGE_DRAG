#!bin/bash
#the conversion process was looped over each sample because the software consumes tons of RAM
VCF=$1
bcftools query -l $VCF >${VCF%%vcf}sample_list

while read SAMPLE; do
    echo "$SAMPLE"
    vcftools --vcf $VCF --indv $SAMPLE --recode --recode-INFO-all --out $SAMPLE
    java -jar /DATA/home/mjahani/bin/NGSEPcore_4.2.0/NGSEPcore.jar VCFConverter -rrBLUP -i ${SAMPLE}.recode.vcf -o ${VCF%%vcf}_${SAMPLE}
    cat ${VCF%%vcf}_${SAMPLE}_rrBLUP.in >>${VCF%%.vcf}_rrBLUP.in
    rm ${VCF%%vcf}_${SAMPLE}_rrBLUP.in ${VCF%%vcf}_${SAMPLE}_rrBLUP_samples.txt ${SAMPLE}.log ${SAMPLE}.recode.vcf
done <${VCF%%vcf}sample_list
