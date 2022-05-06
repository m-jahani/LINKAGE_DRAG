#!/bin/bash

awk '{ if ($1 == 1) { print $5"\t"$2"\t"$3} }' pcadmix.xxx2.regions.txt | bedtools sort >introgression_annuus.bed                                      #extract bed for annuus introgressions
bedtools genomecov -i introgression_annuus.bed -g HA412v2_chromosome.bed -bga | awk '{print $1"\t"$2"\t"$3"\t"$4/574 }' >introgression_annuus.coverage #calculate the base frequency
rm introgression_annuus.bed
sed 's/Ha412HOChr0/chr/g' introgression_annuus.coverage | sed 's/Ha412HOChr/chr/g' | sed '1i chr\tstart\tend\tvalue' >introgression_annuus.coverage.shiny
awk '{ if ($4 > 0) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_annuus.coverage.shiny | sed 's/\t/,/g' >introgression_annuus.coverage.shiny.csv                                                                                      #prepare voverage data for shiny
awk '{ if ($4 >= 0.02) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_annuus.coverage.shiny | bedtools merge | awk '{print $1","$2","$3","$1":"$2"-"$3}' | sed '1i chr,start,end,group' >introgression_annuus.merge.frq0.02.shiny.csv #introgression ranges for shiny rect
awk '{ if ($4 >= 0.02) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_annuus.coverage | bedtools sort | bedtools merge | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' >introgression_annuus.frq0.02.bed
rm introgression_annuus.coverage.shiny introgression_annuus.coverage

awk '{ if ($1 == 2) { print $5"\t"$2"\t"$3} }' pcadmix.xxx2.regions.txt | bedtools sort >introgression_other.bed
bedtools genomecov -i introgression_other.bed -g HA412v2_chromosome.bed -bga | awk '{print $1"\t"$2"\t"$3"\t"$4/574 }' >introgression_other.coverage
rm introgression_other.bed
sed 's/Ha412HOChr0/chr/g' introgression_other.coverage | sed 's/Ha412HOChr/chr/g' | sed '1i chr\tstart\tend\tvalue' >introgression_other.coverage.shiny
awk '{ if ($4 > 0) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_other.coverage.shiny | sed 's/\t/,/g' >introgression_other.coverage.shiny.csv
awk '{ if ($4 >= 0.02) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_other.coverage.shiny | bedtools merge | awk '{print $1","$2","$3","$1":"$2"-"$3}' | sed '1i chr,start,end,group' >introgression_other.merge.frq0.02.shiny.csv
awk '{ if ($4 >= 0.02) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_other.coverage | bedtools sort | bedtools merge | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' >introgression_other.frq0.02.bed
rm introgression_other.coverage.shiny introgression_other.coverage



for i in *.bed 
do
SAVE=SNP_${i%%.bed}
Rscript ASSIGNSNP.R SNP_LIST $i $SAVE
awk '{print $3}' $SAVE  > ${SAVE}_int
rm $SAVE
echo ${SAVE##SNP_introgression_} | cat - ${SAVE}_int > ${SAVE}_int_header
rm ${SAVE}_int
done

paste SNPS_ORDERD SNP_introgression_annuus.frq0.01_int_header SNP_introgression_annuus.frq0.02_int_header SNP_introgression_annuus.frq0.03_int_header SNP_introgression_annuus.frq0.04_int_header SNP_introgression_annuus.frq0.05_int_header SNP_introgression_other.frq0.01_int_header SNP_introgression_other.frq0.02_int_header SNP_introgression_other.frq0.03_int_header SNP_introgression_other.frq0.04_int_header SNP_introgression_other.frq0.05_int_header> all_ingt_freq_snp

rm *_int_header


bcftools query -l Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.maf0.3.286.sample.vcf.recode.vcf > sample_list

while read p; do
  echo "$p"
  vcftools --vcf Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.maf0.3.286.sample.vcf.recode.vcf --indv $p --recode --recode-INFO-all --out $p
  java -jar /DATA/home/mjahani/bin/NGSEPcore_4.2.0/NGSEPcore.jar VCFConverter -rrBLUP  -i ${p}.recode.vcf -o ${p}.recode.vcf
  rm ${p}.log ${p}.recode.vcf_rrBLUP_samples.txt ${p}.recode.vcf
  cat ${p}.recode.vcf_rrBLUP.in >> Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.maf0.3.286.sample_rrBLUP.in
done <sample_list

sed  '/##/d' Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.maf0.3.286.sample.vcf.recode.vcf | cut -f1,2 | awk '{print $1":"$2}' | sed '1d' > Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.maf0.3.286.sample_rrBLUP.SNP_LIST
mv sample_list Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.maf0.3.286.sample_rrBLUP.SAMPLE_LIST

PATH=/DATA/home/mjahani/LINKADE_DRAG/GP
MARKER=${PATH}/Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.maf0.3.286.sample_rrBLUP.in
PHENOTYPE=${PATH}/phenotype_common_georgia_corrected.csv
SAMPLES=${PATH}/Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.maf0.3.286.sample_rrBLUP.SAMPLE_LIST
SNPs=${PATH}/Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.maf0.3.286.sample_rrBLUP.SNP_LIST
SAVE_DIR=${PATH}/result
Rscript GnomicPrediction.R $MARKER $PHENOTYPE $SAMPLES $SNPs $SAVE_DIR

Rscript /DATA/home/mjahani/LINKADE_DRAG/GP/result/GenomicPrediction.R /DATA/home/mjahani/LINKADE_DRAG/GP/Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.maf0.3.286.sample_rrBLUP.in /DATA/home/mjahani/LINKADE_DRAG/GP/phenotype_common_georgia_corrected.csv /DATA/home/mjahani/LINKADE_DRAG/GP/Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.maf0.3.286.sample_rrBLUP.SAMPLE_LIST /DATA/home/mjahani/LINKADE_DRAG/GP/Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.maf0.3.286.sample_rrBLUP.SNP_LIST /DATA/home/mjahani/LINKADE_DRAG/GP/result

Rscript JOIN_FRQ2RESULT.R /DATA/home/mjahani/LINKADE_DRAG/introg_effect/all_ingt_freq_snp /DATA/home/mjahani/LINKADE_DRAG/introg_effect/BLUPS /DATA/home/mjahani/LINKADE_DRAG/introg_effect/BLUPS_intg_frq
Rscript JOIN_FRQ2RESULT.R /DATA/home/mjahani/LINKADE_DRAG/introg_effect/all_ingt_freq_snp /DATA/home/mjahani/LINKADE_DRAG/introg_effect/GWAS_RESULT_Beta /DATA/home/mjahani/LINKADE_DRAG/introg_effect/GWAS_Beta_intg_frq
Rscript JOIN_FRQ2RESULT.R /DATA/home/mjahani/LINKADE_DRAG/introg_effect/all_ingt_freq_snp /DATA/home/mjahani/LINKADE_DRAG/introg_effect/GWAS_RESULT_P /DATA/home/mjahani/LINKADE_DRAG/introg_effect/GWAS_P_intg_frq



bedtools shuffle -i /DATA/home/mjahani/LINKADE_DRAG/introg_effect/annuus.frq0.01.bed -g /DATA/home/mjahani/LINKADE_DRAG/introg_effect/HA412v2_chromosome.bed -excl /DATA/home/mjahani/LINKADE_DRAG/introg_effect/other.frq0.01.bed -noOverlapping > /DATA/home/mjahani/LINKADE_DRAG/introg_effect/annuus.frq0.01_other.frq0.01.bed

bedtools shuffle -i /DATA/home/mjahani/LINKADE_DRAG/introg_effect/annuus.frq0.01.bed -g /DATA/home/mjahani/LINKADE_DRAG/introg_effect/HA412v2_chromosome.bed -noOverlapping > /DATA/home/mjahani/LINKADE_DRAG/introg_effect/test.bed


awk '{ if ($1 == 1) { print $5"\t"$2"\t"$3} }' pcadmix.xxx2.regions.txt | bedtools sort >introgression_annuus.bed                                      #extract bed for annuus introgressions
bedtools genomecov -i introgression_annuus.bed -g HA412v2_chromosome.bed -bga | awk '{print $1"\t"$2"\t"$3"\t"$4/574 }' >introgression_annuus.coverage #calculate the base frequency
rm introgression_annuus.bed
sed 's/Ha412HOChr0/chr/g' introgression_annuus.coverage | sed 's/Ha412HOChr/chr/g' | sed '1i chr\tstart\tend\tvalue' >introgression_annuus.coverage.shiny
awk '{ if ($4 > 0) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_annuus.coverage.shiny | sed 's/\t/,/g' >introgression_annuus.coverage.shiny.csv                                                                                      #prepare voverage data for shiny
awk '{ if ($4 >= 0.05) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_annuus.coverage.shiny | bedtools merge | awk '{print $1","$2","$3","$1":"$2"-"$3}' | sed '1i chr,start,end,group' >introgression_annuus.coverage.merge.shiny.csv #introgression ranges for shiny rect
awk '{ if ($4 >= 0.05) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_annuus.coverage | bedtools sort | bedtools merge | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' >introgression_annuus.frq0.05.bed
rm introgression_annuus.coverage.shiny introgression_annuus.coverage

awk '{ if ($1 == 2) { print $5"\t"$2"\t"$3} }' pcadmix.xxx2.regions.txt | bedtools sort >introgression_other.bed
bedtools genomecov -i introgression_other.bed -g HA412v2_chromosome.bed -bga | awk '{print $1"\t"$2"\t"$3"\t"$4/574 }' >introgression_other.coverage
rm introgression_other.bed
sed 's/Ha412HOChr0/chr/g' introgression_other.coverage | sed 's/Ha412HOChr/chr/g' | sed '1i chr\tstart\tend\tvalue' >introgression_other.coverage.shiny
awk '{ if ($4 > 0) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_other.coverage.shiny | sed 's/\t/,/g' >introgression_other.coverage.shiny.csv
awk '{ if ($4 >= 0.05) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_other.coverage.shiny | bedtools merge | awk '{print $1","$2","$3","$1":"$2"-"$3}' | sed '1i chr,start,end,group' >introgression_other.coverage.merge.shiny.csv
awk '{ if ($4 >= 0.05) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_other.coverage | bedtools sort | bedtools merge | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' >introgression_other.frq0.05.bed
rm introgression_other.coverage.shiny introgression_other.coverage
