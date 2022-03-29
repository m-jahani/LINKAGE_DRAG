#!/bin/bash

awk '{ if ($1 == 1) { print $5"\t"$2"\t"$3} }' pcadmix.xxx2.regions.txt | bedtools sort >introgression_annuus.bed                                      #extract bed for annuus introgressions
bedtools genomecov -i introgression_annuus.bed -g HA412v2_chromosome.bed -bga | awk '{print $1"\t"$2"\t"$3"\t"$4/574 }' >introgression_annuus.coverage #calculate the base frequency
rm introgression_annuus.bed
sed 's/Ha412HOChr0/chr/g' introgression_annuus.coverage | sed 's/Ha412HOChr/chr/g' | sed '1i chr\tstart\tend\tvalue' >introgression_annuus.coverage.shiny
awk '{ if ($4 > 0) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_annuus.coverage.shiny | sed 's/\t/,/g' >introgression_annuus.coverage.shiny.csv                                                                                      #prepare voverage data for shiny
awk '{ if ($4 >= 0.03) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_annuus.coverage.shiny | bedtools merge | awk '{print $1","$2","$3","$1":"$2"-"$3}' | sed '1i chr,start,end,group' >introgression_annuus.coverage.merge.shiny.csv #introgression ranges for shiny rect
awk '{ if ($4 >= 0.03) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_annuus.coverage | bedtools sort | bedtools merge | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' >introgression_annuus.frq0.03.bed
rm introgression_annuus.coverage.shiny introgression_annuus.coverage

awk '{ if ($1 == 2) { print $5"\t"$2"\t"$3} }' pcadmix.xxx2.regions.txt | bedtools sort >introgression_other.bed
bedtools genomecov -i introgression_other.bed -g HA412v2_chromosome.bed -bga | awk '{print $1"\t"$2"\t"$3"\t"$4/574 }' >introgression_other.coverage
rm introgression_other.bed
sed 's/Ha412HOChr0/chr/g' introgression_other.coverage | sed 's/Ha412HOChr/chr/g' | sed '1i chr\tstart\tend\tvalue' >introgression_other.coverage.shiny
awk '{ if ($4 > 0) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_other.coverage.shiny | sed 's/\t/,/g' >introgression_other.coverage.shiny.csv
awk '{ if ($4 >= 0.03) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_other.coverage.shiny | bedtools merge | awk '{print $1","$2","$3","$1":"$2"-"$3}' | sed '1i chr,start,end,group' >introgression_other.coverage.merge.shiny.csv
awk '{ if ($4 >= 0.03) {print $1"\t"$2"\t"$3"\t"$4} }' introgression_other.coverage | bedtools sort | bedtools merge | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' >introgression_other.frq0.03.bed
rm introgression_other.coverage.shiny introgression_other.coverage
