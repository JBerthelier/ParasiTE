#!/bin/bash

###Intragenic TEs

#Sort intronic TEs
mkdir ../ParasiTE_output/intragenic_TEs &&
cd ../ParasiTE_output/intragenic_TEs &&

awk -F "\t" '$10 ~ "=RLC" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RLC-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=RLG" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RLG-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=Caulimoviridae" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > Caulimoviridae-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=RYX" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RYX-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=RIX" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RIX-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=RLX" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RLX-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=RXX" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RXX-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=RSX" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RSX-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=DTT" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTT-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=DTM" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTM-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=DTA" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTA-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=DTC" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTC-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=DTH" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTH-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=DXX" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DXX-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=DTX" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTX-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=DHX" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DHX-intronic_TEs.bed &&
awk -F "\t" '$10 ~ "=TXX" {print}' ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > TXX-intronic_TEs.bed &&
awk -F "\t" '{print $3-$2 }' OFS="\t" ../Results/Annotations_TEs/Intragenic_intronic_TEs.bed > TOTAL-intronic_TEs.bed

#Sort exonic TEs
awk -F "\t" '$10 ~ "=RLC" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RLC-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=RLG" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RLG-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=Caulimoviridae" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > Caulimoviridae-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=RYX" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RYX-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=RIX" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RIX-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=RLX" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RLX-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=RXX" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RXX-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=RSX" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RSX-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=DTT" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTT-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=DTM" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTM-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=DTA" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTA-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=DTC" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTC-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=DTH" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTH-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=DXX" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DXX-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=DTX" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTX-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=DHX" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DHX-exonic_TEs.bed &&
awk -F "\t" '$10 ~ "=TXX" {print}' ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > TXX-exonic_TEs.bed &&
awk -F "\t" '{print $3-$2 }' OFS="\t" ../Results/Annotations_TEs/Intragenic_exonic_TEs.bed > TOTAL-exonic_TEs.bed &&

#Sort ambigous TEs
awk -F "\t" '$10 ~ "=RLC" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RLC-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=RLG" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RLG-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=Caulimoviridae" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > Caulimoviridae-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=RYX" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RYX-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=RIX" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RIX-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=RLX" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RLX-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=RXX" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t"  > RXX-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=RSX" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > RSX-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=DTT" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTT-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=DTM" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTM-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=DTA" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTA-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=DTC" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTC-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=DTH" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTH-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=DXX" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DXX-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=DTX" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DTX-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=DHX" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > DHX-ambigous_TEs.bed &&
awk -F "\t" '$10 ~ "=TXX" {print}' ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed | awk -F "\t" '{print $3-$2 }' OFS="\t" > TXX-ambigous_TEs.bed &&
awk -F "\t" '{print $3-$2 }' OFS="\t" ../Results/Annotations_TEs/Intragenic_ambigous_TEs.bed > TOTAL-ambigous_TEs.bed &&

find -empty -type f -delete



