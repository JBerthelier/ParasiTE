#!/bin/bash

###Overlap TEs

mkdir ../ParasiTE_output/overlap_TEs-length &&
cd ../ParasiTE_output/overlap_TEs-length &&

#Sort overlap downstream TEs
awk -F "\t" '{$12 = ($3-$2)+1} 1' OFS='\t' ../Results/Annotations_TEs_events/Overlap_TEs_downstream_events.bed > Overlap-TE-downstream.figure.bed && 
awk -F "\t" '$10 ~ "=RLC" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > RLC-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=RLG" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > RLG-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=Caulimoviridae" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > Caulimoviridae-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=RIX" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > RIX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=RLX" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > RLX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=RXX" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > RXX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=RSX" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > RSX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DTT" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > DTT-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DTM" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > DTM-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DTA" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > DTA-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DTC" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > DTC-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DTH" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > DTH-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DXX" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > DXX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DTX" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > DTX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DHX" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > DHX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=TXX" {print}' ./Overlap-TE-downstream.figure.bed | awk -F "\t" '{print $12 }' > TXX-overlap-TE-downstream.bed &&
awk -F "\t" '{print $12 }' ./Overlap-TE-downstream.figure.bed > TOTAL-overlap-TE-downstream.bed

##Sort overlap upstream TEs

awk -F "\t" '{$12 = ($3-$2)+1} 1' OFS='\t'  ../Results/Annotations_TEs_events/Overlap_TEs_upstream_events.bed > Overlap-TE-upstream.figure.bed && 
awk -F "\t" '$10 ~ "=RLC" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > RLC-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=RLG" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > RLG-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=Caulimoviridae" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > Caulimoviridae-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=RIX" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > RIX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=RLX" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > RLX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=RXX" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > RXX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=RSX" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > RSX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DTT" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > DTT-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DTM" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > DTM-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DTA" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > DTA-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DTC" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > DTC-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DTH" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > DTH-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DXX" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > DXX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DTX" {print}'  ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > DTX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DHX" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > DHX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=TXX" {print}' ./Overlap-TE-upstream.figure.bed | awk -F "\t" '{print $12 }' > TXX-overlap-TE-upstream.bed &&
awk -F "\t" '{print $12 }' ./Overlap-TE-upstream.figure.bed > TOTAL-overlap-TE-upstream.bed &&

#Sort full TEs
awk -F "\t" '{$12 = ($3-$2)+1} 1' OFS='\t' ../Results/Annotations_TEs_events/Overlap_TEs_full_events.bed > Overlap_TEs_full.figure.bed&& 

awk -F "\t" '$10 ~ "=RLC" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > RLC-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=RLG" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > RLG-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=Caulimoviridae" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > Caulimoviridae-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=RIX" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > RIX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=RLX" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > RLX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=RXX" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > RXX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=RSX" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > RSX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DTT" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > DTT-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DTM" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > DTM-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DTA" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > DTA-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DTC" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > DTC-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DTH" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > DTH-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DXX" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > DXX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DTX" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > DTX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DHX" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > DHX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=TXX" {print}' ./Overlap_TEs_full.figure.bed| awk -F "\t" '{print $12 }' > TXX-overlap-TE-full.bed &&
awk -F "\t" '{print $12 }' ./Overlap_TEs_full.figure.bed> TOTAL-overlap-TE-full.bed

find -empty -type f -delete
