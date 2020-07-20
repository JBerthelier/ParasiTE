#!/bin/bash

###Intragenic TEs

#Sort overlap downstream TEs
cd /home/jeremy/galaxy/tools/Pipeline/ParasiTE/WORK &&
rm -Rf /home/jeremy/galaxy/tools/Pipeline/ParasiTE/WORK/overlap_TEs-overlap-length &&
mkdir /home/jeremy/galaxy/tools/Pipeline/ParasiTE/WORK/overlap_TEs-overlap-length &&
cd /home/jeremy/galaxy/tools/Pipeline/ParasiTE/WORK/overlap_TEs-overlap-length &&
awk -F "\t" '$10 ~ "=RLC" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > RLC-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=RLG" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > RLG-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=Caulimoviridae" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > Caulimoviridae-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=RIX" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > RIX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=RLX" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > RLX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=RXX" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > RXX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=RSX" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > RSX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DTM" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DTM-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DTA" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DTA-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DTC" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DTC-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DTH" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DTH-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DXX" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DXX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DTX" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DTX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=DHX" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DHX-overlap-TE-downstream.bed &&
awk -F "\t" '$10 ~ "=TXX" {print}' ../Overlap-TE-downstream.final.sorted.bed | awk -F "\t" '{print $21 }' > TXX-overlap-TE-downstream.bed &&
awk -F "\t" '{print $21 }' ../Overlap-TE-downstream.final.sorted.bed > TOTAL-overlap-TE-downstream.bed

##Sort overlap upstream TEs
awk -F "\t" '$10 ~ "=RLC" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > RLC-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=RLG" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > RLG-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=Caulimoviridae" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > Caulimoviridae-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=RIX" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > RIX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=RLX" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > RLX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=RXX" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > RXX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=RSX" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > RSX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DTM" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DTM-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DTA" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DTA-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DTC" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DTC-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DTH" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DTH-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DXX" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DXX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DTX" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DTX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=DHX" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > DHX-overlap-TE-upstream.bed &&
awk -F "\t" '$10 ~ "=TXX" {print}' ../Overlap-TE-upstream.final.sorted.bed | awk -F "\t" '{print $21 }' > TXX-overlap-TE-upstream.bed &&
awk -F "\t" '{print $21 }' ../Overlap-TE-upstream.final.sorted.bed > TOTAL-overlap-TE-upstream.bed &&

#Sort full TEs
awk -F "\t" '$10 ~ "=RLC" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > RLC-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=RLG" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > RLG-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=RIX" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > RIX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=RLX" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > RLX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=RXX" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > RXX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=RSX" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > RSX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DTM" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > DTM-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DTA" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > DTA-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DTC" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > DTC-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DTH" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > DTH-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DXX" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > DXX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DTX" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > DTX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=DHX" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > DHX-overlap-TE-full.bed &&
awk -F "\t" '$10 ~ "=TXX" {print}' ../TEs_fulloverlap_genes.details.uniq.sorted.bed | awk -F "\t" '{print $21 }' > TXX-overlap-TE-full.bed &&
awk -F "\t" '{print $21 }' ../TEs_fulloverlap_genes.details.uniq.sorted.bed > TOTAL-overlap-TE-full.bed
