#!/bin/bash

mkdir ../ParasiTE_output/neighbor_TEs &&
cd ../ParasiTE_output/neighbor_TEs &&

awk -F "\t" '$20 ~ "=RLC" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > RLC-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=RLG" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > RLG-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=Caulimoviridae" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > Caulimoviridae-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=RIX" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > RIX-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=RLX" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > RLX-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=RXX" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > RXX-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=RSX" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > RSX-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=DTT" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DTT-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=DTM" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DTM-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=DTA" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DTA-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=DTC" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DTC-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=DTH" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DTH-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=DXX" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DXX-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=DTX" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DTX-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=DHX" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DHX-neighbor_TEs-downstream.bed &&
awk -F "\t" '$20 ~ "=TXX" {print}' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > TXX-neighbor_TEs-downstream.bed &&
awk -F "\t" '{print $21 }' ../neighbor_TEs_wo-overlap-downstream.uniq.length.sorted.bed > TOTAL-neighbor_TEs-downstream.bed &&


#Sort upstream neighbor TEs
awk -F "\t" '$20 ~ "=RLC" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > RLC-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=RLG" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > RLG-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=Caulimoviridae" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > Caulimoviridae-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=RIX" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > RIX-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=RLX" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > RLX-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=RXX" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > RXX-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=RSX" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > RSX-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=DTT" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DTT-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=DTM" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DTM-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=DTA" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DTA-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=DTC" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DTC-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=DTH" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DTH-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=DXX" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DXX-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=DTX" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DTX-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=DHX" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > DHX-neighbor_TEs-upstream.bed &&
awk -F "\t" '$20 ~ "=TXX" {print}' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed | awk -F "\t" '{print $21 }' > TXX-neighbor_TEs-upstream.bed &&
awk -F "\t" '{print $21 }' ../neighbor_TEs_wo-overlap-upstream.uniq.length.sorted.bed > TOTAL-neighbor_TEs-upstream.bed


find -empty -type f -delete
