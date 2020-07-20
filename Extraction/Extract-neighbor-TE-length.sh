#!/bin/bash

#Get downstream neighbor TEs

mkdir ../ParasiTE_output/neighbor_TEs-length &&
cd ../ParasiTE_output/neighbor_TEs-length &&

awk -F "\t" '{$12 = ($3-$2)} 1' OFS='\t' ../Results/Annotations_TEs/Neighbor_TEs_downstream.bed > Neighbor_TEs_downstream.figure.bed && 
awk -F "\t" '$10~ "=RLC" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > RLC-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=RLG" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > RLG-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=RYX" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > RYX-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=Caulimoviridae" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > Caulimoviridae-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=RIX" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > RIX-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=RLX" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > RLX-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=RXX" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > RXX-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=RSX" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > RSX-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=DTT" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > DTT-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=DTM" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > DTM-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=DTA" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > DTA-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=DTC" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > DTC-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=DTH" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > DTH-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=DXX" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > DXX-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=DTX" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > DTX-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=DHX" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > DHX-neighbor-TE-downstream.bed &&
awk -F "\t" '$10~ "=TXX" {print}' ./Neighbor_TEs_downstream.figure.bed | awk -F "\t" '{print $12 }' > TXX-neighbor-TE-downstream.bed &&
awk -F "\t" '{print $12 }' ./Neighbor_TEs_downstream.figure.bed > TOTAL-neighbor-TE-downstream.bed

##Sort overlap upstream TEs

awk -F "\t" '{$12 = ($3-$2)} 1' OFS='\t' ../Results/Annotations_TEs/Neighbor_TEs_upstream.bed > Neighbor_TEs_upstream.figure.bed &&  
awk -F "\t" '$10~ "=RLC" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > RLC-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=RLG" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > RLG-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=RYX" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > RYX-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=Caulimoviridae" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > Caulimoviridae-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=RIX" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > RIX-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=RLX" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > RLX-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=RXX" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > RXX-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=RSX" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > RSX-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=DTT" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > DTT-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=DTM" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > DTM-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=DTA" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > DTA-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=DTC" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > DTC-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=DTH" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > DTH-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=DXX" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > DXX-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=DTX" {print}'  ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > DTX-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=DHX" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > DHX-neighbor-TE-upstream.bed &&
awk -F "\t" '$10~ "=TXX" {print}' ./Neighbor_TEs_upstream.figure.bed | awk -F "\t" '{print $12 }' > TXX-neighbor-TE-upstream.bed &&
awk -F "\t" '{print $12 }' ./Neighbor_TEs_upstream.figure.bed > TOTAL-neighbor-TE-upstream.bed &&


find -empty -type f -delete
