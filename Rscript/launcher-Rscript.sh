#!/bin/bash

echo "launch script-TE-global-stat.r " 
Rscript ../Rscript/script-TE-global-stat.r &&
echo "OK"

echo "launch script-intragenic-stat.r " 
Rscript ../Rscript/script-intragenic-stat.r &&
echo "OK"

echo "launch script-intragenic-stat-length.r " 
Rscript ../Rscript/script-intragenic-stat-length.r &&
echo "OK"

echo " launch script-overlap-stat-family.r "
Rscript ../Rscript/script-overlap-stat-family.r &&
echo "OK"

#echo " launch script-overlap-stat-overlap-lenght.r "
#Rscript ../Rscript/script-overlap-stat-overlap-lenght.r &&
#echo "OK"

echo " launch script-overlap-stat-TElength.r "
Rscript ../Rscript/script-overlap-stat-TElength.r && 
echo "OK"

#echo " launch scrip-neighbor-TE-downstream.r "
#Rscript ../Rscript/scrip-neighbor-TE-downstream.r && 
#echo "OK"

#echo " launch scrip-neighbor-TE-upstream.r "
#Rscript ../Rscript/scrip-neighbor-TE-upstream.r
#echo "OK"

echo " launch scrip-neighbor-stat-family.r "
Rscript ../Rscript/script-neighbor-stat-family.r
echo "OK"

echo " launch scrip-neighbor-stat-length.r "
Rscript ../Rscript/script-neighbor-stat-length.r
echo "OK"
