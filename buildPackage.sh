#!/bin/bash

echo "Building Package..."

cd include

rm *.so *.d
echo HiggsCSandWidth
root -l -b -q buildHiggsCSandWidth.C
echo HiggsCSandWidthFermi
root -l -b -q buildHiggsCSandWidthFermi.C 
echo HiggsCSandWidthSM4
root -l -b -q buildHiggsCSandWidthSM4.C    
cd ../

echo "Done"
echo "Happy Hunting!"

