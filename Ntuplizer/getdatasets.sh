#!/bin/bash
d=$(date +%Y-%m-%d)
olddate=2022-08-01
for i in {0..5} 
do  
    mkdir
    dasgoclient --query "file dataset=/ParkingDoubleElectronLowMass"$i"/Run2022C-PromptReco-v1/MINIAOD instance=prod/global" > batchjobsautomation/stream$i/completedatafiles"_"$d/ParkingDoubleElectronLowMass$i"_"$d"_"new.txt
done

for i in {0..5}
do
    grep -Fxvf ParkingDoubleElectronLowMass$i"_"$olddate.txt ParkingDoubleElectronLowMass$i"_"$d"_"new.txt >  batchjobsautomation/stream$i/newfiles$i.txt
done