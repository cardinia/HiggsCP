#!/bin/bash

YEAR=$1
CHANNEL=$2

if [ "$CHANNEL" == "mt" ]; then
    echo "mutau channel"
    #./NohupParser.sh SingleMuon $YEAR $CHANNEL 
    #./NohupParser.sh EmbeddedMuTau $YEAR $CHANNEL 
else
    echo "etau channel"
    #./NohupParser.sh SingleElectron $YEAR $CHANNEL 
    #./NohupParser.sh EmbeddedElTau $YEAR $CHANNEL 
fi
#./NohupParser.sh DYJets $YEAR $CHANNEL 
#./NohupParser.sh WJets $YEAR $CHANNEL 
#./NohupParser.sh TTbar $YEAR $CHANNEL 
#./NohupParser.sh Diboson $YEAR $CHANNEL 
#./NohupParser.sh EWKZ $YEAR $CHANNEL 
#./NohupParser.sh HToWW $YEAR $CHANNEL 
#./NohupParser.sh GluGluHToUncorrTauTau $YEAR $CHANNEL 
./NohupParser.sh VBFHToUncorrTauTau $YEAR $CHANNEL 
#./NohupParser.sh WHToUncorrTauTau $YEAR $CHANNEL 
#./NohupParser.sh ZHToUncorrTauTau $YEAR $CHANNEL 

