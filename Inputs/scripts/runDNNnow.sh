#!/bin/bash

YEAR=$1

./NohupParser.sh SingleMuon $YEAR 
./NohupParser.sh EmbeddedMuTau $YEAR 
./NohupParser.sh DYJets $YEAR 
./NohupParser.sh WJets $YEAR 
./NohupParser.sh TTbar $YEAR 
./NohupParser.sh Diboson $YEAR 
./NohupParser.sh GluGluHToUncorrTauTau $YEAR 
./NohupParser.sh VBFHToUncorrTauTau $YEAR 

