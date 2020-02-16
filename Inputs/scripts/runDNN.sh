#!/bin/bash

YEAR=$1

./HTC_submit.sh SingleMuon $YEAR
./HTC_submit.sh EmbeddedMuTau $YEAR
./HTC_submit.sh DYJets $YEAR
./HTC_submit.sh WJets $YEAR
./HTC_submit.sh TTbar $YEAR
./HTC_submit.sh SingleTop $YEAR
./HTC_submit.sh Diboson $YEAR
./HTC_submit.sh GluGluHToUncorrTauTau $YEAR
./HTC_submit.sh VBFHToUncorrTauTau $YEAR

