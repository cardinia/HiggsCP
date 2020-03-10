#!/bin/bash

YEAR=$1

./HTC_submit.sh SingleMuon $YEAR 10000
./HTC_submit.sh EmbeddedMuTau $YEAR 20000
./HTC_submit.sh DYJets $YEAR 10000
./HTC_submit.sh WJets $YEAR 10000
./HTC_submit.sh TTbar $YEAR 10000
./HTC_submit.sh Diboson $YEAR 10000
./HTC_submit.sh GluGluHToUncorrTauTau $YEAR 10000
./HTC_submit.sh VBFHToUncorrTauTau $YEAR 10000

