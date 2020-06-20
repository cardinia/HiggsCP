#!/bin/sh 
# $1 - era
# $2 - category 
# $3 - channel 

cat > createcards$1_$2_$3.sh <<EOF1
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
cd ${CMSSW_BASE}/src
cmsenv
cd -
CreateCards card_config_etau_$1.conf $2 $3
EOF1

chmod u+x createcards$1_$2_$3.sh
nohup ./createcards$1_$2_$3.sh > nohup_$1_$2_$3.out &
