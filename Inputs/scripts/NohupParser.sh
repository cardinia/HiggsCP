#!/bin/sh 
# $1 - sample
# $2 - era 


cat > $1_$2.sh <<EOF1
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
cd ${CMSSW_BASE}/src
cmsenv
cd -
CreateDNN $1 $2 mt
EOF1

chmod u+x $1_$2.sh
nohup ./$1_$2.sh > nohup_$1_$2.out &
