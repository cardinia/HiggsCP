#!/bin/sh 
# $1 - sample
# $2 - era 
# $3 - runtime

cat > $1_$2.submit <<EOF
+RequestRuntime=$3

RequestMemory = 2000

executable = $1_$2.sh

transfer_executable = True
universe            = vanilla
getenv              = True
Requirements        = OpSysAndVer == "SL6"

output              = $1_$2.out
error               = $1_$2.error
log                 = $1_$2.log

queue

EOF

cat > $1_$2.sh <<EOF1
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
cd ${CMSSW_BASE}/src
cmsenv
cd -
CreateDNN $1 $2 mt
EOF1

chmod u+x $1_$2.sh
chmod u+x $1_$2.submit
condor_submit $1_$2.submit