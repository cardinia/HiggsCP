#!/bin/sh 
# $1 - sample
# $2 - era 

cat > $1_$2.submit <<EOF
+RequestRuntime=10000

RequestMemory = 2000

executable = $1_$2.sh

transfer_executable = True
universe            = vanilla
getenv              = True
Requirements        = OpSysAndVer == "SL6"

output              = createCards$1_$2_$3.out
error               = createCards$1_$2_$3.error
log                 = createCards$1_$2_$3.log

queue

EOF

cat > createCards$1_$2_$3.sh <<EOF1
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
cd ${CMSSW_BASE}/src
cmsenv
cd -
CreateCards $1 $2 $3
EOF1

chmod u+x createCards$1_$2_$3.sh
chmod u+x createCards$1_$2_$3.submit
condor_submit createCards$1_$2_$3.submit