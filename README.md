# HiggsCP

Instructions for setup:

```
export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/${USER}/.cmsgit-cache
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700

cmsrel CMSSW_10_2_16

cd CMSSW_10_2_16/src
cmsenv

git cms-init

git clone https://github.com/CMS-HTT/Jet2TauFakes.git HTTutilities/Jet2TauFakes
cd HTTutilities/Jet2TauFakes
git checkout v0.2.2
scram b
mkdir data
git clone -b 2017 ssh://git@gitlab.cern.ch:7999/cms-htt/Jet2TauFakesFiles.git data
```
