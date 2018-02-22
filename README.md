# Instructions

--- set up working directory

mkdir WJetsGenMaker16

cd WJetsGenMaker16

source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH=slc6_amd64_gcc493

cmsrel CMSSW_8_0_20

cd CMSSW_8_0_20/src

cmsenv

voms-proxy-init -voms cms // to acsess file in das

--------------------------------------------------------------------------------

-- Clone /WJets13TevGenAnalyzer from the github

// an analyser must be in a subdirectory inside /src

// otherwise it will not be compiled 

mkdir WJets13TeV

cd WJets13TeV

git clone https://github.com/awisecar/WJets13TevGenAnalyzer.git

--------------------------------------------------------------------------------
-- Go to the directory

cd WJets13TevGenAnalyzer

--------------------------------------------------------------------------------
-- Run

scram b

cmsRun wjets13tevanalyzer_cfg.py

--------------------------------------------------------------------------------
-- for submitting crab jobs

source /cvmfs/cms.cern.ch/crab3/crab.sh

crab submit --config=crabConfig_WJetsAMCNLO.py

