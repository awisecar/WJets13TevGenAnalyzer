#!/bin/bash 
cd /afs/cern.ch/user/a/awisecar/WJetsGenMaker16/CMSSW_8_0_20/src/WJets13TeV/WJets13TevGenAnalyzer 

cd /afs/cern.ch/user/a/awisecar/WJetsGenMaker16/CMSSW_8_0_20/src 
export SCRAM_ARCH=slc6_amd64_gcc530 
eval `scramv1 runtime -sh` 
cmsRun /afs/cern.ch/user/a/awisecar/WJetsGenMaker16/CMSSW_8_0_20/src/WJets13TeV/WJets13TevGenAnalyzer/wjetsGENTreeMaker_batch_cfg.py JetNum=0 VariationNum=124 OutFileNum=4 2>&1 

