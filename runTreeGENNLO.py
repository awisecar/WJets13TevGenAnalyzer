import os
import sys
#from ROOT import *
import random
import time
import datetime

### Objective: Make trees out of GEN-level info

cwd = os.getcwd()
print 'Current working directory: ' + cwd + '\n'

dateTo = datetime.datetime.now().strftime("%Y_%m_%d_%H%M%S")
txtdir = 'jobOutputGENTree_NLO_' + dateTo
os.system('mkdir -p ' + txtdir)

cmsswdir = '/afs/cern.ch/user/a/awisecar/WJetsGenMaker16/CMSSW_8_0_20/src'
cfgdir = '/afs/cern.ch/user/a/awisecar/WJetsGenMaker16/CMSSW_8_0_20/src/WJets13TeV/WJets13TevGenAnalyzer'

job = '#!/bin/bash \n'
job += 'cd ' + cwd + ' \n'
#job += 'eval `scramv1 runtime -sh` \n'

################
njobs = 20
doJet = [0]
#doJet = [2]
doVariation = [108, 110, 112, 114, 117, 119, 122, 124]
#doVariation = [112, 124]
################

for jet in doJet:
    for variation in doVariation:
        
        ### names for files and directories
        cfgFileName = 'wjetsGENTreeMaker_batch_cfg.py'
        
        jind = 0
        for i in range(0, njobs):
            print 'Job submission number #' + str(i)
            jind += 1
            
            ### build job script and outfile
            ajob = str(job)
            tjobname = txtdir + '/job_'+str(jet)+'j_0'+str(variation)+'_'+str(jind)+'.sh'
            tjobname_out = txtdir + '/job_'+str(jet)+'j_0'+str(variation)+'_'+str(jind)+'.out'
            tjob = open(tjobname,'w')
            tjob.write(ajob+'\n')
            
            ### commands to write to job script

            comCd = 'cd ' + cmsswdir + ' \n'
            comArch = 'export SCRAM_ARCH=slc6_amd64_gcc530 \n'
            comCmsenv = 'eval `scramv1 runtime -sh` \n'
            comCmsrun = 'cmsRun ' + cfgdir + '/' + cfgFileName + ' JetNum=jetnum VariationNum=variationnum OutFileNum=oFileNumIn 2>&1 \n\n'
            
            ### renaming
            cmsrun = comCmsrun.replace('oFileNumIn', str(jind))
            cmsrun = cmsrun.replace('jetnum', str(jet))
            cmsrun = cmsrun.replace('variationnum', str(variation))

            ### write out rest of job
            tjob.write(comCd)
            tjob.write(comArch)
            tjob.write(comCmsenv)
            tjob.write(cmsrun)
            
            tjob.close()
            os.system('chmod 755 ' + tjobname)
            bsub = 'bsub -q 1nd  -o ' + tjobname_out + ' -J ' +  tjobname + ' < ' + tjobname + ' '
            print bsub
            
            os.system(bsub)
            os.system('sleep 1')

