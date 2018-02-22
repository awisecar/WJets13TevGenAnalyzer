from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'WJets13Tev_amcnlo_Syst'
config.General.transferOutputs = True
config.General.transferLogs = True


config.JobType.pluginName = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName = 'wjets13tevanalyzer_cfg.py'


# 8 TeV sample
#config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_8TeV-amcatnloFXFX-pythia8/Summer12DR53X-PU_S10_TuneCUETP8M1_START53_V19-v3/AODSIM'

# large sample
#config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext2-v1/MINIAODSIM'
# small sample
config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'



#config.Data.inputDBS = 'global' # The default is global.
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 300
#config.Data.totalUnits = 2
#config.Data.outLFNDirBase = '/store/user/ahortian/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'WJets13Tev_amcnlo_Syst'


# This string is used to construct the output dataset name
#config.Data.outputDatasetTag = 'CRAB3_Analysis_test1'

# These values only make sense for processing data
#    Select input data based on a lumi mask
#config.Data.lumiMask = 'Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#    Select input data based on run-ranges
#config.Data.runRange = '190456-194076'

# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'
