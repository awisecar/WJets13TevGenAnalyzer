import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
									  #'file:myfile.root'
									  # MLM
#									  '/store/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0C765598-8BD1-E511-BF63-20CF3027A566.root',
									  
									  # amcAtNLO small sample
									  '/store/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/00782B93-95BF-E511-86D4-0025905C446C.root'
									  # amcAtNLO
#									  '/store/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext2-v1/10000/002B75B9-10DB-E511-915B-02163E017888.root',
#									  '/store/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext2-v1/10000/006E9C0C-12DB-E511-981C-00261894397D.root',
#									  '/store/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext2-v1/10000/0097C77B-12DB-E511-AAC0-02163E011796.root',
#									  '/store/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext2-v1/10000/00B697CF-10DB-E511-A9D9-782BCB407CFD.root',
#									  '/store/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext2-v1/10000/00BC9D60-0FDB-E511-A038-00215AD4D670.root',
#									  '/store/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext2-v1/10000/00D29FE6-0CDB-E511-8602-BCAEC53F6D4E.root',
									  
    )
)

process.demo = cms.EDAnalyzer('WJets13TevAnalyzer')




process.printEventNumber = cms.OutputModule("AsciiOutputModule")
process.TFileService = cms.Service("TFileService",
								   fileName = cms.string('wjets_amcnlo_gen.root')
								   )

process.p = cms.Path(process.demo)
