import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/work/a/awisecar/MadGraph5Gen/genOutputTotal_NLO/WJToLNu_0j_5f_NLO_FxFx_as_0124_cff_py_LHE_GEN_tempEnd.root'
        'file:/afs/cern.ch/user/a/awisecar/WJToLNu_1j_5f_NLO_FxFx_as_0108_GEN.root'
									  
    )
)

process.demo = cms.EDAnalyzer('WJets13TevAnalyzer')




process.printEventNumber = cms.OutputModule("AsciiOutputModule")
process.TFileService = cms.Service("TFileService",
								   fileName = cms.string('wjets_amcnlo_gen_DATE.root')
								   )

process.p = cms.Path(process.demo)
