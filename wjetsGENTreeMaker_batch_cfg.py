#Used for making trees out of GEN-level information
#See accompanying batch submission script
#A. Wisecarver - 26.3.2018

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

opt = VarParsing.VarParsing ('analysis')
opt.register('JetNum', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, 'input jet multiplicity')
opt.register('VariationNum', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, 'input alpha-s variation')
opt.register('OutFileNum', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, 'label output file number')
opt.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(20))

inputDir = '/store/user/awisecar/WJets_2018.02.26/GENoutput-NLO/'
fileDir = 'WJToLNu_MG5_NLO_FxFx_' + str(opt.JetNum) + 'j_as_0' + str(opt.VariationNum) + '/'
fileName = 'WJToLNu_' + str(opt.JetNum) + 'j_5f_NLO_FxFx_as_0' + str(opt.VariationNum) + '_GEN_' + str(opt.OutFileNum) + '.root'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         inputDir + fileDir + fileName
    )
)

process.demo = cms.EDAnalyzer('WJets13TevAnalyzer')

process.printEventNumber = cms.OutputModule("AsciiOutputModule")

treeDir = '/eos/cms/store/user/awisecar/WJets_2018.02.26/GENtree-NLO/'
outputFileName = 'WJToLNu_MG5_NLO_FxFx_' + str(opt.JetNum) + 'j_as_0' + str(opt.VariationNum) + '_' + str(opt.OutFileNum) + '.root'

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
        treeDir + fileDir + outputFileName
    )
)

process.p = cms.Path(process.demo)
