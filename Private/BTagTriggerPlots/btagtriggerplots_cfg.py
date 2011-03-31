import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
                                'file:jet2011A_aod.root'

    )
)

process.trigMuTwentyRelToDijetTwenty = cms.EDAnalyzer('BTagTriggerPlots',
                              jetSrc=cms.string("selectedPatJets"),
                              muSrc=cms.string("selectedPatMuons"),
                              trigSrc=cms.string("hltTriggerSummaryAOD"),
                              trigName=cms.string("hltBDiJet20Central"),
                               monTrigName=cms.string("hltBSoftMuonDiJet20Mu5SelL3FilterByDR")
                               )

process.trigMuSixtyRelToDijetTwenty = process.trigMuTwentyRelToDijetTwenty.clone()
process.trigMuSixtyRelToDijetTwenty.monTrigName=cms.string("hltBSoftMuonDiJet60Mu7SelL3FilterByDR")
process.trigMuSixtyRelToMuTwenty = process.trigMuSixtyRelToDijetTwenty.clone()
process.trigMuSixtyRelToMuTwenty.trigName=cms.string("hltBSoftMuonDiJet20Mu5SelL3FilterByDR")
process.trigMuEightyRelToDijetTwenty = process.trigMuTwentyRelToDijetTwenty.clone()
process.trigMuEightyRelToDijetTwenty.monTrigName=cms.string("hltBSoftMuonDiJet80Mu9SelL3FilterByDR")
process.trigMuHundredRelToDijetTwenty = process.trigMuTwentyRelToDijetTwenty.clone()
process.trigMuHundredRelToDijetTwenty.monTrigName=cms.string("hltBSoftMuonDiJet100Mu9SelL3FilterByDR")
process.trigMuEightyReltoMuTwenty = process.trigMuSixtyRelToMuTwenty.clone()
process.trigMuEightyReltoMuTwenty.monTrigName=cms.string("hltBSoftMuonDiJet80Mu9SelL3FilterByDR")
process.trigMuHundredReltoMuTwenty = process.trigMuSixtyRelToMuTwenty.clone()
process.trigMuHundredReltoMuTwenty.monTrigName=cms.string("hltBSoftMuonDiJet100Mu9SelL3FilterByDR")




process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('ntuple.root')
                                   )

process.p = cms.Path(process.trigMuHundredRelToDijetTwenty*process.trigMuTwentyRelToDijetTwenty*process.trigMuSixtyRelToDijetTwenty*process.trigMuEightyRelToDijetTwenty*process.trigMuSixtyRelToMuTwenty*process.trigMuEightyReltoMuTwenty*process.trigMuHundredReltoMuTwenty)
