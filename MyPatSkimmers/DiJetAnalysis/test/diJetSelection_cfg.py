from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *


useData=False

###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

import sys
 
# get the 7 TeV jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *


if useData == False :
    # Make sure to NOT apply L2L3Residual to MC
    corrections = ['L2Relative', 'L3Absolute']
    # global tag for 384 MC
    process.GlobalTag.globaltag = cms.string('START311_V2::All')
else :
    # Make sure to apply L2L3Residual to data
    corrections = ['L2Relative', 'L3Absolute', 'L2L3Residual']
    # global tag for 386 data
    process.GlobalTag.globaltag = cms.string('GR_R_311_V2::All')


# add the flavor history
process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryPaths_cfi")

# PF from RECO and not using PF2PAT
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF',corrections),
                 doType1MET   = False,
                 doL1Cleaning = False,
                 doL1Counters = False,                 
                 genJetCollection = cms.InputTag("ak5GenJets"),
                 doJetID      = False,
                 jetIdLabel   = "ak5"
                 )

from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')


process.selectedPatJetsAK5PF.cut = cms.string('pt > 20 & abs(eta) < 2.4')
process.selectedPatJets.cut = cms.string('pt > 20 & abs(eta) < 2.4')

process.patJets.tagInfoSources = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAOD")
    )
process.patJetsAK5PF.tagInfoSources = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAK5PF")
    )

process.selectedPatMuons.cut = cms.string("pt > 10")
process.patMuons.usePV = False
process.patMuons.embedTrack = True

process.selectedPatElectrons.cut = cms.string("pt > 10")
process.patElectrons.usePV = False
process.patElectrons.embedTrack = True


process.patJets.userData.userFunctions = cms.vstring( "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
                                                      "tagInfoSecondaryVertex('secondaryVertex').secondaryVertex(0).p4().mass() : 0")
process.patJets.userData.userFunctionLabels = cms.vstring('secvtxMass')


process.patJetsAK5PF.userData.userFunctions = cms.vstring( "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
                                                      "tagInfoSecondaryVertex('secondaryVertex').secondaryVertex(0).p4().mass() : 0")
process.patJetsAK5PF.userData.userFunctionLabels = cms.vstring('secvtxMass')

process.patJets.embedPFCandidates = True
process.patJets.embedCaloTowers = True
process.patJetsAK5PF.embedCaloTowers = True
process.patJetsAK5PF.embedPFCandidates = True

# prune gen particles
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


process.prunedGenParticles = cms.EDProducer("GenParticlePruner",
                                            src = cms.InputTag("genParticles"),
                                            select = cms.vstring(
                                                "drop  *",
                                                "keep status = 3", #keeps all particles from the hard matrix element
                                                "+keep (abs(pdgId) = 11 | abs(pdgId) = 13) & status = 1" #keeps all stable muons and electrons and their (direct) mothers.
                                                )
                                            )

process.patseq = cms.Sequence(
#    process.scrapingVeto*
#    process.primaryVertexFilter*
#    process.HBHENoiseFilter*
    process.patDefaultSequence* 
    process.flavorHistorySeq *
    process.prunedGenParticles
)

process.p1 = cms.Path(
    process.patseq
    )

process.out.SelectEvents.SelectEvents = cms.vstring('p1')

process.maxEvents.input = -1
process.options.wantSummary = True
process.out.dropMetaData = cms.untracked.string("DROPPED")

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent

process.out.outputCommands = [
    'keep GenRunInfoProduct_generator_*_*',
    'keep GenEventInfoProduct_generator_*_*',
    'keep *_flavorHistoryFilter_*_*',
    'keep *_prunedGenParticles_*_*',
#    'keep *_decaySubset_*_*',
#    'keep *_initSubset_*_*',
    'drop *_cleanPat*_*_*',
    'keep *_selectedPat*_*_*',
    'keep *_patMETs*_*_*',
#    'keep recoPFCandidates_particleFlow_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_offlinePrimaryVertices_*_*',
    'keep recoTracks_generalTracks_*_*',
    'drop patPFParticles_*_*_*',
#    'keep patTriggerObjects_patTrigger_*_*',
#    'keep patTriggerFilters_patTrigger_*_*',
    'keep patTriggerPaths_patTrigger_*_*',
    'keep patTriggerEvent_patTriggerEvent_*_*',
#    'keep *_cleanPatPhotonsTriggerMatch_*_*',
#    'keep *_cleanPatElectronsTriggerMatch_*_*',
#    'keep *_cleanPatMuonsTriggerMatch_*_*',
#    'keep *_cleanPatTausTriggerMatch_*_*',
#    'keep *_cleanPatJetsTriggerMatch_*_*',
#    'keep *_patMETsTriggerMatch_*_*',
    'drop *_MEtoEDMConverter_*_*'
    ]
