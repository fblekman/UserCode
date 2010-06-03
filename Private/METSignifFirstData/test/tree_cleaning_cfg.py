import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

process = cms.Process('NTUPLING')

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration/StandardSequences/MagneticField_cff")
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag ='GR_R_35X_V8B::All'


### PAT:
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.metTools import *

process.load('Configuration.StandardSequences.Services_cff')
process.add_( cms.Service( "TFileService",
fileName = cms.string("NtupleOutput.root"),
                           closeFileFast = cms.untracked.bool(True)  ) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source (
    "PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Commissioning10/MinimumBias/RAW-RECO/May6thPDSkim_GOODCOLL-v1/0003/26E87BA3-D05C-DF11-8FE0-001BFCDBD100.root'
        ),
    
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    
    secondaryFileNames = cms.untracked.vstring())

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.default.limit = 100

# summary
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# jet corrections, does not work in 358patch3
#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

# PhysicsDeclared filter
process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

# BPTX & BSC triggers filter
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

# Primary vertex filter
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"),
    filter = cms.bool(True)
)

# Scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

######################################################
####### MET CLEANING RECOMMENDATION STARTS HERE#######
######################################################

#is it MC or DATA
isMC = False
useHBHEcleaning = True
useHBHEfilter = True

HFPMTcleaningversion = 4   # version 1 = (loose), version 2 = (medium), version 3 = (tight)
# VERSION 4 is the currently recommended version, as of 28 May 2010.

process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
if useHBHEfilter == True:
    process.p = cms.Path(process.HBHENoiseFilter)
    
# New SeverityLevelComputer that forces RecHits with UserDefinedBit0 set to be excluded from new rechit collection
import JetMETAnalysis.HcalReflagging.RemoveAddSevLevel as RemoveAddSevLevel
process.hcalRecAlgos=RemoveAddSevLevel.RemoveFlag(process.hcalRecAlgos,"HFLongShort")

# UserDefinedBit0 is used by both the HF and HBHE reflaggers
process.hcalRecAlgos=RemoveAddSevLevel.AddFlag(process.hcalRecAlgos,"UserDefinedBit0",10)

# HF RecHit reflagger
process.load("JetMETAnalysis/HcalReflagging/HFrechitreflaggerJETMET_cff")
if HFPMTcleaningversion==1:
    process.hfrecoReflagged = process.HFrechitreflaggerJETMETv1.clone()
elif HFPMTcleaningversion==2:
    process.hfrecoReflagged = process.HFrechitreflaggerJETMETv2.clone()
elif HFPMTcleaningversion==3:
    process.hfrecoReflagged = process.HFrechitreflaggerJETMETv3.clone()
elif HFPMTcleaningversion==4:
    if (isMC==False):
        process.hfrecoReflagged = process.HFrechitreflaggerJETMETv4.clone()
    else:
        process.hfrecoReflagged = process.HFrechitreflaggerJETMETv2.clone()
elif HFPMTcleaningversion==5:
    if (isMC==False):
        process.hfrecoReflagged = process.HFrechitreflaggerJETMETv5.clone()
    else:
        process.hfrecoReflagged = process.HFrechitreflaggerJETMETv3.clone()


# HBHE RecHit reflagger
process.load("JetMETAnalysis/HcalReflagging/hbherechitreflaggerJETMET_cfi")
process.hbherecoReflagged = process.hbherechitreflaggerJETMET.clone()
process.hbherecoReflagged.debug=0

# Use the reflagged HF RecHits to make the CaloTowers
process.towerMaker.hfInput = "hfrecoReflagged"
process.towerMakerWithHO.hfInput = "hfrecoReflagged"

# Path and EndPath definitions

if (useHBHEcleaning==False):
    process.reflagging_step = cms.Path(process.hfrecoReflagged)
else:
    process.reflagging_step = cms.Path(process.hfrecoReflagged+process.hbherecoReflagged)
    # Need to specify that new HBHE collection should be fed to calotower maker
    process.towerMaker.hbheInput = "hbherecoReflagged"
    process.towerMakerWithHO.hbheInput = "hbherecoReflagged"

# Instead of rejecting the event, add a flag indicating the HBHE noise 
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
process.hbheflag = cms.Path(process.HBHENoiseFilterResultProducer)

######################################################
####### MET CLEANING RECOMMENDATION ENDS HERE#######
######################################################


process.calometntuple= cms.EDAnalyzer(
    "METNtupler",
    usePAT=cms.untracked.bool(True),
    src = cms.InputTag("met"),
	jetsrc = cms.InputTag("ak5CaloJets"),
	JetIDParams  = cms.PSet(useRecHits      = cms.bool(True),
							hbheRecHitsColl = cms.InputTag("hbhereco"),
							hoRecHitsColl   = cms.InputTag("horeco"),
							hfRecHitsColl   = cms.InputTag("hfreco"),
							ebRecHitsColl   = cms.InputTag("ecalRecHit", "EcalRecHitsEB"),
							eeRecHitsColl   = cms.InputTag("ecalRecHit", "EcalRecHitsEE")
							)
    )


process.pfmetntuple = process.calometntuple.clone()
process.pfmetntuple.usePAT=False
process.pfmetntuple.src = cms.InputTag("pfMet::RECO")
process.pfmetntuple.usePF = cms.untracked.bool(True)

process.combinedNtuple = cms.EDFilter("CombinedNtupler",
                                      filter = cms.untracked.bool(False),
                                      isMC		    = cms.bool(True),
                                      genmetTag	    = cms.InputTag("genMetCalo"),
                                      calometTag	    = cms.InputTag("met"), 
                                      pfmetTag	    = cms.InputTag("pfMet"), 
                                      tcmetTag	    = cms.InputTag("tcMet"), 
                                      calometoldTag	    = cms.untracked.InputTag("met","","RECO"), 
                                      pfmetoldTag	    = cms.untracked.InputTag("pfMet","","RECO"), 
                                      
                                      t1met05Tag	    = cms.untracked.InputTag("t1met05"), 
                                      t1met10Tag	    = cms.untracked.InputTag("t1met10"), 
                                      t1met20Tag	    = cms.untracked.InputTag("t1met20"), 
                                      t2met05Tag	    = cms.untracked.InputTag("t2met05"), 
                                      t2met10Tag	    = cms.untracked.InputTag("t2met10"), 
                                      t2met20Tag	    = cms.untracked.InputTag("t2met20"), 
                                      
                                      mht05Tag	    = cms.untracked.InputTag("mht05"), 
                                      mht10Tag	    = cms.untracked.InputTag("mht10"), 
                                      mht20Tag	    = cms.untracked.InputTag("mht20"), 
                                      
                                      mht05emfTag	    = cms.untracked.InputTag("mht05emf"), 
                                      mht10emfTag	    = cms.untracked.InputTag("mht10emf"), 
                                      mht20emfTag	    = cms.untracked.InputTag("mht20emf"), 
                                      
                                      l2l3mht05Tag	    = cms.untracked.InputTag("l2l3mht05"), 
                                      l2l3mht10Tag	    = cms.untracked.InputTag("l2l3mht10"), 
                                      l2l3mht20Tag	    = cms.untracked.InputTag("l2l3mht20"), 
                                      
                                      l2l3mht05newTag	    = cms.untracked.InputTag("l2l3mht05new"), 
                                      l2l3mht10newTag	    = cms.untracked.InputTag("l2l3mht10new"), 
                                      l2l3mht20newTag	    = cms.untracked.InputTag("l2l3mht20new"), 
                                      
                                      l2l3mht05emfTag	    = cms.untracked.InputTag("l2l3mht05emf"), 
                                      l2l3mht10emfTag	    = cms.untracked.InputTag("l2l3mht10emf"), 
                                      l2l3mht20emfTag	    = cms.untracked.InputTag("l2l3mht20emf"), 
                                      
                                      jptmht05emfTag	    = cms.untracked.InputTag("jptmht05emf"), 
                                      jptmht10emfTag	    = cms.untracked.InputTag("jptmht10emf"), 
                                      jptmht20emfTag	    = cms.untracked.InputTag("jptmht20emf"), 
                                      
                                      electronTag	    = cms.untracked.InputTag("gsfElectrons"),
                                      patElectronTag        = cms.InputTag("cleanPatElectrons"),
                                      jetTag		    = cms.InputTag("ak5CaloJets"),
                                      patJetTag             = cms.InputTag("cleanPatJets"),
                                      jetIDTag	    = cms.InputTag("ak5JetID"),
                                      TriggerResultsTag   = cms.untracked.InputTag("TriggerResults","","HLT"),
                                      TriggerPath	    = cms.untracked.string("HLT_Ele15_LW_L1R"),
###                                    jetCorrectorTag	    = cms.untracked.string("ak5CaloL2L3"), # not avaialble in the release in 258patch3
                                      genjetTag	    = cms.InputTag("ak5GenJets"),
                                      TriggerNamesForNtuple = cms.untracked.vstring("HLT_Ele15_LW_L1R","HLT_Ele15_SW_L1R","HLT_Photon10_L1R","HLT_MinBiasBSC","HLT_Jet15U")
                                      )
# use the value of isMC as set in this file:
process.combinedNtuple.isMC = isMC

process.cleanupFilter = cms.Path(
    process.hltLevel1GTSeed*
    process.primaryVertexFilter*
    process.scrapingVeto*
    process.HBHENoiseFilter
    )

process.filterSequence = cms.Sequence(

    process.hltLevel1GTSeed*
    process.primaryVertexFilter*
    process.scrapingVeto*
    process.HBHENoiseFilter
)

# pool output module, uses above filter:

#control event content
from PhysicsTools.PatAlgos.patEventContent_cff import *
if (isMC==False):
    removeMCMatching(process, ['All'])
    
process.out = cms.OutputModule("PoolOutputModule",
                               # additional content is added below this module def
                               outputCommands = cms.untracked.vstring("drop *","keep *_pfMet*_*_*","keep *_particleFlow_*_*"),
                               dataset = cms.untracked.PSet(dataTier = cms.untracked.string('RAW-RECO'),
                                                            filterName = cms.untracked.string('cleanupFilter')
                                                            ),
                               
                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('cleanupFilter')
                                                                 )
                               )
process.out.outputCommands += patEventContent 
process.out.outputCommands += patTriggerEventContent
process.out.outputCommands += patExtraAodEventContent

process.out.fileName = cms.untracked.string('EDMOutput.root')


# Path and EndPath definitions
process.rereco_step = cms.Path(process.filterSequence*process.caloTowersRec*(process.recoJets*process.recoJetIds+process.recoTrackJets)*process.recoJetAssociations*process.btagging*process.metreco) # re-reco jets and MET


process.fix_and_pat=cms.Path(
    process.filterSequence*
    # not all PF stuff is present to run this one
#    process.highlevelreco*
    process.recoJetAssociations*
    process.btagging*
    process.patDefaultSequence
    )


process.ntuple_step=cms.Path(
    process.filterSequence*
    (
    process.calometntuple*
    process.pfmetntuple
    *process.combinedNtuple
    
    )
)

process.out_step = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.reflagging_step, process.rereco_step , process.hbheflag, process.fix_and_pat,process.ntuple_step, process.cleanupFilter, process.out_step)
