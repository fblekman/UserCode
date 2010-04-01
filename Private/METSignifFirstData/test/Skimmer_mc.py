#
# using: 
# Revision: 1.149 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: promptCollisionReco -s RAW2DIGI,L1Reco,RECO,DQM,ALCA:SiStripCalZeroBias --datatier RECO --eventcontent RECO --conditions CRAFT09_R_V4::All --scenario pp --no_exec --data --magField AutoFromDBCurrent -n 100
import FWCore.ParameterSet.Config as cms

process = cms.Process('METSIGNIF')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/RawToDigi_Data_cff')
process.load('Configuration/StandardSequences/L1Reco_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('DQMOffline/Configuration/DQMOffline_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('Configuration/StandardSequences/AlCaRecoStreams_cff')
process.load('Configuration/EventContent/AlCaRecoOutput_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.4 $'),
    annotation = cms.untracked.string('rereco nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
	)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True) 
)
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/FED831D5-F03B-DF11-9C75-0030487D05B0.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/FED078E6-E63B-DF11-B17F-001D09F24763.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/FE290C8E-F83B-DF11-8C7D-001D09F2423B.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/FE28014B-E83B-DF11-9B06-001D09F24493.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/FCDF1B49-E83B-DF11-B789-001D09F2841C.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/FC7EB40D-F63B-DF11-8DFA-001D09F24EE3.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/FC2A72FC-E83B-DF11-AFE2-001D09F2432B.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/FC263E7C-F63B-DF11-B2E5-001D09F26509.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/FAF4FC91-EC3B-DF11-9E16-001D09F24664.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/FA9CAE79-F13B-DF11-A70B-001D09F251BD.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/FA0AF7DF-F23B-DF11-A48B-000423D9863C.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/F8171699-EC3B-DF11-B21C-001D09F28EC1.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/F6A93EE4-F23B-DF11-B0B3-001D09F24691.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/F64F0B1D-F73B-DF11-9535-000423D98B5C.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/F489FD08-FA3B-DF11-BB79-001D09F251CC.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/F2823533-FB3B-DF11-9CDA-000423D98AF0.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/F2348807-FA3B-DF11-B296-001D09F2924F.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/F030FEF3-ED3B-DF11-824C-000423D9989E.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/EE6BC8C3-F73B-DF11-8265-000423D99EEE.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/ECE1B065-EF3B-DF11-A910-001D09F28E80.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/ECA4C56A-EA3B-DF11-AF13-001D09F23A3E.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/EC911185-F13B-DF11-AC24-001D09F28F1B.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/EC632C91-F83B-DF11-A0E6-001D09F25208.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/EAA75F3C-FB3B-DF11-A7E8-001D09F251CC.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/EA83F1FE-E83B-DF11-A804-001D09F24EE3.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/EA0FD15F-EF3B-DF11-B32C-001D09F29114.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E88F90C6-EB3B-DF11-9703-001D09F231C9.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E885D5E3-F23B-DF11-9608-000423D986A8.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E8759D1F-F53B-DF11-B39D-001D09F24EAC.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E6BCCD7B-F63B-DF11-BC39-000423D6CA6E.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E68EB53A-FB3B-DF11-8007-001D09F24682.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E667695F-EF3B-DF11-8A0F-000423D991F0.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E65B1BFC-E83B-DF11-9686-001D09F28EA3.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E6407E33-FB3B-DF11-A9DB-0030487CD6E6.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E623F945-F93B-DF11-B532-000423D6CAF2.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E4CF833C-FB3B-DF11-8248-001D09F29524.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E496D7FC-E83B-DF11-A7D1-001D09F24DDF.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E43EF52D-ED3B-DF11-861E-001D09F2B30B.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E4311843-F93B-DF11-B640-000423D9870C.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E2D4F81A-EB3B-DF11-9AF8-001D09F2424A.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E0DE7DAD-E93B-DF11-90B0-001D09F24691.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E0C2B3E6-E63B-DF11-9623-001D09F2527B.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E0A64036-F43B-DF11-B031-001D09F29533.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E06B2608-F03B-DF11-AF3E-000423D999CA.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/E0592D2F-ED3B-DF11-98B0-001D09F252F3.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/DEA5480D-FA3B-DF11-A3B8-001D09F291D7.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/DE7A5F96-F33B-DF11-BB50-000423D6BA18.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/DE65B86A-EA3B-DF11-A290-001D09F251BD.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/DE568531-ED3B-DF11-B491-001D09F29597.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/DE5320FF-F53B-DF11-A816-0030487C8CB8.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/DCDAD660-EF3B-DF11-8EA1-001D09F24DA8.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/DCC7B2AC-E93B-DF11-9FD0-001D09F28D54.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/DC5B897B-F63B-DF11-8353-001617E30CD4.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/DACCE280-F13B-DF11-8678-000423D9517C.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/D8D58E37-D4632331-F23B-DF11-9EED-000423D6B48C.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/D43107C7-EB3B-DF11-8043-0015C5FDE067.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/D40DD232-ED3B-DF11-B229-0019B9F72BAA.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/D2F0F3B3-F03B-DF11-A360-001D09F27003.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/D2B0BA81-F13B-DF11-9B65-001D09F295FB.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/D2B04463-EA3B-DF11-8C1C-000423D99F1E.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/D27D619A-EC3B-DF11-83ng10/ExpressPhysics/FEVT/v7/000/132/440/D2037809-FA3B-DF11-A821-001D09F2527B.root',
	          '/storning10/ExpressPhysics/FEVT/v7/000/132/440/BA4106FF-F53B-DF11-BDB0-0019B9F709A4.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/BA0C4C6A-E83B-DF11-97E3-001D09F2512C.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/B81E8A81-F13B-DF11-BA9B-001D09F28EA3.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/B6F76493-EC3B-DF11-ACFA-001D09F253C0.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/B6674D6A-EA3B-DF11-899F-001D09F231C9.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/B661C0BE-E73B-DF11-9048-000423D98AF0.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/B60DFCAA-FA3B-DF11-880E-000423D174FE.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/B4FF1BE9-F23B-DF11-92FF-001D09F25456.root',
	          '/store/express/Commissioning10/ExpressPhysics/FEVT/v7/000/132/440/B4BB7EC2-EB3B-DF11-9653-001D09F2A690.root'
    )
)

			
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(
#   '123596:2-123596:9999',    
#   '123615:70-123615:9999',   
#   '123732:62-123732:9999',  
#   '123815:8-123815:9999',   
#   '123818:2-123818:42',     
#   '123908:2-123908:12',     
#   '124008:1-124008:1',      
#   '124009:1-124009:68',     
#   '124020:12-124020:94',     
#   '124022:66-124022:179',    
#   '124023:38-124023:9999',   
#   '124024:2-124024:83',     
#   '124025:5-124025:13',    
#   '124027:24-124027:9999',   
#   '124030:2-124030:9999',  
#   '124120:1-124120:9999',   
#   '124275:3-124275:30'
#   )

#process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*", "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT")

# Output definition
process.FEVT = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
#    outputCommands = process.RECOEventContent.outputCommands,
								outputCommands = cms.untracked.vstring("drop *","keep *_TriggerResults_*_*","keep *_hcalhoise_*_*",
																	   "keep *_hltL1GtObjectMap_*_*","keep *_gtDigis_*_*",
																	   "keep *_hltTriggerSummaryAOD_*_*","keep *_ak5CaloJets_*_*","keep *_ak5PFJets_*_*",
																	   "keep *_ecalRecHit_*_*","keep *_hbhereco_*_*","keep *_hfreco_*_*","keep *_horeco_*_*",
																	   "keep *_cosmicMuons_*_*","keep *_ctfPixelLess_*_*","keep *_generalTracks_*_*",
																	   "keep *_pixelTracks_*_*","keep *_towerMaker_*_*",
																	   "keep *_hfEMClusters_*_*","keep *_hybridSuperClusters_*_*",
																	   "keep *_met_*_*","keep *_tcMet_*_*","keep *_particleFlow_*_*",
																	   "keep *_particleFlowClusterECAL_*_*","keep *_particleFlowClusterHCAL_*_*",
																	   "keep *_pfMet_*_*","keep *_offlinePrimaryVertices_*_*",
																	   "keep *_impactParameterTagInfos_*_*","keep *_pixelVertices_*_*",
																	   "keep *_generalV0Candidates_*_*",
																	   # additional MC info:
																	   "keep *_generator_*_*","keep *_genParticles_*_*"),
    fileName = cms.untracked.string('SKIMMED_MC_withNewMETandGTdigis.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RECO'),
        filterName = cms.untracked.string('')
    )
)


# Include UserDefinedBit0 in HcalSeverityLevelComputer
process.hcalRecAlgos.SeverityLevels.append(cms.PSet(RecHitFlags = cms.vstring('UserDefinedBit0'),
                                                    ChannelStatus = cms.vstring(''),
                                                    Level = cms.int32(9998)))

# HF RecHit re-flagger
process.load("RecoLocalCalo.HcalRecAlgos.hcalrechitreflagger_cfi")
process.hfrecoReflagged = process.hcalrechitReflagger.clone()
#process.hfrecoReflagged.hf_Algo3test = True  # S9/S1 algorithm (current default)
#process.hfrecoReflagged.hf_Algo2test = False # PET algorithm

# Use the re-flagged HF RecHits to make the CaloTowers
process.towerMaker.hfInput = cms.InputTag("hfrecoReflagged")
process.towerMakerWithHO.hfInput = cms.InputTag("hfrecoReflagged")


# Other statements
process.GlobalTag.globaltag = 'START3X_V25B::All'


#####################################################################################################
####
####  Top level replaces for handling strange scenarios of early collisions
####

## TRACKING:
## Skip events with HV off
process.newSeedFromTriplets.ClusterCheckPSet.MaxNumberOfPixelClusters=2000
process.newSeedFromPairs.ClusterCheckPSet.MaxNumberOfCosmicClusters=10000
process.secTriplets.ClusterCheckPSet.MaxNumberOfPixelClusters=1000
process.fifthSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters = 5000
process.fourthPLSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters=10000

###### FIXES TRIPLETS FOR LARGE BS DISPLACEMENT ######

### pixelTracks
#---- replaces ----
process.pixelTracks.RegionFactoryPSet.ComponentName = 'GlobalRegionProducerFromBeamSpot' # was GlobalRegionProducer
process.pixelTracks.OrderedHitsFactoryPSet.GeneratorPSet.useFixedPreFiltering = True     # was False
#---- new parameters ----
process.pixelTracks.RegionFactoryPSet.RegionPSet.nSigmaZ  = cms.double(4.06) # was originHalfLength = 15.9; translated assuming sigmaZ ~ 3.8
process.pixelTracks.RegionFactoryPSet.RegionPSet.beamSpot = cms.InputTag("offlineBeamSpot")

### 0th step of iterative tracking
#---- replaces ----
process.newSeedFromTriplets.RegionFactoryPSet.ComponentName = 'GlobalRegionProducerFromBeamSpot' # was GlobalRegionProducer
process.newSeedFromTriplets.OrderedHitsFactoryPSet.GeneratorPSet.useFixedPreFiltering = True     # was False
#---- new parameters ----
process.newSeedFromTriplets.RegionFactoryPSet.RegionPSet.nSigmaZ   = cms.double(4.06)  # was originHalfLength = 15.9; translated assuming sigmaZ ~ 3.8
process.newSeedFromTriplets.RegionFactoryPSet.RegionPSet.beamSpot = cms.InputTag("offlineBeamSpot")

### 2nd step of iterative tracking
#---- replaces ----
process.secTriplets.RegionFactoryPSet.ComponentName = 'GlobalRegionProducerFromBeamSpot' # was GlobalRegionProducer
process.secTriplets.OrderedHitsFactoryPSet.GeneratorPSet.useFixedPreFiltering = True     # was False
#---- new parameters ----
process.secTriplets.RegionFactoryPSet.RegionPSet.nSigmaZ  = cms.double(4.47)  # was originHalfLength = 17.5; translated assuming sigmaZ ~ 3.8
process.secTriplets.RegionFactoryPSet.RegionPSet.beamSpot = cms.InputTag("offlineBeamSpot")

## Primary Vertex
process.offlinePrimaryVerticesWithBS.PVSelParameters.maxDistanceToBeam = 2
process.offlinePrimaryVerticesWithBS.TkFilterParameters.maxNormalizedChi2 = 20
process.offlinePrimaryVerticesWithBS.TkFilterParameters.minSiliconHits = 6
process.offlinePrimaryVerticesWithBS.TkFilterParameters.maxD0Significance = 100
process.offlinePrimaryVerticesWithBS.TkFilterParameters.minPixelHits = 1
process.offlinePrimaryVerticesWithBS.TkClusParameters.zSeparation = 10
process.offlinePrimaryVertices.PVSelParameters.maxDistanceToBeam = 2
process.offlinePrimaryVertices.TkFilterParameters.maxNormalizedChi2 = 20
process.offlinePrimaryVertices.TkFilterParameters.minSiliconHits = 6
process.offlinePrimaryVertices.TkFilterParameters.maxD0Significance = 100
process.offlinePrimaryVertices.TkFilterParameters.minPixelHits = 1
process.offlinePrimaryVertices.TkClusParameters.zSeparation = 10

## ECAL 
process.ecalRecHit.ChannelStatusToBeExcluded = [ 1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14, 78, 142 ]

##Preshower
process.ecalPreshowerRecHit.ESGain = 2
process.ecalPreshowerRecHit.ESBaseline = 0
process.ecalPreshowerRecHit.ESMIPADC = 55

##only for 34X
process.ecalPreshowerRecHit.ESRecoAlgo = cms.untracked.int32(1)

## HCAL temporary fixes
process.hfreco.firstSample  = 3
process.hfreco.samplesToAdd = 4

process.zdcreco.firstSample = 4
process.zdcreco.samplesToAdd = 3

## EGAMMA
process.ecalDrivenElectronSeeds.SCEtCut = cms.double(1.0)
process.ecalDrivenElectronSeeds.applyHOverECut = cms.bool(False)
process.ecalDrivenElectronSeeds.SeedConfiguration.z2MinB = cms.double(-0.9)
process.ecalDrivenElectronSeeds.SeedConfiguration.z2MaxB = cms.double(0.9)
process.ecalDrivenElectronSeeds.SeedConfiguration.r2MinF = cms.double(-1.5)
process.ecalDrivenElectronSeeds.SeedConfiguration.r2MaxF = cms.double(1.5)
process.ecalDrivenElectronSeeds.SeedConfiguration.rMinI = cms.double(-2.)
process.ecalDrivenElectronSeeds.SeedConfiguration.rMaxI = cms.double(2.)
process.ecalDrivenElectronSeeds.SeedConfiguration.DeltaPhi1Low = cms.double(0.3)
process.ecalDrivenElectronSeeds.SeedConfiguration.DeltaPhi1High = cms.double(0.3)
process.ecalDrivenElectronSeeds.SeedConfiguration.DeltaPhi2 = cms.double(0.3)
process.gsfElectrons.applyPreselection = cms.bool(False)
process.photons.minSCEtBarrel = 1.
process.photons.minSCEtEndcap =1.
process.photonCore.minSCEt = 1.
process.conversionTrackCandidates.minSCEt =1.
process.conversions.minSCEt =1.
process.trackerOnlyConversions.AllowTrackBC = cms.bool(False)
process.trackerOnlyConversions.AllowRightBC = cms.bool(False)
process.trackerOnlyConversions.MinApproach = cms.double(-.25)
process.trackerOnlyConversions.DeltaCotTheta = cms.double(.07)
process.trackerOnlyConversions.DeltaPhi = cms.double(.2)

###
###  end of top level replacements
###
###############################################################################################

# produce L1 trigger object maps (temporary fix for HLT mistake 
# in event content definition of RAW datatier for stream A)
import L1Trigger.GlobalTrigger.gtDigis_cfi
process.hltL1GtObjectMap = L1Trigger.GlobalTrigger.gtDigis_cfi.gtDigis.clone()
process.hltL1GtObjectMap.GmtInputTag = cms.InputTag( "gtDigis" )
process.hltL1GtObjectMap.GctInputTag = cms.InputTag( "gctDigis" )
process.hltL1GtObjectMap.CastorInputTag = cms.InputTag( "castorL1Digis" )
process.hltL1GtObjectMap.ProduceL1GtDaqRecord = cms.bool( False )
process.hltL1GtObjectMap.ProduceL1GtEvmRecord = cms.bool( False )
process.hltL1GtObjectMap.ProduceL1GtObjectMapRecord = cms.bool( True )
process.hltL1GtObjectMap.WritePsbL1GtDaqRecord = cms.bool( False )
process.hltL1GtObjectMap.ReadTechnicalTriggerRecords = cms.bool( True )
process.hltL1GtObjectMap.EmulateBxInEvent = cms.int32( 1 )
process.hltL1GtObjectMap.AlternativeNrBxBoardDaq = cms.uint32( 0 )
process.hltL1GtObjectMap.AlternativeNrBxBoardEvm = cms.uint32( 0 )
process.hltL1GtObjectMap.BstLengthBytes = cms.int32( -1 )
process.hltL1GtObjectMap.TechnicalTriggersInputTags = cms.VInputTag( 'simBscDigis' )
process.hltL1GtObjectMap.RecordLength = cms.vint32( 3, 0 )



# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1GtObjectMap_step = cms.Path(process.hltL1GtObjectMap)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction_withPixellessTk)
process.dqmoffline_step = cms.Path(process.DQMOffline)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.FEVT)


process.pathALCARECOHcalCalIsoTrk = cms.Path(process.seqALCARECOHcalCalIsoTrk*process.ALCARECOHcalCalIsoTrackDQM)
process.pathALCARECOTkAlMuonIsolated = cms.Path(process.seqALCARECOTkAlMuonIsolated*process.ALCARECOTkAlMuonIsolatedDQM)
process.pathALCARECOHcalCalDijets = cms.Path(process.seqALCARECOHcalCalDijets*process.ALCARECOHcalCalDiJetsDQM)
process.pathALCARECOMuAlCalIsolatedMu = cms.Path(process.seqALCARECOMuAlCalIsolatedMu*process.ALCARECOMuAlCalIsolatedMuDQM*process.ALCARECODTCalibrationDQM)
process.pathALCARECOSiStripCalZeroBias = cms.Path(process.seqALCARECOSiStripCalZeroBias*process.ALCARECOSiStripCalZeroBiasDQM)
process.pathALCARECOTkAlMinBias = cms.Path(process.seqALCARECOTkAlMinBias*process.ALCARECOTkAlMinBiasDQM)
process.pathALCARECOMuAlOverlaps = cms.Path(process.seqALCARECOMuAlOverlaps*process.ALCARECOMuAlOverlapsDQM)

process.reflagging_step = cms.Path(process.hfrecoReflagged)
process.rereco_step = cms.Path(process.caloTowersRec*(process.recoJets*process.recoJetIds+process.recoTrackJets)*process.recoJetAssociations*process.metreco) # re-

process.metreco_step = cms.Path(process.caloTowersRec+process.metreco+process.pfMet)

# Schedule definition
process.schedule = cms.Schedule(process.reflagging_step,process.rereco_step,process.out_step)

                                #process.reconstruction_step,process.dqmoffline_step,process.pathALCARECOSiStripCalZeroBias,process.pathALCARECOTkAlMinBias,process.pathALCARECOTkAlMuonIsolated,process.pathALCARECOMuAlCalIsolatedMu,process.pathALCARECOMuAlOverlaps,process.pathALCARECOHcalCalIsoTrk,process.pathALCARECOHcalCalDijets,process.endjob_step,process.out_step,process.ALCARECOStreamTkAlMinBiasOutPath,process.ALCARECOStreamTkAlMuonIsolatedOutPath,process.ALCARECOStreamMuAlOverlapsOutPath,process.ALCARECOStreamMuAlCalIsolatedMuOutPath,process.ALCARECOStreamHcalCalIsoTrkOutPath,process.ALCARECOStreamHcalCalDijetsOutPath,process.ALCARECOStreamSiStripCalZeroBiasOutPath)
