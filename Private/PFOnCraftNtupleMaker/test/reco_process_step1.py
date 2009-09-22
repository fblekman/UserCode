import FWCore.ParameterSet.Config as cms

process = cms.Process("COSMICRECO")

# Messaging
process.load("FWCore.MessageService.MessageLogger_cfi")

# DQM services
#process.load("DQMServices.Core.DQM_cfg")

# DB Configuration
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load("CondCore.DBCommon.CondDBCommon_cfi")

# Geometry
process.load("Configuration.StandardSequences.Geometry_cff")
# Magnetic Field
process.load("Configuration.StandardSequences.MagneticField_38T_cff")#with MF
#process.load("Configuration.StandardSequences.MagneticField_cff")#0T
#process.load("Configuration.GlobalRuns.ForceZeroTeslaField_cff")#0T

# reconstruction sequence for Cosmics
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")


# Global Tag
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.connect = "frontier://FrontierProd/CMS_COND_21X_GLOBALTAG" # FOR 21X REPRO
process.GlobalTag.connect = "frontier://FrontierProd/CMS_COND_31X_GLOBALTAG" # FOR 31X REPRO

#process.GlobalTag.globaltag = "CRAFT_ALL_V4::All"
########tag behind for MC dataset
#process.GlobalTag.globaltag = 'COSMMC_21X_V1::All'
#process.GlobalTag.globaltag = 'COSMMC_22X_V1::All'
#process.GlobalTag.globaltag = "CRAFT_ALL_V11::All"
######## TAG FOR 31X REPROCESSING
#process.GlobalTag.globaltag = "GR09_31X_V6P::All"   #CRAFT09
process.GlobalTag.globaltag = "CRAFT09_R_V2::All"
#process.GlobalTag.globaltag = "CRAFT0831X_V1::All"   #CRAFT08


#############################
#  LATEST GAIN CALIBRATION
process.test = cms.ESSource("PoolDBESSource",
                                        DBParameters = cms.PSet(
                                           messageLevel = cms.untracked.int32(0),
					   #authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')
                                           authenticationPath = cms.untracked.string('0')
                                        ),
                                        connect = cms.string("sqlite_file:/afs/cern.ch/user/m/mucib/public/prova.db"),
                                        toGet = cms.VPSet(cms.PSet(record = cms.string("SiPixelGainCalibrationOfflineRcd"),
                                                                   tag = cms.string("GainCalib_TEST_offline"))
                                                          )
                                        )
process.es_prefer_test = cms.ESPrefer("PoolDBESSource","test")



## Load and Configure OfflineValidation
process.load("Alignment.OfflineValidation.TrackerOfflineValidation_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
   #comment this line when ALCARECO
   #inputCommands = cms.untracked.vstring('keep *',"drop *_*_*_FU"),
   #lastRun = cms.untracked.uint32(109624),
   #timetype = cms.string('runnumber'),
   #firstRun = cms.untracked.uint32(109011),
   #interval = cms.uint32(1),

#replace 'myfile.root' with the source file you want to use
     fileNames = cms.untracked.vstring(
    "infile.root"
    )
)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('myfile_localreco.root')
                               )

process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.TrackRefitterP5.src = 'TRACKINPUTTAG'
process.TrackRefitterP5.TrajectoryInEvent = True

process.load("RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilderWithoutRefit_cfi")



process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = 'INFO'#32X version
#process.MessageLogger.cerr.threshold = 'Info' #22X version

#RECONSTRUCTION, Hybrid, local reco comes from Configurations/StandardSequences/reconstructionCosmics_cff,
process.load("Configuration.StandardSequences.ReconstructionCosmics_cff")
#process.load("Configuration.StandardSequences.Reconstruction_cff")
# fix all missing labels


process.p = cms.Path(process.RawToDigi*process.reconstructionCosmics)#*process.globalreco*process.highlevelreco)

    
process.TrackerDigiGeometryESModule.applyAlignment = True

process.ep = cms.EndPath( process.out )



