# Auto generated configuration file
# using: 
# Revision: 1.156 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: FREYACONFIGDATA -s RECO --processName=METSIGNIF --conditions=FrontierConditions_GlobalTag,IDEAL_34X:All --data --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('METSIGNIF')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('FREYACONFIGDATA nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('OtherCMS', 
                                    'StdException', 
                                    'Unknown', 
                                    'BadAlloc', 
                                    'BadExceptionType', 
                                    'ProductNotFound', 
                                    'DictionaryNotFound', 
                                    'InsertFailure', 
                                    'Configuration', 
                                    'LogicError', 
                                    'UnimplementedFeature', 
                                    'InvalidReference', 
                                    'NullPointerError', 
                                    'NoProductSpecified', 
                                    'EventTimeout', 
                                    'EventCorruption', 
                                    'ScheduleExecutionFailure', 
                                    'EventProcessorFailure', 
                                    'FileInPathError', 
                                    'FileOpenError', 
                                    'FileReadError', 
                                    'FatalRootError', 
                                    'MismatchedInputFiles', 
                                    'ProductDoesNotSupportViews', 
                                    'ProductDoesNotSupportPtr', 
                                    'NotFound')
    )
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Jan23Skim-v1/0015/FA138089-F809-DF11-ADAD-002618943945.root'
    )
)

process.metsignifntup1 = cms.EDAnalyzer('METSignifFirstData'
)
process.metsignifntup2 = cms.EDAnalyzer('METSignifFirstData'
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("histo.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )
				

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = cms.untracked.vstring( "drop *","keep *_*_*_METSIGNIF","keep *_gtDigis_*_*","keep *_particleFlow_*_*","keep *_offlinePrimaryVertices_*_*","keep *_generalTracks_*_*"),

#    process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('FREYACONFIGMC_RECO.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'STARTUP3X_V13:All'
#process.GlobalTag.globaltag = 'GR09_R_34X_V3::All'


# Path and EndPath definitions
process.reconstruction_step = cms.Path(process.reconstruction)
process.metreco_step = cms.Path(process.metsignifntup1+process.ecalDigis+process.ecalPreshowerDigis+process.hcalDigis+process.calolocalreco+process.caloTowersRec+process.metreco+process.recoPFMET+process.metsignifntup2)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.metreco_step,process.endjob_step,process.out_step)
