# Auto generated configuration file
# using: 
# Revision: 1.155 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: reco -s RAW2DIGI,RECO --conditions FrontierConditions_GlobalTag,STARTUP_V4::All --eventcontent RECOSIM
import FWCore.ParameterSet.Config as cms

process = cms.Process('DIGIS')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/RawToDigi_Data_cff')
process.load('Configuration/StandardSequences/L1Reco_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('DQMOffline/Configuration/DQMOffline_cff')
process.load('Configuration/StandardSequences/AlCaRecoStreams_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.155 $'),
    annotation = cms.untracked.string('reco nevts:1'),
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
    fileNames = cms.untracked.vstring('/store/data/BeamCommissioning09/ZeroBias/RAW/v1/000/123/151/FC28B891-0EDE-DE11-B6D6-001D09F25217.root',
    '/store/data/BeamCommissioning09/ZeroBias/RAW/v1/000/123/151/F2E0B5AB-F3DD-DE11-BC66-001D09F2932B.root',
    '/store/data/BeamCommissioning09/ZeroBias/RAW/v1/000/123/151/ECDDE02B-FCDD-DE11-B06B-001D09F2AF1E.root',
    '/store/data/BeamCommissioning09/ZeroBias/RAW/v1/000/123/151/B6975215-FADD-DE11-8031-001D09F28EC1.root',
    '/store/data/BeamCommissioning09/ZeroBias/RAW/v1/000/123/151/A6300CE4-FCDD-DE11-8AF8-001D09F25456.root',
    '/store/data/BeamCommissioning09/ZeroBias/RAW/v1/000/123/151/72B97676-FBDD-DE11-A454-001D09F23944.root',
    '/store/data/BeamCommissioning09/ZeroBias/RAW/v1/000/123/151/509653C1-FADD-DE11-A408-001D09F2438A.root',
    '/store/data/BeamCommissioning09/ZeroBias/RAW/v1/000/123/151/40F845AD-0BDE-DE11-AE51-001D09F2910A.root',
    '/store/data/BeamCommissioning09/ZeroBias/RAW/v1/000/123/151/22E2AC6C-05DE-DE11-86C7-001D09F282F5.root',
    '/store/data/BeamCommissioning09/ZeroBias/RAW/v1/000/123/151/0EB6E8E2-FCDD-DE11-9191-001D09F24DDF.root',
    '/store/data/BeamCommissioning09/ZeroBias/RAW/v1/000/123/151/0A39B1A4-F3DD-DE11-9C83-000423D94E70.root'
                                      )
)

process.pixelDigiFilter = cms.EDFilter("PixelEventFilter",
                                       saveEventTail = cms.untracked.uint32(10),
                                       minDigisToReset = cms.untracked.uint32(50)
                                       )
# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
#    outputCommands = process.RECOSIMEventContent.outputCommands,
                                  outputCommands = cms.untracked.vstring("keep *","drop *_source_*_*"),
    fileName = cms.untracked.string('reco_RAW2DIGI_RECO.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('pixelDigiFilter')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'FIRSTCOLL::All'

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.siPixelDigis*process.pixelDigiFilter)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.endjob_step,process.out_step)
