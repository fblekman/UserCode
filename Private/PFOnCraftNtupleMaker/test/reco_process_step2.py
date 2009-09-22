# Auto generated configuration file
# using: 
# Revision: 1.138 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: -s RECO --filein=/tmp/fblekman/myfile_z.root --no_exec --conditions=FrontierConditions_GlobalTag,CRAFT09_R_V2::All
import FWCore.ParameterSet.Config as cms

process = cms.Process('PFRECO')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.138 $'),
    annotation = cms.untracked.string('-s nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:myfile_localreco.root')
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = cms.untracked.vstring("drop *","keep *_*_*_COSMICRECO","keep *_particleFlow_*_PFRECO","keep recoTracks_*_*_PFRECO","keep recoTrackExtras_*_*_PFRECO","keep recoPF*_*_*_PFRECO"),
    fileName = cms.untracked.string('myfile_pfreco.root'),
                                  
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'CRAFT09_R_V2::All'

# Path and EndPath definitions
process.reconstruction_step= cms.Path(process.reconstruction)

process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step,process.out_step)
