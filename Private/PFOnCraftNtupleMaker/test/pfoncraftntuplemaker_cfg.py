import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    "rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_v75.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged001.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged002.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged003.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged004.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged005.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged006.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged007.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged008.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged009.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged010.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged011.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged012.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged013.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged014.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged015.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged016.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged017.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged018.root",
#"rfio:/castor/cern.ch/user/f/fblekman/pftest_CMSSW326/pruned_pfreco_merged019.root"

   
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('MyCalibratableTree.root')
)
import RecoParticleFlow.PFAnalyses.pflowCalibratable_cfi as calibratable
process.demo = cms.EDAnalyzer('PFOnCraftNtupleMaker',
                              calibratable.EventDelegate

)


process.p = cms.Path(process.demo)
