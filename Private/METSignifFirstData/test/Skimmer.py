# Auto generated configuration file
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
    version = cms.untracked.string('$Revision: 1.1 $'),
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
                '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0001/DC74A007-122D-DF11-88BC-0017A4770414.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0001/5EF46CD4-6C2D-DF11-BD31-001A4BA98052.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/F6D7AC19-052B-DF11-98C4-0017A4771020.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/F00638E4-122B-DF11-A54C-0017A4770C10.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/EC8FED12-0E2B-DF11-A585-0017A4770434.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/EC865CF0-1A2B-DF11-91F4-001CC47A52B6.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/E8E141AF-0B2B-DF11-994A-0017A4770000.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/E67D16B8-162B-DF11-9396-00237DA12CA0.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/E482B811-0E2B-DF11-B28E-00237DA15C00.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/E0F5805E-112B-DF11-B954-00237DA1494E.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/E07631D7-A22B-DF11-B7DB-0017A477040C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/DE5D7D61-0B2B-DF11-8523-0017A4770420.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/DCEE6C1B-052B-DF11-ACD7-001E0B5FB53E.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/DC8888DF-122B-DF11-B4B0-0017A4770C04.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/DAF23CF6-022B-DF11-805C-00237DA14F92.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/DAA107B9-162B-DF11-BD3D-00237DA14F86.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/D6909EB8-082B-DF11-A92B-0017A477102C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/D4AF3FAD-052B-DF11-980D-0017A4770C38.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/CEC7EB2F-102B-DF11-B8E9-00237DA10D06.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/CC56C9E1-0C2B-DF11-A22F-0017A4771038.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/CAF67E00-812B-DF11-BB06-0017A4770028.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/C88C706D-162B-DF11-A246-001E0B472C96.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/C88851E0-0C2B-DF11-8AEB-00237DA13FB6.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/C63D37F6-172B-DF11-AD6F-0017A4770000.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/C49170C9-042B-DF11-A49D-0017A4771014.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/C064A769-FF2A-DF11-A1C2-0017A4771010.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/B8148BA5-142B-DF11-BE29-0017A4770C28.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/B69963C1-012B-DF11-A6FA-0022649F01AA.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/AE1F7AB8-162B-DF11-8BFE-0017A4770C08.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/ACC34D6C-1F2B-DF11-A924-0017A4771034.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/A6EBA5BD-082B-DF11-B1B7-0017A477102C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/A2EE2359-0E2B-DF11-BCD2-0017A477042C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/A0C7A872-162B-DF11-B78B-0017A477101C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/9E7971C1-082B-DF11-9451-0017A4770818.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/9E229FC2-082B-DF11-BA33-0017A4771020.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/96D370EF-1A2B-DF11-87A2-0017A4770C10.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/92ED0DB9-162B-DF11-8D0D-0017A477002C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/9216FD34-0A2B-DF11-9C14-00237D9F2120.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/90426AB0-0B2B-DF11-97C1-0017A4770014.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/8A98F3B7-162B-DF11-92E6-0017A4770018.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/84A01AC9-042B-DF11-AA3E-00237DA1CD7E.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/808DA530-102B-DF11-99E0-0017A4770028.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/7E8B31AD-0B2B-DF11-B5AB-00237DA15C5E.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/7E6733B4-082B-DF11-BD9D-0017A4770418.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/7C5889F6-172B-DF11-AFD2-0017A4770418.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/768316B9-162B-DF11-ADA5-002264055CE4.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/6C20D812-0E2B-DF11-986B-001E0B482944.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/6C1D72F7-052B-DF11-9155-0017A4771020.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/685CE7DF-032B-DF11-8AB3-0017A4770C30.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/64CED4DF-0C2B-DF11-AC38-0017A4770828.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/5EF1A5F3-052B-DF11-889D-001E0B5FA528.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/5C218968-022B-DF11-B5C6-0017A4770418.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/5AD9EB11-0E2B-DF11-B86A-0017A477003C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/5AB854AF-052B-DF11-A4B6-0017A4770434.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/5A051558-0E2B-DF11-9A77-0017A4770C10.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/4E75BE6F-162B-DF11-B191-001CC416B30C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/4AE212AF-052B-DF11-9846-001E0B48D108.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/4A1B56D5-082B-DF11-AAB2-001E0B5F5898.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/48A97A19-052B-DF11-9467-0017A4771010.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/46885E0E-0E2B-DF11-9AB0-0017A4770C08.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/42C69A79-042B-DF11-AB05-001CC4164B40.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/42AF0CC9-042B-DF11-A941-0017A477080C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/404F38C4-082B-DF11-8B72-0017A477103C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/3EB7F0B9-082B-DF11-BA05-0017A4770434.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/3A9CDE2E-102B-DF11-A8B4-001E0B5F95B2.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/3282C100-232B-DF11-83F4-00215AAC88D2.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/3277A1BE-082B-DF11-A4B9-00237DA41368.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/2E41FE62-142B-DF11-A0CC-0017A4770404.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/2CDFEE65-022B-DF11-B0C7-0017A4771028.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/2CBEC434-102B-DF11-8D49-00237DA15C66.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/2CB19E14-0E2B-DF11-9985-001E0B482944.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/2CA335B4-082B-DF11-8089-00237DA41368.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/2C3CD6AD-022B-DF11-ADEC-0017A4770414.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/2A1B6A6F-162B-DF11-99E0-00237DA1ED1C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/2A1B18B6-082B-DF11-A7AD-0017A4770034.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/26B4F45E-112B-DF11-BB98-00237DA16C5E.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/26440611-0E2B-DF11-81E9-0017A4770404.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/2446BD6F-162B-DF11-9BFF-0017A4770410.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/241807C3-012B-DF11-B88E-0017A477041C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/2206BC35-1B2B-DF11-88FA-00237DA1ED1C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/206F0BE8-0C2B-DF11-8DB7-0017A477102C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/1EA1AC19-052B-DF11-BCF4-00237DA1ED4A.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/1E142BAA-142B-DF11-9D37-00237DA14F86.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/1A375DF6-172B-DF11-BAF6-0017A4770C10.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/180C07AE-0B2B-DF11-841D-00237DA14F92.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/1803DCB7-082B-DF11-999B-0017A4770C04.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/1682D500-172B-DF11-9B9F-0017A4770C1C.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/12D5B36A-022B-DF11-9DDF-0017A4771024.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/129A56F7-172B-DF11-A506-0017A4770014.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/10854CE1-122B-DF11-B18B-0017A4771000.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/0E6B1933-102B-DF11-B1B9-0017A4770000.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/08EF28F3-172B-DF11-B8B3-001B78E1096E.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/08CEC719-052B-DF11-9FA7-0017A4770400.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/06EBF0B6-052B-DF11-8C0A-0017A4771028.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/06384AB9-162B-DF11-BB86-0017A4771028.root',
				        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Mar3rdSkim_v2/0000/0232E895-002B-DF11-BEF4-001CC416B30C.root'
				
    )
)

			
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(
   '123596:2-123596:9999',    
   '123615:70-123615:9999',   
   '123732:62-123732:9999',  
   '123815:8-123815:9999',   
   '123818:2-123818:42',     
   '123908:2-123908:12',     
   '124008:1-124008:1',      
   '124009:1-124009:68',     
   '124020:12-124020:94',     
   '124022:66-124022:179',    
   '124023:38-124023:9999',   
   '124024:2-124024:83',     
   '124025:5-124025:13',    
   '124027:24-124027:9999',   
   '124030:2-124030:9999',  
   '124120:1-124120:9999',   
   '124275:3-124275:30'
   )

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*", "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT")

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
																	   "keep *_generalV0Candidates_*_*"),
    fileName = cms.untracked.string('/tmp/fblekman/SKIMMED_DATA_withNewMETandGTdigis.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RECO'),
        filterName = cms.untracked.string('')
    )
)



# Other statements
process.GlobalTag.globaltag = 'GR09_R_34X_V4::All'


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

process.metreco_step = cms.Path(process.caloTowersRec+process.metreco+process.pfMet)

# Schedule definition
process.schedule = cms.Schedule(process.metreco_step,process.endjob_step,process.out_step)
                                #process.reconstruction_step,process.dqmoffline_step,process.pathALCARECOSiStripCalZeroBias,process.pathALCARECOTkAlMinBias,process.pathALCARECOTkAlMuonIsolated,process.pathALCARECOMuAlCalIsolatedMu,process.pathALCARECOMuAlOverlaps,process.pathALCARECOHcalCalIsoTrk,process.pathALCARECOHcalCalDijets,process.endjob_step,process.out_step,process.ALCARECOStreamTkAlMinBiasOutPath,process.ALCARECOStreamTkAlMuonIsolatedOutPath,process.ALCARECOStreamMuAlOverlapsOutPath,process.ALCARECOStreamMuAlCalIsolatedMuOutPath,process.ALCARECOStreamHcalCalIsoTrkOutPath,process.ALCARECOStreamHcalCalDijetsOutPath,process.ALCARECOStreamSiStripCalZeroBiasOutPath)
