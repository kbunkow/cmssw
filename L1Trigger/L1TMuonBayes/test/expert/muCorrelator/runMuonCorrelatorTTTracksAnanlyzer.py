# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from sqlalchemy.sql.expression import false

import os
import sys
import commands

makeTTracks = True 

processName = "L1TMuonEmulation"
if makeTTracks :
    processName = "L1TMuonEmulation"
else :
    processName = "L1TMuonAnalysis"
    
process = cms.Process(processName)

process.load("FWCore.MessageLogger.MessageLogger_cfi")

verbose = True



if verbose: 
    process.MessageLogger = cms.Service("MessageLogger",
       #suppressInfo       = cms.untracked.vstring('AfterSource', 'PostModule'),
       destinations   = cms.untracked.vstring(
                                               #'detailedInfo',
                                               #'critical',
                                               #'cout',
                                               #'cerr',
                                               'muCorrelatorEventPrint'
                    ),
       categories        = cms.untracked.vstring('l1tMuBayesEventPrint'),
       muCorrelatorEventPrint = cms.untracked.PSet(    
                         extension = cms.untracked.string('.txt'),                
                         threshold = cms.untracked.string('DEBUG'),
                         default = cms.untracked.PSet( limit = cms.untracked.int32(0) ), 
                         #INFO   =  cms.untracked.int32(0),
                         #DEBUG   = cms.untracked.int32(0),
                         l1tMuBayesEventPrint = cms.untracked.PSet( limit = cms.untracked.int32(100000000) )
                       ),
       debugModules = cms.untracked.vstring('L1TMuonBayesMuCorrelatorTrackProducer', 'OmtfTTAnalyzer', 'simOmtfDigis', 'omtfTTAnalyzer', 'simBayesMuCorrelatorTrackProducer') 
       #debugModules = cms.untracked.vstring('*')
    )

    #process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
if not verbose:
    process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
    process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False), 
                                         #SkipEvent = cms.untracked.vstring('ProductNotFound') 
                                     )


#######################################TTTracks################################################
GEOMETRY = "D17"

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D41_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

Source_Files = cms.untracked.vstring(
#        "/store/relval/CMSSW_10_0_0_pre1/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/94X_upgrade2023_realistic_v2_2023D17noPU-v2/10000/06C888F3-CFCE-E711-8928-0CC47A4D764C.root"
         #"/store/relval/CMSSW_9_3_2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/0681719F-AFA6-E711-87C9-0CC47A4C8E14.root"
         #"file:///eos/user/k/kbunkow/cms_data/0681719F-AFA6-E711-87C9-0CC47A4C8E14.root"
         #'file:///eos/user/k/kbunkow/cms_data/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100_L1TnoPU_F4EEAE55-C937-E811-8C29-48FD8EE739D1_dump1000Events.root'
         #'file:///eos/user/k/kbunkow/cms_data/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100_L1TnoPU_F4EEAE55-C937-E811-8C29-48FD8EE739D1_dump1000Events.root'
         #"file:///eos/cms/store/group/upgrade/sandhya/SMP-PhaseIIFall17D-00001.root"
         #'file:///afs/cern.ch/work/k/kbunkow/private/omtf_data/SingleMu_15_p_1_1_qtl.root' 
         #'file:///eos/user/k/kbunkow/cms_data/mc/PhaseIIFall17D/HSCPppstau_M_871_noPU_18156A80-66EC-E811-AE02-0CC47AFCC62A_dump2000Events.root'
         #'file:///eos/user/k/kbunkow/cms_data/mc/PhaseIIFall17D/HSCPppstau_M_200_noPU_FE3D8AD6-B6D0-E811-8FBD-141877412793_dump2000Events.root'
         #"/store/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/00000/00B70D7F-333F-E811-9095-0CC47A4D9A10.root",
         #'/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU300_93X_upgrade2023_realistic_v5-v1/110000/FEFE2FE6-32E8-E811-9995-0242AC130002.root'
         #'/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/FEAC79F7-5B5C-E811-B830-0025905C2D9A.root'
         #'file:///eos/user/a/akalinow/Data/SingleMu/9_3_14_FullEta_v2/SingleMu_26_m_1.root'
         #'file:///eos/user/a/akalinow/Data/9_3_14_HSCP_v5/HSCPppstau_M_494_TuneZ2star_13TeV_pythia6_cff_py_GEN_SIM_DIGI_L1_L1TrackTrigger_DIGI2RAW_HLT_1.root'
         #'file:///afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_5_0_pre1/src/L1Trigger/L1TMuonBayes/test/SingleMu_PU200_32DF01CC-A342-E811-9FE7-48D539F3863E_dump500Events.root'
         #'file:///eos/user/k/kbunkow/cms_data/mc/PhaseIIFall17D/SingleMu_PU200_32DF01CC-A342-E811-9FE7-48D539F3863E_dump500Events.root'
         #'file:///eos/user/k/kbunkow/cms_data/mc/PhaseIIFall17D/ZMM_EE29AF8E-51AF-E811-A2BD-484D7E8DF0D3_dump1000Events.root'
         #'file:///eos/user/a/akalinow/Data/SingleMu/SingleMuFlatPt_50GeVto10GeV_cfi_py_GEN_SIM_DIGI_L1_L1TrackTrigger_DIGI2RAW_HLT.root'
         #'file:///eos/user/k/kbunkow/cms_data/mc/PhaseIITDRSpring19DR/GluGluHToZZTo4L_NoPU_FB98030B-16C5-9842-9698-8371EB8D8B01_dump1000Ev.root'
         #'file:///eos/user/k/kbunkow/cms_data/mc/PhaseIITDRSpring19DR/HSCPppstau_M_200_NoPU_F9EFA962-CDD8-C643-8F62-8F75384875F0_dump4000Ev.root'
         #'file:///eos/user/k/kbunkow/cms_data/mc/PhaseIITDRSpring19DR/JPsiToMuMu_Pt0to100_NoPU_FDA71CB6-4C3B-4540-99EB-803077C6EC2D_dump4000Ev.root'
        #'file:///eos/user/k/kbunkow/cms_data/mc/PhaseIITDRSpring19DR/PhaseIITDRSpring19DR_HSCPppstau_M_871__noPU_v32_F9357CE3-E1BD-C64C-8F43-895CFA3A0AFC_dump1000Ev.root'
        #'file:///eos/user/k/kbunkow/cms_data/mc/PhaseIITDRSpring19DR/PhaseIITDRSpring19DR_HSCPppstau_M_200__noPU_v32_A91AA4D8-5187-5544-8304-365404899406_dump1000Ev.root'
        #"file:///eos/user/k/kbunkow/cms_data/mc/PhaseIITDRSpring19DR/PhaseIITDRSpring19DR_Mu_FlatPt2to100_noPU_v31_E0D5C6A5-B855-D14F-9124-0B2C9B28D0EA_dump4000Ev.root"
        "file:///eos/user/k/kbunkow/cms_data/mc/PhaseIITDRSpring19DR/HSCPppstau_M_200_PU200_v3-v1_ACF9C8E2-0570-6A4A-983A-E2B230F6FCAA_dump300Ev.root"
        # "file:///eos/user/k/kbunkow/cms_data/mc/PhaseIITDRSpring19DR/HSCPppstau_M_871_PU200_v3-v2_1ADE9D9E-8C0C-1948-A405-5DFDA1AF5172_dump100Ev.root"
        #"file:///afs/cern.ch/work/k/kbunkow/public/CMSSW/cmssw_10_x_x_l1tOfflinePhase2/CMSSW_10_6_1_patch2/src/L1Trigger/L1TMuonBayes/test/expert/muCorrelator/outCollections_HSCPppstau_M_200_PU200_v2.root"
)


process.source = cms.Source("PoolSource", fileNames = Source_Files,
        inputCommands=cms.untracked.vstring(
        'keep *',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
        'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016s_simEmtfDigis__HLT',
          'drop l1tHGCalClusterBXVector_hgcalTriggerPrimitiveDigiProducer_cluster2D_HLT',
          'drop *HGCal*_*_*_*',
          'drop *hgcal*_*_*_*',
          'drop *Ecal*_*_*_*',
          'drop *Hcal*_*_*_*',
          'drop *Calo*_*_*_*',
          
          'drop *_*HGCal*_*_*',
          'drop *_*hgcal*_*_*',
          'drop *_*Ecal*_*_*',
          'drop *_*Hcal*_*_*',
          'drop *_*Calo*_*_*',
          
          'drop *_*_*HGCal*_*',
          'drop *_*_*hgcal*_*',
          'drop *_*_*Ecal*_*',
          'drop *_*_*Hcal*_*',
          'drop *_*_*Calo*_*')
)

process.TFileService = cms.Service("TFileService", fileName = cms.string('muCorrelatorTTAnalysis1_HSCPppstau_M_200_PU200.root'), closeFileFast = cms.untracked.bool(True))
#process.TFileService = cms.Service("TFileService", fileName = cms.string('muCorrelatorTTAnalysis1_JPsiToMuMu_Pt0to100_NoPU.root'), closeFileFast = cms.untracked.bool(True))
#process.TFileService = cms.Service("TFileService", fileName = cms.string('muCorrelatorTTAnalysis1_HSCPppstau_M_200_NoPU.root'), closeFileFast = cms.untracked.bool(True))
#process.TFileService = cms.Service("TFileService", fileName = cms.string('muCorrelatorTTAnalysis1_HSCPppstau_M_200_NoPU.root'), closeFileFast = cms.untracked.bool(True))

#from L1Trigger.Configuration.customiseUtils import L1TrackTriggerTracklet, L1TTurnOffHGCalTPs_v9, configureCSCLCTAsRun2
#process = L1TrackTriggerTracklet(process)
#process = L1TTurnOffHGCalTPs_v9(process)
#process = configureCSCLCTAsRun2(process) #has no efect when csc digis alreaady taken form the data file


if makeTTracks :
    ############################################################
    # remake L1 stubs and/or cluster/stub truth ??
    ############################################################
    
    process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
    from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
    process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff")
    
    #if GEOMETRY == "D10": 
    #    TTStubAlgorithm_official_Phase2TrackerDigi_.zMatchingPS = cms.bool(False)
    
    if GEOMETRY != "TkOnly": 
        from SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff import *
        TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")
    
    process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)
    process.TTClusterStubTruth = cms.Path(process.TrackTriggerAssociatorClustersStubs)
    
    
    ############################################################
    # L1 tracking
    ############################################################
    
    #from L1Trigger.TrackFindingTracklet.Tracklet_cfi import *
    #if GEOMETRY == "D10": 
    #    TTTracksFromTracklet.trackerGeometry = cms.untracked.string("flat")
    #TTTracksFromTracklet.asciiFileName = cms.untracked.string("evlist.txt")
    
    process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")
    process.TTTracks = cms.Path(process.L1TrackletTracks)
    process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators)


#######################################TTTracks################################################

##This overrides the tracker geometry and the TTTriger does not work!!!!!!!!!!!!
# # PostLS1 geometry used
# process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
# process.load('Configuration.Geometry.GeometryExtended2015_cff')
# ############################
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
# from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


####Event Setup Producer
# process.load('L1Trigger.L1TMuonOverlap.fakeOmtfParams_cff')
# process.esProd = cms.EDAnalyzer("EventSetupRecordDataGetter",
#    toGet = cms.VPSet(
#       cms.PSet(record = cms.string('L1TMuonOverlapParamsRcd'),
#                data = cms.vstriL1TMuonBayesrlapParams'))
#                    ),
#    verbose = cms.untracked.bool(False)
# )


####OMTF Emulator
process.load('L1Trigger.L1TMuonBayes.simBayesMuCorrelatorTrackProducer_cfi')

#process.TFileService = cms.Service("TFileService", fileName = cms.string('muCorrelatorHists.root'), closeFileFast = cms.untracked.bool(True))

process.simBayesMuCorrelatorTrackProducer.ttTracksSource = cms.string("L1_TRACKER")
#process.simBayesMuCorrelatorTrackProducer.ttTracksSource = cms.string("SIM_TRACKS") #TODO !!!!!!!
process.simBayesMuCorrelatorTrackProducer.pdfModuleType = cms.string("PdfModuleWithStats") #TODO
#process.simBayesMuCorrelatorTrackProducer.pdfModuleFile = cms.FileInPath("L1Trigger/L1TMuonBayes/test/pdfModule.xml") #TODO!!!!!!!!!!!!!!!!!!!!!!!!!!11
#process.simBayesMuCorrelatorTrackProducer.pdfModuleFile = cms.FileInPath("L1Trigger/L1TMuonBayes/test/pdfModuleSimTracks100FilesWithiRPC.xml")
#process.simBayesMuCorrelatorTrackProducer.timingModuleFile  = cms.FileInPath("L1Trigger/L1TMuonBayes/test/muTimingModule100FilesWithiRPC.xml")
#process.simBayesMuCorrelatorTrackProducer.timingModuleFile  = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/muTimingModuleTest.xml")
#process.simBayesMuCorrelatorTrackProducer.pdfModuleFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/pdfModuleSimTracks100FilesSigma1p3.xml")  

process.simBayesMuCorrelatorTrackProducer.generateTiming = cms.bool(False)
process.simBayesMuCorrelatorTrackProducer.useStubsFromAdditionalBxs = cms.int32(3)
#process.simBayesMuCorrelatorTrackProducer.bxRangeMin = cms.int32(-3)
#process.simBayesMuCorrelatorTrackProducer.bxRangeMax = cms.int32(3)

process.dumpED = cms.EDAnalyzer("EventContentAnalyzer")
process.dumpES = cms.EDAnalyzer("PrintEventSetupContent")

process.L1TMuonSeq = cms.Sequence( #process.esProd +         
                                   process.simBayesMuCorrelatorTrackProducer 
                                   #+ process.dumpED
                                   #+ process.dumpES
)

process.L1TMuonPath = cms.Path(process.L1TMuonSeq)

if makeTTracks:
    process.out = cms.OutputModule("PoolOutputModule", 
        fileName = cms.untracked.string("outCollections.root"),
        outputCommands=cms.untracked.vstring(
            'drop *', 
            #'keep l1BayesMuCorrelatorTrackBXVector_simBayesMuCorrelatorTrackProducer_AllTracks_L1TMuonEmulation'
            'keep *_*_*_L1TMuonEmulation',
            'keep TrackingParticles_mix_MergedTrackTruth_HLT'
            #'keep l1tRegionalMuonCandBXVector_simOmtfDigis_OMTF_HLT'
            )
    )
    process.output_step = cms.EndPath(process.out)

############################################################

analysisType = "efficiency" # or rate
  
for a in sys.argv :
    if a == "efficiency" or a ==  "rate" or a == "withTrackPart" :
        analysisType = a
        break;
    
print "analysisType=" + analysisType

process.omtfTTAnalyzer= cms.EDAnalyzer("MuCorrelatorAnalyzer", 
                                 outRootFile = cms.string("muCorrelatorTTAnalysis1.root"),
                                 etaCutFrom = cms.double(0.), #OMTF eta range
                                 etaCutTo = cms.double(2.4),
                                          
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(verbose),      # printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   # save *all* L1 tracks, not just truth matched to primary particle
                                       SaveStubs = cms.bool(False),      # save some info for *all* stubs
                                       LooseMatch = cms.bool(True),     # turn on to use "loose" MC truth association
                                       L1Tk_nPar = cms.int32(4),         # use 4 or 5-parameter L1 track fit ??
                                       L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
                                       TP_minNStub = cms.int32(4),       # require TP to have >= X number of stubs associated with it
                                       TP_minNStubLayer = cms.int32(4),  # require TP to have stubs in >= X layers/disks
                                       TP_minPt = cms.double(2.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.4),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       TP_maxRho = cms.double(30.0),     # for efficiency analysis, to not inlude the muons from the far decays 
                                       L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),               ## TTTrack input
                                       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), ## MCTruth input 
                                       # other input collections
                                       L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
                                       MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
                                       MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
                                       TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       
                                       muCandQualityCut = cms.int32(12),
                                       analysisType = cms.string(analysisType)
                                        )
process.omtfTTAnalyzerPath = cms.Path(process.omtfTTAnalyzer)

###########################################################

############################################################

# use this if you want to re-run the stub making
#process.schedule = cms.Schedule(process.TTClusterStub,process.TTClusterStubTruth,process.TTTracksWithTruth, process.L1TMuonPath, process.omtfTTAnalyzerPath)

# use this if cluster/stub associators not available
# process.TTClusterStub, process.TTTracks, 

if makeTTracks:
    process.schedule = cms.Schedule(process.TTTracksWithTruth, process.L1TMuonPath, process.omtfTTAnalyzerPath) #TODO default
else :
    process.schedule = cms.Schedule(process.omtfTTAnalyzerPath)  
    
#process.schedule = cms.Schedule(process.TTTracks, process.TTTracksWithTruth, process.L1TMuonPath, process.omtfTTAnalyzerPath)
#process.schedule = cms.Schedule(process.TTClusterStub, process.TTClusterStubTruth, process.TTTracksWithTruth, process.TTTracks, process.L1TMuonPath, process.omtfTTAnalyzerPath)

# use this to only run tracking + track associator
#process.schedule = cms.Schedule(process.TTTracksWithTruth,process.ana)

#if makeTTracks: TODO
#    process.schedule.extend([process.output_step])
