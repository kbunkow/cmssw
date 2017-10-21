import FWCore.ParameterSet.Config as cms

omtfParamsSource = cms.ESSource(
    "EmptyESSource",
    recordName = cms.string('L1TMuonOverlapParamsRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

###OMTF ESProducer. Fills CondFormats from XML files.
omtfParams = cms.ESProducer(
    "L1TMuonOverlapParamsESProducer",
    patternsXMLFiles = cms.VPSet(
        #cms.PSet(patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuonOverlap/test/expert/Patterns_0x00020007_scaledEndcaps2.xml")),
        cms.PSet(patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_0x00020007.xml")),#this are the patterns used in 2017 in p5
        #cms.PSet(patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_0x0003.xml")),
        #cms.PSet(patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuonOverlap/test/GPs08Norm1_78.xml")),        
    ),
    configXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/hwToLogicLayer_0x0004.xml"),
    #configXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/hwToLogicLayer_0x0006.xml"),
)




