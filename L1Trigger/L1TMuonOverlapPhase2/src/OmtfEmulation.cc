/*
 * OmtfEmulation.cpp
 *
 *  Created on: May 20, 2020
 *      Author: kbunkow
 */

#include <memory>

#include "L1Trigger/L1TMuonOverlapPhase2/interface/OmtfEmulation.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/InputMakerPhase2.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/PtAssignmentNNRegression.h"

#include "DataFormats/L1TMuonPhase2/interface/Constants.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <memory>

OmtfEmulation::OmtfEmulation(const edm::ParameterSet& edmParameterSet,
                             MuStubsInputTokens& muStubsInputTokens,
                             MuStubsPhase2InputTokens& muStubsPhase2InputTokens)
    : OMTFReconstruction(edmParameterSet, muStubsInputTokens), muStubsPhase2InputTokens(muStubsPhase2InputTokens) {}

void OmtfEmulation::beginJob() {
  if (edmParameterSet.exists("usePhase2DTPrimitives") && edmParameterSet.getParameter<bool>("usePhase2DTPrimitives")) {
    inputMaker = std::make_unique<InputMakerPhase2>(edmParameterSet,
                                                    muStubsInputTokens,
                                                    muStubsPhase2InputTokens,
                                                    omtfConfig.get(),
                                                    std::make_unique<OmtfPhase2AngleConverter>());
  } else {
    inputMaker = std::make_unique<OMTFinputMaker>(
        edmParameterSet, muStubsInputTokens, omtfConfig.get(), std::make_unique<OmtfAngleConverter>());
  }

  //.....................rrrrrrrrccccdddddd
  //.....................765432109876543210
  firedLayersToQuality[0b000000110000000011] = 1;
  firedLayersToQuality[0b000000100000000011] = 1;
  firedLayersToQuality[0b000000010000000011] = 1;
  firedLayersToQuality[0b000000110000000001] = 1;
  firedLayersToQuality[0b000001000000001100] = 1;
  firedLayersToQuality[0b000011000000001100] = 1;
  firedLayersToQuality[0b000010000000001100] = 1;
  firedLayersToQuality[0b000011000000000100] = 1;
  firedLayersToQuality[0b000000011000000001] = 1;
  firedLayersToQuality[0b001000010000000001] = 1;

  firedLayersToQuality[0b000100000000110000] = 1;
  firedLayersToQuality[0b001100000000010000] = 1;

  firedLayersToQuality[0b010000110000000001] = 8;
  firedLayersToQuality[0b000000111110000001] = 8;
  firedLayersToQuality[0b000000001000000011] = 8;
  firedLayersToQuality[0b000000111000000001] = 8;
  firedLayersToQuality[0b000000101000000001] = 8;
  firedLayersToQuality[0b010000011000000001] = 8;
  firedLayersToQuality[0b010000100000000001] = 8;
  firedLayersToQuality[0b000000110100000001] = 8;
  firedLayersToQuality[0b000000100100000001] = 8;
  firedLayersToQuality[0b001000100000000001] = 8;
  firedLayersToQuality[0b010000010000000001] = 8;
  firedLayersToQuality[0b001000110000000001] = 8;
  firedLayersToQuality[0b001000110000000000] = 8;
  firedLayersToQuality[0b000000010100000001] = 8;
  firedLayersToQuality[0b000010100000000001] = 8;
  firedLayersToQuality[0b000000100010000001] = 8;
  firedLayersToQuality[0b001010010000000101] = 8;
  firedLayersToQuality[0b100000000000000011] = 8;
  firedLayersToQuality[0b011011000000000000] = 8;
  firedLayersToQuality[0b000010110000000001] = 8;
  firedLayersToQuality[0b001001110000000001] = 8;
  firedLayersToQuality[0b000010100000000101] = 8;
  firedLayersToQuality[0b000011110000000001] = 8;
  firedLayersToQuality[0b000011110000001101] = 8;
  firedLayersToQuality[0b000011100000000101] = 8;
  firedLayersToQuality[0b000011110000000101] = 8;
  firedLayersToQuality[0b000100000001110000] = 8;
  firedLayersToQuality[0b000001110000001101] = 8;
  firedLayersToQuality[0b000000110110000001] = 8;
  firedLayersToQuality[0b000001110000000001] = 8;
  firedLayersToQuality[0b001000010001000001] = 8;
  firedLayersToQuality[0b000001100000000101] = 8;
  firedLayersToQuality[0b000001100000000001] = 8;
  firedLayersToQuality[0b000001110000000101] = 8;
  firedLayersToQuality[0b001001110001000001] = 8;
  firedLayersToQuality[0b000010110000000101] = 8;
  firedLayersToQuality[0b000000010001000001] = 8;
  firedLayersToQuality[0b000000100110000001] = 8;
  firedLayersToQuality[0b001001100000001100] = 8;
  firedLayersToQuality[0b000001010000000001] = 8;
  firedLayersToQuality[0b000010100000000011] = 8;
  firedLayersToQuality[0b000000100001000001] = 8;
  firedLayersToQuality[0b001000110001000001] = 8;
  firedLayersToQuality[0b000000010010000001] = 8;
  firedLayersToQuality[0b000001010000000101] = 8;
  firedLayersToQuality[0b100000100110000001] = 8;
  firedLayersToQuality[0b000010010000000101] = 8;
  firedLayersToQuality[0b000000110010000001] = 8;
  firedLayersToQuality[0b000000000000110100] = 8;
  firedLayersToQuality[0b000000010000000101] = 8;
  firedLayersToQuality[0b000000110001000001] = 8;
  firedLayersToQuality[0b000000010000001100] = 8;
  firedLayersToQuality[0b000010110000001101] = 8;
  firedLayersToQuality[0b000011010000001101] = 8;
  firedLayersToQuality[0b000000100000010001] = 8;
  firedLayersToQuality[0b000000110000000101] = 8;
  firedLayersToQuality[0b000001100000000111] = 8;
  firedLayersToQuality[0b000000100000000101] = 8;
  firedLayersToQuality[0b010000010010000001] = 8;
  firedLayersToQuality[0b000001100000001101] = 8;
  firedLayersToQuality[0b000011100000000111] = 8;
  firedLayersToQuality[0b000000010110000001] = 8;
  firedLayersToQuality[0b000011110000000111] = 8;
  firedLayersToQuality[0b000000011100000000] = 8;
  firedLayersToQuality[0b001000010000000011] = 8;
  firedLayersToQuality[0b000001110000000011] = 8;
  firedLayersToQuality[0b000100000000110000] = 8;
  firedLayersToQuality[0b000111100000110100] = 8;
  firedLayersToQuality[0b010000010010000000] = 8;
  firedLayersToQuality[0b100000010100000000] = 8;
  firedLayersToQuality[0b001000100000000011] = 8;
  firedLayersToQuality[0b000011100000001101] = 8;
  firedLayersToQuality[0b100000011100000000] = 8;
  firedLayersToQuality[0b110000011110000001] = 8;
  //firedLayersToQuality[0b000000000000110011] = 8;
  //firedLayersToQuality[0b000000100110000011] = 8;
  //firedLayersToQuality[0b110000000100000000] = 8;
  //firedLayersToQuality[0b001011110001001101] = 8;
  //firedLayersToQuality[0b010000100001000011] = 8;
  //firedLayersToQuality[0b000001100000001100] = 8;
  //firedLayersToQuality[0b000001110001000011] = 8;
  //firedLayersToQuality[0b011000000010000000] = 8;
  //firedLayersToQuality[0b001000110100000011] = 8;
  //firedLayersToQuality[0b010001000011000000] = 8;
  //firedLayersToQuality[0b100000000110000000] = 8;
  //firedLayersToQuality[0b000000000000111100] = 8;
}

void OmtfEmulation::addObservers(const MuonGeometryTokens& muonGeometryTokens,
                                 const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord>& magneticFieldEsToken,
                                 const edm::ESGetToken<Propagator, TrackingComponentsRecord>& propagatorEsToken) {
  if (observers.empty()) {  //assuring it is done only at the first run
    OMTFReconstruction::addObservers(muonGeometryTokens, magneticFieldEsToken, propagatorEsToken);
    /*    if(edmParameterSet.exists("patternsPtAssignment") && edmParameterSet.getParameter<bool>("patternsPtAssignment")) {
      //std::string rootFileName = edmParameterSet.getParameter<std::string>("dumpHitsFileName");
      .emplace_back(std::make_unique<PatternsPtAssignment>(edmParameterSet, omtfConfig.get(), omtfProcGoldenPat->getPatterns(), ""));
    }*/
  }

  //addObservers is called in OMTFReconstruction::beginRun after the omtfProc is constructed, therefore here we can used omtfProc
  if (edmParameterSet.exists("neuralNetworkFile") && !ptAssignment) {
    edm::LogImportant("OMTFReconstruction") << "constructing PtAssignmentNNRegression" << std::endl;
    std::string neuralNetworkFile = edmParameterSet.getParameter<edm::FileInPath>("neuralNetworkFile").fullPath();
    ptAssignment = std::make_unique<PtAssignmentNNRegression>(edmParameterSet, omtfConfig.get(), neuralNetworkFile);
  }

  auto omtfProcGoldenPat = dynamic_cast<OMTFProcessor<GoldenPattern>*>(omtfProc.get());
  if (omtfProcGoldenPat) {
    omtfProcGoldenPat->setPtAssignment(ptAssignment.get());
    //omtfProcGoldenPat can be constructed from scratch each run, so ptAssignment is set herer every run
  }

  //TODO un-comment when convertToOuputScalesPhase2 is implemented
  omtfProc->setOutpuConversionFunction([&](unsigned int iProcessor, l1t::tftype mtfType, const AlgoMuons& gbCandidates) {
    return this->convertToOuputScalesPhase2(iProcessor, mtfType, gbCandidates);
  }); 
}

void OmtfEmulation::getQualityFromFiredLayers(FinalMuon& finalMuon) {
  auto it = firedLayersToQuality.find(finalMuon.getAlgoMuon()->getFiredLayerBits());
  if (it != firedLayersToQuality.end()) {
    finalMuon.setQuality(it->second);
  } else {
    finalMuon.setQuality(12);  //default value
  }
};

FinalMuons OmtfEmulation::convertToOuputScalesPhase2(unsigned int iProcessor,
                                                     l1t::tftype mtfType,
                                                     const AlgoMuons& gbCandidates) {
  FinalMuons finalMuons;
  for (auto& myCand : gbCandidates) {
    FinalMuon finalMuon(myCand);

    if (myCand->getPdfSumConstr() > 0 && myCand->getFiredLayerCntConstr() >= 3)
      finalMuon.setPt(round(myCand->getPtConstr()* 0.5 / Phase2L1GMT::LSBpt));
    else if (myCand->getPtUnconstr() > 0)
      finalMuon.setPt(round(1 * 0.5 / Phase2L1GMT::LSBpt));
    else
      finalMuon.setPt(0);

    if (finalMuon.getPt() == 0)
      continue;
    
    finalMuon.setSign(myCand->getChargeConstr()<0 ? 1 : 0);

    int etaValue = myCand->getEtaHw();
    if (mtfType == l1t::omtf_pos) {
      finalMuon.setEta(etaValue* 0.010875 / Phase2L1GMT::LSBeta);
    }
    else   {
      finalMuon.setEta((-1)*etaValue * 0.010875 / Phase2L1GMT::LSBeta);
    }

    int phiValue = myCand->getPhi();
    if (phiValue >= int(this->omtfConfig->nPhiBins()))
      phiValue -= this->omtfConfig->nPhiBins();
    //new the else cond.
    //phiValue = floor(phiValue * 437. / (1 << 12));
    float LSBphi_omtf = 2 * M_PI / this->omtfConfig->nPhiBins();
    int globPhi = iProcessor * 1800 + phiValue; // this is the phi in the OMTF in global coordinates 
    // first processor starts at CMS phi = 15 degrees (225 in int)... Handle wrap-around with %. Add 5400 to make sure the number is positive
    globPhi = (globPhi + 5400) % this->omtfConfig->nPhiBins(); 
    // convert to GMT Phi (2*pi / 2^13)
    finalMuon.setPhi(round(phiValue * LSBphi_omtf / Phase2L1GMT::LSBphi));

    if (myCand->getPtUnconstr() >= 0) {
      finalMuon.setPtUnconstr(round(myCand->getPtUnconstr() * 1.0 / Phase2L1GMT::LSBpt) );
    } else {
      finalMuon.setPtUnconstr(0);
    }
    if (ptAssignment) {
      finalMuon.setPt(round(myCand->getPtNNConstr() * 0.5 / Phase2L1GMT::LSBpt));
      finalMuon.setPtUnconstr(round(myCand->getPtNNUnconstr()* 1.0 / Phase2L1GMT::LSBpt));
      finalMuon.setSign(myCand->getChargeNNConstr() < 0 ? 1 : 0);
      finalMuon.setQuality(myCand->getQualityNN());
    }
    finalMuons.push_back(finalMuon);
  }
  return finalMuons;
}


l1t::SAMuonCollection OmtfEmulation::getSAMuons(unsigned int iProcessor,
                                                l1t::tftype mtfType,
                                                FinalMuons& finalMuons,
                                                bool uncostrainedPt) {
  l1t::SAMuonCollection saMuons;

  for (auto& finalMuon : finalMuons) {
    unsigned int qual = finalMuon.getQuality();
    int charge = finalMuon.getSign();
    unsigned int pt = uncostrainedPt ? finalMuon.getPtUnconstr() : finalMuon.getPt();
    int eta = finalMuon.getEta();
    int phi = finalMuon.getPhi();
    int z0 = 0;  // No tracks info
    int d0 = finalMuon.getHwD0();
  
    // Calculate Lorentz Vector
    math::PtEtaPhiMLorentzVector p4(pt * Phase2L1GMT::LSBpt, eta * Phase2L1GMT::LSBeta, phi * Phase2L1GMT::LSBphi, 0.0);
    l1t::SAMuon saMuon(p4, charge, pt, eta, phi, z0, d0, qual);
    saMuon.setTF(mtfType);
    //samuon.setWord(word);

    if (saMuon.hwPt() > 0) {
      saMuons.push_back(saMuon);
    }
  }

  return saMuons;
}

std::unique_ptr<l1t::SAMuonCollection> OmtfEmulation::run(
    const edm::Event& iEvent,
    const edm::EventSetup& evSetup,
    std::unique_ptr<l1t::RegionalMuonCandBxCollection>& candidates) {
  LogTrace("l1tOmtfEventPrint") << "\n" << __FUNCTION__ << ":" << __LINE__ << " iEvent " << iEvent.id().event() << endl;
  inputMaker->loadAndFilterDigis(iEvent);

  for (auto& obs : observers) {
    obs->observeEventBegin(iEvent);
  }

  std::unique_ptr<l1t::SAMuonCollection> saMuons = std::make_unique<l1t::SAMuonCollection>();
  candidates->setBXRange(bxMin, bxMax);

  ///The order is important: first put omtf_pos candidates, then omtf_neg.
  for (int bx = bxMin; bx <= bxMax; bx++) {
    for (unsigned int iProcessor = 0; iProcessor < omtfConfig->nProcessors(); ++iProcessor) {
      FinalMuons finalMuons = omtfProc->run(iProcessor, l1t::tftype::omtf_pos, bx, inputMaker.get(), observers);

      l1t::SAMuonCollection procSAMuons = getSAMuons(iProcessor, l1t::tftype::omtf_pos, finalMuons, false);

      //fill outgoing collection
      for (auto& saMuon : procSAMuons) {
        saMuons->push_back(saMuon);
      }

      std::vector<l1t::RegionalMuonCand> candMuons =
          omtfProc->getRegionalMuonCands(iProcessor, l1t::tftype::omtf_pos, finalMuons);
      for (auto& candMuon : candMuons) {
        candidates->push_back(bx, candMuon);
      }
    }

    for (unsigned int iProcessor = 0; iProcessor < omtfConfig->nProcessors(); ++iProcessor) {
      FinalMuons finalMuons = omtfProc->run(iProcessor, l1t::tftype::omtf_neg, bx, inputMaker.get(), observers);

      l1t::SAMuonCollection procSAMuons = getSAMuons(iProcessor, l1t::tftype::omtf_neg, finalMuons, false);
      //fill outgoing collection
      for (auto& saMuon : procSAMuons) {
        saMuons->push_back(saMuon);
      }

      std::vector<l1t::RegionalMuonCand> candMuons =
          omtfProc->getRegionalMuonCands(iProcessor, l1t::tftype::omtf_neg, finalMuons);
      for (auto& candMuon : candMuons) {
        candidates->push_back(bx, candMuon);
      }
    }

    //edm::LogInfo("OMTFReconstruction") <<"OMTF:  Number of candidates in BX="<<bx<<": "<<candidates->size(bx) << std::endl;;
  }

  LogTrace("l1tOmtfEventPrint") << __FUNCTION__ << ":" << __LINE__ << endl;
  for (auto& obs : observers) {
    obs->observeEventEnd(iEvent, candidates);
  }

  return saMuons;
}
