/*
 * OMTFInputDumper.cc
 *
 *  Created on: Apr 8, 2025
 *      Author: kbunkow, folguera
 */

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/OMTFInputDumper.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "TFile.h"
#include "TTree.h"

OMTFInputDumper::OMTFInputDumper(const edm::ParameterSet& edmCfg,
                                 const OMTFConfiguration* omtfConfig,
                                 CandidateSimMuonMatcher* candidateSimMuonMatcher)
    : EmulationObserverBase(edmCfg, omtfConfig), candidateSimMuonMatcher(candidateSimMuonMatcher), inputInProcs(omtfConfig->processorCnt())
    {
  edm::LogVerbatim("l1tOmtfEventPrint") << " omtfConfig->nTestRefHits() " << omtfConfig->nTestRefHits()
                                        << " event.omtfGpResultsPdfSum.num_elements() " << endl;

  initializeTTree();

  if (edmCfg.exists("dumpKilledOmtfCands"))
    if (edmCfg.getParameter<bool>("dumpKilledOmtfCands"))
      dumpKilledOmtfCands = true;

  if (edmCfg.exists("candidateSimMuonMatcherType")) {
    if (edmCfg.getParameter<std::string>("candidateSimMuonMatcherType") == "propagation")
      usePropagation = true;
    else if (edmCfg.getParameter<std::string>("candidateSimMuonMatcherType") == "matchSimple")
      usePropagation = false;

    edm::LogImportant("l1tOmtfEventPrint")
        << " CandidateSimMuonMatcher: candidateSimMuonMatcherType "
        << edmCfg.getParameter<std::string>("candidateSimMuonMatcherType") << std::endl;
  }

  edm::LogVerbatim("l1tOmtfEventPrint") << " OMTFInputDumper created. dumpKilledOmtfCands " << dumpKilledOmtfCands
                                        << std::endl;
}

OMTFInputDumper::~OMTFInputDumper() {}

void OMTFInputDumper::initializeTTree() {
  edm::Service<TFileService> fs;

  rootTree = fs->make<TTree>("OMTFHitsTree", "");

  rootTree->Branch("eventNum", &omtfEvent.eventNum);
  rootTree->Branch("muonEvent", &omtfEvent.muonEvent);

  rootTree->Branch("muonPt", &omtfEvent.muonPt);
  rootTree->Branch("muonEta", &omtfEvent.muonEta);
  rootTree->Branch("muonPhi", &omtfEvent.muonPhi);
  rootTree->Branch("muonPropEta", &omtfEvent.muonPropEta);
  rootTree->Branch("muonPropPhi", &omtfEvent.muonPropPhi);
  rootTree->Branch("muonCharge", &omtfEvent.muonCharge);

  rootTree->Branch("muonDxy", &omtfEvent.muonDxy);
  rootTree->Branch("muonRho", &omtfEvent.muonRho);

  rootTree->Branch("omtfPt", &omtfEvent.omtfPt);
  rootTree->Branch("omtfUPt", &omtfEvent.omtfUPt);
  rootTree->Branch("omtfEta", &omtfEvent.omtfEta);
  rootTree->Branch("omtfPhi", &omtfEvent.omtfPhi);
  rootTree->Branch("omtfCharge", &omtfEvent.omtfCharge);

  rootTree->Branch("omtfHwEta", &omtfEvent.omtfHwEta);

  rootTree->Branch("omtfProcessor", &omtfEvent.omtfProcessor);
  rootTree->Branch("omtfScore", &omtfEvent.omtfScore);
  rootTree->Branch("omtfQuality", &omtfEvent.omtfQuality);
  rootTree->Branch("omtfRefLayer", &omtfEvent.omtfRefLayer);
  rootTree->Branch("omtfRefHitNum", &omtfEvent.omtfRefHitNum);

  rootTree->Branch("omtfFiredLayers", &omtfEvent.omtfFiredLayers);  //<<<<<<<<<<<<<<<<<<<<<<!!!!TODOO

  rootTree->Branch("killed", &omtfEvent.killed);

  rootTree->Branch("stubNo", &omtfEvent.stubNo);
  rootTree->Branch("stubLayer", &omtfEvent.stubLayer);
  rootTree->Branch("stubQuality", &omtfEvent.stubQuality);
  rootTree->Branch("stubZ", &omtfEvent.stubZ);
  rootTree->Branch("stubValid", &omtfEvent.stubValid);
  rootTree->Branch("stubEta", &omtfEvent.stubEta);
  rootTree->Branch("stubPhiDist", &omtfEvent.stubPhiDist);
  rootTree->Branch("stubPhi", &omtfEvent.stubPhi);
  rootTree->Branch("stubIsRefLayer", &omtfEvent.stubIsRefLayer);
  rootTree->Branch("stubTiming", &omtfEvent.stubTiming);
  rootTree->Branch("stubDetId", &omtfEvent.stubDetId);
  rootTree->Branch("stubBx", &omtfEvent.stubBx);
  rootTree->Branch("stubPhiB", &omtfEvent.stubPhiB);
  rootTree->Branch("stubR", &omtfEvent.stubR);
  rootTree->Branch("stubDeltaR", &omtfEvent.stubDeltaR);
  rootTree->Branch("stubType", &omtfEvent.stubType);

  //Now also inputStubs: 
  rootTree->Branch("inputStubNo", &omtfEvent.inputStubNo);
  rootTree->Branch("inputStubLayer", &omtfEvent.inputStubLayer);
  rootTree->Branch("inputStubQuality", &omtfEvent.inputStubQuality);
  rootTree->Branch("inputStubZ", &omtfEvent.inputStubZ);
  rootTree->Branch("inputStubProc", &omtfEvent.inputStubProc);
  rootTree->Branch("inputStubEta", &omtfEvent.inputStubEta);
  rootTree->Branch("inputStubPhi", &omtfEvent.inputStubPhi);
  rootTree->Branch("inputStubIsRefLayer", &omtfEvent.inputStubIsRefLayer);
  rootTree->Branch("inputStubTiming", &omtfEvent.inputStubTiming);
  rootTree->Branch("inputStubDetId", &omtfEvent.inputStubDetId);
  rootTree->Branch("inputStubBx", &omtfEvent.inputStubBx);
  rootTree->Branch("inputStubPhiB", &omtfEvent.inputStubPhiB);
  rootTree->Branch("inputStubR", &omtfEvent.inputStubR);
  rootTree->Branch("inputStubType", &omtfEvent.inputStubType);
  rootTree->Branch("inputStubIsMatched", &omtfEvent.inputStubIsMatched);

  rootTree->Branch("deltaEta", &omtfEvent.deltaEta);
  rootTree->Branch("deltaPhi", &omtfEvent.deltaPhi);

  ptGenPos = fs->make<TH1I>("ptGenPos", "ptGenPos, eta at vertex 0.8 - 1.24", 400, 0, 200);  //TODO
  ptGenNeg = fs->make<TH1I>("ptGenNeg", "ptGenNeg, eta at vertex 0.8 - 1.24", 400, 0, 200);
}
void OMTFInputDumper::observeEventBegin(const edm::Event& iEvent) {
  clearOmtfInputStubs();
  clearOmtfStubs();
  for (auto& input: inputInProcs)
    input.reset();
}
void OMTFInputDumper::observeProcesorEmulation(unsigned int iProcessor,
                                               l1t::tftype mtfType,
                                               const std::shared_ptr<OMTFinput>& input,
                                               const AlgoMuons& algoCandidates,
                                               const AlgoMuons& gbCandidates,
                                               const FinalMuons& finalMuons) {

  unsigned int procIndx = omtfConfig->getProcIndx(iProcessor, mtfType);
  inputInProcs[procIndx] = input;
}

void OMTFInputDumper::observeEventEnd(const edm::Event& iEvent,
                                      std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) {
  /*
  int muonCharge = 0;
  if (simMuon) {
    if (std::abs(simMuon->momentum().eta()) < 0.8 || std::abs(simMuon->momentum().eta()) > 1.24)
      return;

    muonCharge = (std::abs(simMuon->type()) == 13) ? simMuon->type() / -13 : 0;
    if (muonCharge > 0)
      ptGenPos->Fill(simMuon->momentum().pt());
    else
      ptGenNeg->Fill(simMuon->momentum().pt());
  }

  if (simMuon == nullptr || !omtfCand->isValid())  //no sim muon or empty candidate
    return;

  omtfEvent.muonPt = simMuon->momentum().pt();
  omtfEvent.muonEta = simMuon->momentum().eta();

  //TODO add cut on ete if needed
    if(std::abs(event.muonEta) < 0.8 || std::abs(event.muonEta) > 1.24)
    return;

  omtfEvent.muonPhi = simMuon->momentum().phi();
  omtfEvent.muonCharge = muonCharge;  //TODO
   */

  std::vector<MatchingResult> matchingResults = candidateSimMuonMatcher->getMatchingResults();
  LogTrace("l1tOmtfEventPrint") << "\nOMTFInputDumper::observeEventEnd matchingResults.size() "
                                << matchingResults.size() << std::endl;

  //candidateSimMuonMatcher should use the  trackingParticles, because the simTracks are not stored for the pile-up events

  //for some events there are more than one matchingResults,
  //Usually at least one them has  genPt 0, which means no simMuon was matched, so candidate is ghost (or fake)
  //so better is to to drop such event, as it is not sure if the correct simMuon was matched to the candidate.
  //So we assume here that when the propagation is not used it is a single mu sample and this filter has sense
  //the propagation is used for multi-muon sample, so then this filter cannot be used
  //TODO add a flag to enable this filter? Disable it if not needed
  if (!usePropagation && matchingResults.size() > 1) {  //omtfConfig->cleanStubs() &&
    edm::LogVerbatim("l1tOmtfEventPrint")
        << "\nOMTFInputDumper::observeEventEnd matchingResults.size() " << matchingResults.size() << std::endl;

    for (auto& matchingResult : matchingResults) {
      edm::LogVerbatim("l1tOmtfEventPrint") << "matchingResult: genPt " << matchingResult.genPt;
      if (matchingResult.procMuon)
        edm::LogVerbatim("l1tOmtfEventPrint")
            << " procMuon.PtConstr " << matchingResult.procMuon->getPtConstr() << " processor "
            << matchingResult.muonCand->processor() << " hwPhi " << matchingResult.muonCand->hwPhi();
      else
        edm::LogVerbatim("l1tOmtfEventPrint") << " no procMuon" << std::endl;
    }
    edm::LogVerbatim("l1tOmtfEventPrint") << "dropping the event!!!\n" << std::endl;
    return;
  }

  for (auto& matchingResult : matchingResults) {
    omtfEvent.eventNum = iEvent.id().event();

    if (matchingResult.trackingParticle) {
      auto trackingParticle = matchingResult.trackingParticle;

      if (matchingResult.result == MatchingResult::ResultType::propagationFailed)
        omtfEvent.muonEvent = -2;
      else
        omtfEvent.muonEvent = trackingParticle->eventId().event();

      omtfEvent.muonPt = trackingParticle->pt();
      omtfEvent.muonEta = trackingParticle->momentum().eta();
      omtfEvent.muonPhi = trackingParticle->momentum().phi();
      omtfEvent.muonPropEta = matchingResult.propagatedEta;
      omtfEvent.muonPropPhi = matchingResult.propagatedPhi;
      omtfEvent.muonCharge = (std::abs(trackingParticle->pdgId()) == 13) ? trackingParticle->pdgId() / -13 : 0;

      if (trackingParticle->parentVertex().isNonnull()) {
        omtfEvent.muonDxy = trackingParticle->dxy();
        omtfEvent.muonRho = trackingParticle->parentVertex()->position().Rho();
      }

      omtfEvent.deltaEta = matchingResult.deltaEta;
      omtfEvent.deltaPhi = matchingResult.deltaPhi;

      LogTrace("l1tOmtfEventPrint") << "OMTFInputDumper::observeEventEnd trackingParticle: eventId "
                                    << trackingParticle->eventId().event() << " pdgId " << std::setw(3)
                                    << trackingParticle->pdgId() << " trackId "
                                    << trackingParticle->g4Tracks().at(0).trackId() << " pt " << std::setw(9)
                                    << trackingParticle->pt()  //<<" Beta "<<simMuon->momentum().Beta()
                                    << " eta " << std::setw(9) << trackingParticle->momentum().eta() << " phi "
                                    << std::setw(9) << trackingParticle->momentum().phi() << std::endl;

      if (std::abs(omtfEvent.muonEta) > 0.8 && std::abs(omtfEvent.muonEta) < 1.24) {
        if (omtfEvent.muonCharge > 0)
          ptGenPos->Fill(omtfEvent.muonPt);
        else
          ptGenNeg->Fill(omtfEvent.muonPt);
      }
    } else if (matchingResult.simTrack) {
      auto simTrack = matchingResult.simTrack;
      if (matchingResult.result == MatchingResult::ResultType::propagationFailed)
        omtfEvent.muonEvent = -2;
      else
        omtfEvent.muonEvent = simTrack->eventId().event();

      omtfEvent.muonPt = simTrack->momentum().pt();
      omtfEvent.muonEta = simTrack->momentum().eta();
      omtfEvent.muonPhi = simTrack->momentum().phi();
      omtfEvent.muonPropEta = matchingResult.propagatedEta;
      omtfEvent.muonPropPhi = matchingResult.propagatedPhi;
      omtfEvent.muonCharge = simTrack->charge();

      if (!simTrack->noVertex() && matchingResult.simVertex) {
        const math::XYZTLorentzVectorD& vtxPos = matchingResult.simVertex->position();
        omtfEvent.muonDxy = (-vtxPos.X() * simTrack->momentum().py() + vtxPos.Y() * simTrack->momentum().px()) /
                            simTrack->momentum().pt();
        omtfEvent.muonRho = vtxPos.Rho();
      }

      omtfEvent.deltaEta = matchingResult.deltaEta;
      omtfEvent.deltaPhi = matchingResult.deltaPhi;

      LogTrace("l1tOmtfEventPrint") << "OMTFInputDumper::observeEventEnd simTrack: eventId "
                                    << simTrack->eventId().event() << " pdgId " << std::setw(3)
                                    << simTrack->type()  //<< " trackId " << simTrack->g4Tracks().at(0).trackId()
                                    << " pt " << std::setw(9)
                                    << simTrack->momentum().pt()  //<<" Beta "<<simMuon->momentum().Beta()
                                    << " eta " << std::setw(9) << simTrack->momentum().eta() << " phi " << std::setw(9)
                                    << simTrack->momentum().phi() << std::endl;

      if (std::abs(omtfEvent.muonEta) > 0.8 && std::abs(omtfEvent.muonEta) < 1.24) {
        if (omtfEvent.muonCharge > 0)
          ptGenPos->Fill(omtfEvent.muonPt);
        else
          ptGenNeg->Fill(omtfEvent.muonPt);
      }
    } else {
      omtfEvent.muonEvent = -1;

      omtfEvent.muonPt = 0;

      omtfEvent.muonEta = 0;
      omtfEvent.muonPhi = 0;

      omtfEvent.muonPropEta = 0;
      omtfEvent.muonPropPhi = 0;

      omtfEvent.muonCharge = 0;  //TODO

      omtfEvent.muonDxy = 0;
      omtfEvent.muonRho = 0;
    }

    auto addOmtfCand = [&](AlgoMuonPtr& procMuon) {
      //the charge is only for the constrained measurement. The constrained measurement is always defined for a valid candidate
      if (procMuon->getPdfSumConstr() > 0 && procMuon->getFiredLayerCntConstr() >= 3)
        omtfEvent.omtfPt = omtfConfig->hwPtToGev(procMuon->getPtConstr());
      else if (procMuon->getPtUnconstr() > 0)
        //if myCand->getPdfSumConstr() == 0, the myCand->getPtConstr() might not be 0, see the end of GhostBusterPreferRefDt::select
        //but hwPt=0 means empty candidate, hwPt=1 maens pt=0,
        //but omtfPt = 0 means empty candidate
        //therefore here we set omtfPt=0.5 GeV, as the PtUnconstr > 0
        //N.B it is different than in the OMTFProcessor<GoldenPatternType>::convertToOuputScalesPhase1, where hwPt=1
        omtfEvent.omtfPt = 0.5;
      else
        omtfEvent.omtfPt = omtfConfig->hwPtToGev(0);

      //for candidate with no unconstrained measurement, hardware upt = 0
      //so then omtfEvent.omtfUPt is -1
      omtfEvent.omtfUPt = omtfConfig->hwUPtToGev(procMuon->getPtUnconstr());
      omtfEvent.omtfEta = omtfConfig->hwEtaToEta(procMuon->getEtaHw());
      omtfEvent.omtfPhi = procMuon->getPhi();
      omtfEvent.omtfCharge = procMuon->getChargeConstr();
      omtfEvent.omtfScore = procMuon->getPdfSum();

      omtfEvent.omtfHwEta = procMuon->getEtaHw();

      omtfEvent.omtfFiredLayers = procMuon->getFiredLayerBits();
      omtfEvent.omtfRefLayer = procMuon->getRefLayer();
      omtfEvent.omtfRefHitNum = procMuon->getRefHitNumber();

      clearOmtfStubs();

      //TODO choose, which gpResult should be dumped
      //auto& gpResult = procMuon->getGpResultConstr();
      auto& gpResult = (procMuon->getGpResultUnconstr().getPdfSumUnconstr() > procMuon->getGpResultConstr().getPdfSum())
                           ? procMuon->getGpResultUnconstr()
                           : procMuon->getGpResultConstr();

      /*
        edm::LogVerbatim("l1tOmtfEventPrint")<<"OMTFInputDumper:;observeEventEnd muonPt "<<event.muonPt<<" muonCharge "<<event.muonCharge
            <<" omtfPt "<<event.omtfPt<<" RefLayer "<<event.omtfRefLayer<<" omtfPtCont "<<event.omtfPtCont
            <<std::endl;  */

      for (unsigned int iLogicLayer = 0; iLogicLayer < gpResult.getStubResults().size(); ++iLogicLayer) {
        auto& stubResult = gpResult.getStubResults()[iLogicLayer];

        //TODO it is to have the hit if it is below the quality cut
        /*if (omtfConfig->isBendingLayer(iLogicLayer) && !stubResult.getMuonStub()) {
          auto& stubResult = gpResult.getStubResults()[iLogicLayer-1];
        }*/

        if (stubResult.getMuonStub()) {  //&& stubResult.getValid() //TODO!!!!!!!!!!!!!!!!1
          omtfEvent.stubNo++;
          omtfEvent.stubLayer.push_back(iLogicLayer);
          omtfEvent.stubQuality.push_back(stubResult.getMuonStub()->qualityHw);
          omtfEvent.stubValid.push_back(stubResult.getValid());

          unsigned int refLayerLogicNum = omtfConfig->getRefToLogicNumber()[procMuon->getRefLayer()];

          omtfEvent.stubPhi.push_back(stubResult.getMuonStub()->phiHw);
          omtfEvent.stubEta.push_back(stubResult.getMuonStub()->etaHw);
          omtfEvent.stubPhiDist.push_back(stubResult.getDeltaPhi());
          omtfEvent.stubIsRefLayer.push_back(iLogicLayer == refLayerLogicNum);
          omtfEvent.stubTiming.push_back(stubResult.getMuonStub()->timing);
          omtfEvent.stubDetId.push_back(stubResult.getMuonStub()->detId);
          omtfEvent.stubBx.push_back(stubResult.getMuonStub()->bx);
          omtfEvent.stubPhiB.push_back(stubResult.getMuonStub()->phiBHw);
          omtfEvent.stubR.push_back(stubResult.getMuonStub()->r);
          if (refLayerLogicNum == iLogicLayer)
            omtfEvent.stubDeltaR.push_back(stubResult.getMuonStub()->r - 413);
          else
            omtfEvent.stubDeltaR.push_back(stubResult.getMuonStub()->r - gpResult.getStubResults()[refLayerLogicNum].getMuonStub()->r);
          omtfEvent.stubType.push_back(stubResult.getMuonStub()->type);

          DetId detId(stubResult.getMuonStub()->detId);
          if (detId.subdetId() == MuonSubdetId::CSC) {
            CSCDetId cscId(detId);
            omtfEvent.stubZ.push_back(cscId.chamber() % 2);
          }
          //edm::LogVerbatim("l1tOmtfEventPrint")<<" hit.layer "<<(int)hit.layer<<" hit.phiDist "<<hit.phiDist<<" hit.rawData "<<hit.rawData << std::endl;
        }
      }

      LogTrace("l1tOmtfEventPrint") << "OMTFInputDumper::observeEventEnd adding omtfCand : " << std::endl;
      auto finalCandidate = matchingResult.muonCand;
      LogTrace("l1tOmtfEventPrint") << " hwPt " << finalCandidate->hwPt() << " hwSign " << finalCandidate->hwSign()
                                    << " hwQual " << finalCandidate->hwQual() << " hwEta " << std::setw(4)
                                    << finalCandidate->hwEta() << std::setw(4) << " hwPhi " << finalCandidate->hwPhi()
                                    << "    eta " << std::setw(9) << (finalCandidate->hwEta() * 0.010875)
                                    << " isKilled " << procMuon->isKilled() << " tRefLayer " << procMuon->getRefLayer()
                                    << " RefHitNumber " << procMuon->getRefHitNumber() << std::endl;
    };

    if (matchingResult.muonCand && matchingResult.procMuon->getPtConstr() > 0 &&
        matchingResult.muonCand->hwQual() >= 1) {
      //TODO set the quality, quality 0 has the candidates with eta > 1.3(?) EtaHw >= 121
      //&& matchingResult.genPt < 20

      omtfEvent.omtfQuality = matchingResult.muonCand->hwQual();  //procMuon->getQ();
      omtfEvent.killed = false;
      omtfEvent.omtfProcessor = matchingResult.muonCand->processor();

      if (matchingResult.muonCand->trackFinderType() == l1t::omtf_neg) {
        omtfEvent.omtfProcessor *= -1;
      }
      addOmtfInputStubsFromProc(matchingResult.muonCand->processor(), matchingResult.muonCand->trackFinderType(),
                                matchingResult.procMuon);
      addOmtfCand(matchingResult.procMuon);
      rootTree->Fill();
      clearOmtfStubs();
      clearOmtfInputStubs();

      if (dumpKilledOmtfCands) {
        for (auto& killedCand : matchingResult.procMuon->getKilledMuons()) {
          omtfEvent.omtfQuality = 0;
          omtfEvent.killed = true;
          if (killedCand->isKilled() == false) {
            edm::LogVerbatim("l1tOmtfEventPrint") << " killedCand->isKilled() == false !!!!!!!!";
          }
          addOmtfCand(killedCand);
          rootTree->Fill();
          clearOmtfStubs();
          clearOmtfInputStubs();
        }
      }
    } else if (omtfEvent.muonPt > 0) {  //checking if there was a simMuon
      LogTrace("l1tOmtfEventPrint") << "OMTFInputDumper::observeEventEnd no matching omtfCand" << std::endl;

      omtfEvent.omtfPt = 0;
      omtfEvent.omtfUPt = 0;
      omtfEvent.omtfEta = 0;
      omtfEvent.omtfPhi = 0;
      omtfEvent.omtfCharge = 0;
      omtfEvent.omtfScore = 0;

      omtfEvent.omtfHwEta = 0;

      omtfEvent.omtfFiredLayers = 0;
      omtfEvent.omtfRefLayer = 0;
      omtfEvent.omtfRefHitNum = 0;
      omtfEvent.omtfProcessor = 10;

      omtfEvent.omtfQuality = 0;
      omtfEvent.killed = false;

      clearOmtfStubs();
      clearOmtfInputStubs();
      rootTree->Fill();
    }
  }
  evntCnt++;
}
void OMTFInputDumper::addOmtfInputStubsFromProc(int iProc, l1t::tftype mtfType, AlgoMuonPtr& muon){
  int procIndx = omtfConfig->getProcIndx(iProc, mtfType);
  
  if (!inputInProcs[procIndx]) return; 
  auto& omtfInput = *inputInProcs[procIndx];

  for (unsigned int iLayer = 0; iLayer < omtfConfig->nLayers(); ++iLayer) {
    for (unsigned int iInput = 0; iInput < omtfInput.getMuonStubs()[iLayer].size(); ++iInput) {
      auto layer = iLayer;
      if  (omtfConfig->isBendingLayer(iLayer)) continue; //layer = iLayer - 1;

      auto stub = omtfInput.getMuonStub(layer, iInput);
      if (stub) {  
        omtfEvent.inputStubNo++;
        omtfEvent.inputStubLayer.push_back(stub->logicLayer);
        omtfEvent.inputStubProc.push_back(procIndx);
        omtfEvent.inputStubPhi.push_back(stub->phiHw);
        omtfEvent.inputStubPhiB.push_back(stub->phiBHw);
        omtfEvent.inputStubEta.push_back(stub->etaHw);
        omtfEvent.inputStubQuality.push_back(stub->qualityHw);
        omtfEvent.inputStubBx.push_back(stub->bx);
        omtfEvent.inputStubTiming.push_back(stub->timing);
        omtfEvent.inputStubDetId.push_back(stub->detId);
        omtfEvent.inputStubType.push_back(stub->type);
        omtfEvent.inputStubR.push_back(stub->r);
        omtfEvent.inputStubIsMatched.push_back(isMatchedStub(stub, muon));

        DetId detId(stub->detId);
        if (detId.subdetId() == MuonSubdetId::CSC) {
          CSCDetId cscId(detId);
          omtfEvent.inputStubZ.push_back(cscId.chamber() % 2);
        }


      }
    }
  }
}

bool OMTFInputDumper::isMatchedStub(const MuonStubPtr& stub, AlgoMuonPtr& procMuon) {
    if (stub->type == MuonStub::Type::EMPTY)
      return false;
  
    auto& gpResult = (procMuon->getGpResultUnconstr().getPdfSumUnconstr() > procMuon->getGpResultConstr().getPdfSum())
                      ? procMuon->getGpResultUnconstr() : procMuon->getGpResultConstr();
  
    for (unsigned int iLogicLayer = 0; iLogicLayer < gpResult.getStubResults().size(); ++iLogicLayer) {
      auto& stubResult = gpResult.getStubResults()[iLogicLayer];   
      if (stubResult.getMuonStub() && stubResult.getMuonStub()->detId == stub->detId) {
        return true;
      }
    }
    return false;
}
  

void OMTFInputDumper::clearOmtfStubs() {
  omtfEvent.stubNo = 0;
  omtfEvent.stubLayer.clear();
  omtfEvent.stubQuality.clear();
  omtfEvent.stubValid.clear();
  omtfEvent.stubZ.clear();
  omtfEvent.stubPhi.clear();
  omtfEvent.stubEta.clear();
  omtfEvent.stubPhiDist.clear();
  omtfEvent.stubIsRefLayer.clear();
  omtfEvent.stubTiming.clear();
  omtfEvent.stubDetId.clear();
  omtfEvent.stubBx.clear();
  omtfEvent.stubPhiB.clear();
  omtfEvent.stubR.clear();
  omtfEvent.stubDeltaR.clear();
  omtfEvent.stubType.clear();
}

void OMTFInputDumper::clearOmtfInputStubs(){
  omtfEvent.inputStubNo = 0;
  omtfEvent.inputStubLayer.clear();
  omtfEvent.inputStubQuality.clear();
  omtfEvent.inputStubProc.clear();
  omtfEvent.inputStubZ.clear();
  omtfEvent.inputStubPhi.clear();
  omtfEvent.inputStubEta.clear();
  omtfEvent.inputStubIsRefLayer.clear();
  omtfEvent.inputStubTiming.clear();
  omtfEvent.inputStubDetId.clear();
  omtfEvent.inputStubBx.clear();
  omtfEvent.inputStubPhiB.clear();
  omtfEvent.inputStubR.clear();
  omtfEvent.inputStubType.clear();
  omtfEvent.inputStubIsMatched.clear();
  
}
void OMTFInputDumper::endJob() {
  edm::LogVerbatim("l1tOmtfEventPrint") << " evntCnt " << evntCnt << endl;
  rootTree->Write();
}
