/*
 * OMTFInputDumper.h
 *
 *  Created on: Apr 8, 2025
 *      Author: kbunkow, folguera
 */

#ifndef L1T_OmtfP1_TOOLS_OMTFInputDumper_H_
#define L1T_OmtfP1_TOOLS_OMTFInputDumper_H_

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/EmulationObserverBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/CandidateSimMuonMatcher.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "TMap.h"
#include "TArrayI.h"
#include "TFile.h"
#include "TH2.h"

#include <functional>

class TTree;

struct OmtfEventInput {
public:
  unsigned int eventNum = 0;

  //muonPt = 0 means that no muon was matched to the candidate
  short muonEvent = -1;
  float muonPt = 0, muonEta = 0, muonPhi = 0, muonPropEta = 0, muonPropPhi = 0;
  char muonCharge = 0;
  float muonDxy = 0;
  float muonRho = 0;

  float omtfPt = 0, omtfEta = 0, omtfPhi = 0, omtfUPt = 0;
  char omtfCharge = 0;
  char omtfProcessor = 0;
  short omtfScore = 0;

  short omtfHwEta = 0;

  char omtfQuality = 0;
  char omtfRefLayer = 0;
  char omtfRefHitNum = 0;

  unsigned int omtfFiredLayers = 0;

  bool killed = false;

  float deltaPhi = 0, deltaEta = 0;

  //float omtfPtCont = 0;
  int stubNo = 0;
  std::vector<unsigned int> stubLayer;
  std::vector<unsigned int> stubQuality;
  std::vector<int> stubZ, stubValid, stubEta, stubPhi,stubPhiB, stubR, stubPhiDist, stubDeltaR;
  std::vector<int> stubIsRefLayer;
  std::vector<int> stubBx, stubTiming;
  std::vector<int> stubDetId;
  std::vector<int> stubType;

  int inputStubNo = 0;
  std::vector<unsigned int> inputStubLayer;
  std::vector<unsigned int> inputStubQuality;
  std::vector<int> inputStubZ, inputStubEta, inputStubPhi,inputStubPhiB, inputStubR;
  std::vector<int> inputStubIsRefLayer, inputStubProc;
  std::vector<int> inputStubBx, inputStubTiming;
  std::vector<int> inputStubDetId;
  std::vector<int> inputStubType;
  std::vector<int> inputStubIsMatched;
  
};

class OMTFInputDumper : public EmulationObserverBase {
public:
  OMTFInputDumper(const edm::ParameterSet& edmCfg,
                  const OMTFConfiguration* omtfConfig,
                  CandidateSimMuonMatcher* candidateSimMuonMatcher);

  ~OMTFInputDumper() override;

  void observeProcesorEmulation(unsigned int iProcessor,
                                l1t::tftype mtfType,
                                const std::shared_ptr<OMTFinput>&,
                                const AlgoMuons& algoCandidates,
                                const AlgoMuons& gbCandidates,
                                const FinalMuons& finalMuons) override;
  
  void observeEventBegin(const edm::Event& iEvent) override;
  void observeEventEnd(const edm::Event& iEvent,
                       std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) override;

  // Methods for inputStubs...
  void addOmtfInputStubsFromProc(int iProc, l1t::tftype mtfType) {};
  void addOmtfInputStubsFromProc(int iProc, l1t::tftype mtfType, AlgoMuonPtr& procMuon);
  bool isMatchedStub(const MuonStubPtr& stub, AlgoMuonPtr& procMuon);
  bool isRefLayer(const MuonStubPtr& stub, AlgoMuonPtr& procMuon);
  void clearOmtfStubs();
  void clearOmtfInputStubs();
  
  void endJob() override;

private:
  void initializeTTree();

  CandidateSimMuonMatcher* candidateSimMuonMatcher = nullptr;

  TTree* rootTree = nullptr;

  OmtfEventInput omtfEvent;
  std::vector<std::shared_ptr<OMTFinput> > inputInProcs;
  unsigned int evntCnt = 0;

  TH1I* ptGenPos = nullptr;
  TH1I* ptGenNeg = nullptr;

  //std::vector<TH2*> hitVsPt;

  bool dumpKilledOmtfCands = false;

  bool usePropagation = false;
};

#endif /* L1T_OmtfP1_TOOLS_OMTFInputDumper_H_ */
