#include <iostream>
#include <strstream>
#include <iomanip>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/DataRecord/interface/L1TMuonOverlapParamsRcd.h"
#include "CondFormats/L1TObjects/interface/L1TMuonOverlapParams.h"

#include "L1Trigger/L1TMuonOverlap/plugins/OMTFHitAnalyzer.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinputMaker.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfigMaker.h"
#include "L1Trigger/L1TMuonOverlap/interface/XMLConfigWriter.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternPdf4D.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"

#include "SimDataFormats/Track/interface/SimTrack.h"

#include "Math/VectorUtil.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"

#include <TH2F.h>
#include "TFile.h"


OMTFHitAnalyzer::OMTFHitAnalyzer(const edm::ParameterSet& cfg):
theConfig(cfg),
g4SimTrackSrc(cfg.getParameter<edm::InputTag>("g4SimTrackSrc")){

  inputTokenDTPh = consumes<L1MuDTChambPhContainer>(theConfig.getParameter<edm::InputTag>("srcDTPh"));
  inputTokenDTTh = consumes<L1MuDTChambThContainer>(theConfig.getParameter<edm::InputTag>("srcDTTh"));
  inputTokenCSC = consumes<CSCCorrelatedLCTDigiCollection>(theConfig.getParameter<edm::InputTag>("srcCSC"));
  inputTokenRPC = consumes<RPCDigiCollection>(theConfig.getParameter<edm::InputTag>("srcRPC"));
  inputTokenSimHit = consumes<edm::SimTrackContainer>(theConfig.getParameter<edm::InputTag>("g4SimTrackSrc"));

  myInputMaker = new OMTFinputMaker();

  makeGoldenPatterns = theConfig.getParameter<bool>("makeGoldenPatterns");
  makeConnectionsMaps = theConfig.getParameter<bool>("makeConnectionsMaps");
  mergeXMLFiles = theConfig.getParameter<bool>("mergeXMLFiles");

  myOMTFConfig = 0;

  ptDist = new TH1I("ptDist", "ptDist", 200, -0.5, 200-0.5);
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
OMTFHitAnalyzer::~OMTFHitAnalyzer(){

  delete myOMTFConfig;
  //delete myOMTFConfigMaker;
}
/////////////////////////////////////////////////////
void OMTFHitAnalyzer::configureProcesor(const OMTFConfiguration * omtfConfig,
    const L1TMuonOverlapParams* omtfPatterns, unsigned int ptCode, int charge, unsigned int patNum) {
  //myOMTF->setConfigurataion(myOMTFConfig);

  //myResults.assign(omtfConfig->nTestRefHits(),OMTFProcessor::resultsMap()); FIXME is it needed???

/*  const l1t::LUT* chargeLUT =  omtfPatterns->chargeLUT();
  const l1t::LUT* etaLUT =  omtfPatterns->etaLUT();
  const l1t::LUT* ptLUT =  omtfPatterns->ptLUT();

  unsigned int nGPs = omtfConfig->nGoldenPatterns();
  unsigned int address = 0;*/
  unsigned int iEta = 0;
  unsigned int iPt = ptCode;
  int iCharge = charge;
  /*for(unsigned int iGP=0; iGP < nGPs; ++iGP) {
    address = iGP;
    iEta = etaLUT->data(address);
    iCharge = chargeLUT->data(address)==0 ? -1:1;
    iPt = ptLUT->data(address);
    std::cout<<"iGP "<<iGP<<" etaLUT "<<iEta<<" iCharge "<<iCharge<<" iPt "<<iPt<<std::endl;
    if(iPt == ptCode && iCharge == charge) { //TODO check the ptCode
      Key aKey(iEta,iPt,iCharge,iGP);
      std::cout<<"adding GoldenPatternPdf4D "<<aKey<<std::endl;
      GoldenPatternPdf4D *aGP = new GoldenPatternPdf4D(aKey, omtfConfig);
      aGP->reset();
      omtfProc->addGP(aGP);
    }
  }*/

  Key aKey(iEta,iPt,iCharge, patNum);
  std::cout<<"adding GoldenPatternPdf4D "<<aKey<<std::endl;
  GoldenPatternPdf4D* aGP = new GoldenPatternPdf4D(aKey, omtfConfig);
  std::cout<<"adding GoldenPatternPdf4D "<<aGP->key()<<std::endl;
  //aGP->reset();
  Pdf4DGeneratorProcessor::GoldenPatternVec gps;
  gps.emplace_back(aGP);

  myOMTF.reset(new Pdf4DGeneratorProcessor(myOMTFConfig, gps) );

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFHitAnalyzer::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {

  const L1TMuonOverlapParamsRcd& omtfParamsRcd = iSetup.get<L1TMuonOverlapParamsRcd>();

  edm::ESHandle<L1TMuonOverlapParams> omtfParamsHandle;
  omtfParamsRcd.get(omtfParamsHandle);

  const L1TMuonOverlapParams* omtfParams = omtfParamsHandle.product();

  if (!omtfParams) {
    edm::LogError("L1TMuonOverlapTrackProducer") << "Could not retrieve parameters from Event Setup" << std::endl;
  }

  ///Initialise XML writer with default pdf.
  myWriter = new XMLConfigWriter(myOMTFConfig);

  ///For making the patterns use extended pdf width in phi, as pdf are later shifted by the mean value
  ///For low pt muons non shifted pdfs would go out of the default pdf range.
  L1TMuonOverlapParams omtfParamsMutable = *omtfParams;
  std::vector<int> generalParams = *omtfParamsMutable.generalParams();
  nPdfAddrBits = omtfParams->nPdfAddrBits();

  //if(!mergeXMLFiles) //the LUTs are then too big in the GoldenPatternPdf4D
  generalParams[L1TMuonOverlapParams::GENERAL_ADDRBITS] = nPdfAddrBits + 2; //2*nPdfAddrBits;
  omtfParamsMutable.setGeneralParams(generalParams);

  myOMTFConfig->configure(&omtfParamsMutable);
  int ptCode = theConfig.getParameter<int>("ptCode"); //assuming that here the legacy PAC pt code is given
  int charge = theConfig.getParameter<int>("charge");

  //converting to the uGMT ptCode
  double pt = RPCConst::ptFromIpt(ptCode);
  unsigned int patNum = myOMTFConfig->getPatternNum(pt, charge);
  ptCode = omtfParams->ptLUT()->data(patNum );
  configureProcesor(myOMTFConfig, omtfParams, ptCode, charge, patNum);
  //myOMTFConfigMaker = new OMTFConfigMaker(myOMTFConfig);

  for(int iPt = 0; iPt <= 31; iPt++) {
    double pt = RPCConst::ptFromIpt(iPt);
    unsigned int patNum = myOMTFConfig->getPatternNum(pt, charge);
    ptCode = omtfParams->ptLUT()->data(patNum );
    std::cout<<" ipt "<<std::setw(3)<<iPt<<" pt "<<std::setw(3)<<pt<<" [GeV]"<<std::setw(3)<<" patNum "<<std::setw(3)<<patNum<<" ptCodeOmtf "<<std::setw(3)<<ptCode<<endl;
  }

  for(unsigned int iPat = 0; iPat < myOMTFConfig->nGoldenPatterns(); iPat++) {
    std::cout<<"pat num"<<std::setw(3)<<iPat<<" ptFrom "<<std::setw(3)<<myOMTFConfig->getPatternPtRange(iPat).ptFrom
        <<" ptFrom "<<std::setw(3)<<myOMTFConfig->getPatternPtRange(iPat).ptTo<<std::endl;
  }

  std::cout<<"OMTFHitAnalyzer::beginRun: myOMTFConfig "<<*myOMTFConfig;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFHitAnalyzer::beginJob(){

  myOMTFConfig = new OMTFConfiguration();
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////  
void OMTFHitAnalyzer::endJob(){
  std::cout<<"myOMTFConfig "<<*myOMTFConfig;

  //if(makeGoldenPatterns && !makeConnectionsMaps)
  {
    //myWriter->initialiseXMLDocument("OMTF");
/*
    for(auto itGP: myGPmap){
      if(!itGP->hasCounts())
        continue;
      itGP->normalise(nPdfAddrBits);
    }
*/

    ///Put back default value of the pdf width.
/*    L1TMuonOverlapParams omtfParamsMutable = *myOMTFConfig->getRawParams();
    std::vector<int> generalParams = *omtfParamsMutable.generalParams();
    generalParams[L1TMuonOverlapParams::GENERAL_ADDRBITS] = nPdfAddrBits;
    omtfParamsMutable.setGeneralParams(generalParams);
    myOMTFConfig->configure(&omtfParamsMutable);*/

    cout<<__FUNCTION__<<":"<<__LINE__<<"myOMTF->getPatterns().size()"<<myOMTF->getPatterns().size()<<" "<<endl;
    for(auto& itGP: myOMTF->getPatterns()) {
       /*      if(itGP.first.thePtCode==iPt &&
          itGP.first.theCharge==theConfig.getParameter<int>("charge")) {
        //std::cout<<*itGP.second<<std::endl; FIXME
        myWriter->writeGPData(*((GoldenPattern*)(itGP.second)), dummyGP, dummyGP, dummyGP);
      }*/

      cout<<__FUNCTION__<<":"<<__LINE__<<"itGP->key().thePt "<<itGP->key().thePt<<" "<<endl;
      if(itGP->key().thePt == 0)
        continue;

      std::ostringstream fileName;
      fileName<<"bendingDistr_ptCode_"<<itGP->key().thePt
          <<"_ch_"<<itGP->key().theCharge<<".root";

      cout<<"out fileName"<<fileName.str()<<endl;
      TFile* outfile = new TFile(fileName.str().c_str(), "RECREATE");
      cout<<"out fileName "<<fileName.str()<<" outfile->GetName() "<<outfile->GetName()<<endl;
      ptDist->Write();

      for(unsigned int iLayer = 0; iLayer<myOMTFConfig->nLayers(); ++iLayer) {
        for(unsigned int iRefLayer=0; iRefLayer<myOMTFConfig->nRefLayers(); ++iRefLayer) {
          std::ostringstream histName;

          histName<<"ipt_"<<itGP->key().thePt<<"_ch"<<itGP->key().theCharge<<"_layer_"<<iLayer<<"_refLayer_"<<iRefLayer<<" ";

          unsigned int refLayerPhiBSize = itGP->getPdf()[iLayer][iRefLayer].size();
          unsigned int layerPhiSize = (itGP)->getPdf()[iLayer][iRefLayer][0].size();
          cout<<"creating hist "<<histName.str()<<" refLayerPhiBSize "<<refLayerPhiBSize<<" layerPhiSize "<<layerPhiSize<<std::endl;

          if(refLayerPhiBSize == 1 ) {
            TH1F *h1 = new TH1F(histName.str().c_str(), histName.str().c_str(), layerPhiSize, -0.5, layerPhiSize-0.5);
            for(unsigned int iLayerPhi=0; iLayerPhi < layerPhiSize; iLayerPhi++) {
              h1->Fill(iLayerPhi, itGP->getPdf()[iLayer][iRefLayer][0][iLayerPhi]);
            }
            h1->Write();
          }
          else {
            TH2F *h2 = new TH2F(histName.str().c_str(), histName.str().c_str(), refLayerPhiBSize, -0.5, refLayerPhiBSize-0.5,
                layerPhiSize, -0.5, layerPhiSize-0.5);
            h2->GetXaxis()->SetTitle("refLayerPhiB");
            h2->GetYaxis()->SetTitle("layerDelatPhi");
            for(unsigned int iRefLayerPhiB = 0; iRefLayerPhiB < refLayerPhiBSize; iRefLayerPhiB++) {
              for(unsigned int iLayerPhi=0; iLayerPhi < layerPhiSize; iLayerPhi++) {
                h2->Fill(iRefLayerPhiB, iLayerPhi, itGP->getPdf()[iLayer][iRefLayer][iRefLayerPhiB][iLayerPhi]);
              }
            }
            h2->Write();
          }
        }
      }
      outfile->Close();
    }
  }


}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& evSetup){

  if(mergeXMLFiles) return;

  ///Get the simulated muon parameters
  const SimTrack* aSimMuon = findSimMuon(iEvent,evSetup);
  if(!aSimMuon){
    edm::LogError("OMTFHitAnalyzer")<<"No SimMuon found in the event!";
    return;
  }

  //std::cout<<"new event ";
  //cout<<"aSimMuon->momentum().pt "<<aSimMuon->momentum().pt()<<std::endl;

  ptDist->Fill(aSimMuon->momentum().pt());

  myInputMaker->initialize(evSetup, myOMTFConfig);

  edm::Handle<L1MuDTChambPhContainer> dtPhDigis;
  edm::Handle<L1MuDTChambThContainer> dtThDigis;
  edm::Handle<CSCCorrelatedLCTDigiCollection> cscDigis;
  edm::Handle<RPCDigiCollection> rpcDigis;

  ///Filter digis by dropping digis from selected (by cfg.py) subsystems
  if(!theConfig.getParameter<bool>("dropDTPrimitives")){
    iEvent.getByToken(inputTokenDTPh,dtPhDigis);
    iEvent.getByToken(inputTokenDTTh,dtThDigis);
  }
  if(!theConfig.getParameter<bool>("dropRPCPrimitives"))
    iEvent.getByToken(inputTokenRPC,rpcDigis);
  if(!theConfig.getParameter<bool>("dropCSCPrimitives"))
    iEvent.getByToken(inputTokenCSC,cscDigis);

  //l1t::tftype mtfType = l1t::tftype::bmtf;
  l1t::tftype mtfType = l1t::tftype::omtf_pos;
  //l1t::tftype mtfType = l1t::tftype::emtf_pos;

  ///Loop over all processors, each covering 60 deg in phi
  for(unsigned int iProcessor=0;iProcessor<6;++iProcessor) {
    ///Input data with phi ranges shifted for each processor, so it fits 11 bits range
    OMTFinput myInput = myInputMaker->buildInputForProcessor(dtPhDigis.product(),
        dtThDigis.product(),
        cscDigis.product(),
        rpcDigis.product(),
        iProcessor,
        mtfType);

    myOMTF->fillCounts(iProcessor, myInput, aSimMuon);

/*    std::ostringstream ostr;
    bool wasHit = false;
    for(unsigned int iLayer=0; iLayer < myOMTFConfig->nLayers(); ++iLayer){
      const OMTFinput::vector1D& layerHits = myInput.getLayerData(iLayer, false);
      ostr<<"layer "<<iLayer<<" ";
      for(unsigned int i = 0; i < layerHits.size(); i++) {
        if(layerHits[i] != 5400) {
          wasHit = true;
          ostr<<std::setw(4)<<layerHits[i]<<" ";
        }
        else
          ostr<<"     ";
      }
      ostr<<std::endl;
    }
    if(wasHit) {
      //std::cout<<"iProcessor "<<iProcessor<<std::endl<<ostr.str()<<std::endl;
    }*/

  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////  
const SimTrack * OMTFHitAnalyzer::findSimMuon(const edm::Event &ev, const edm::EventSetup &es, const SimTrack * previous){

  const SimTrack * result = 0;
  edm::Handle<edm::SimTrackContainer> simTks;
  ev.getByToken(inputTokenSimHit,simTks);

  for (std::vector<SimTrack>::const_iterator it=simTks->begin(); it< simTks->end(); it++) {
    const SimTrack & aTrack = *it;
    if ( !(aTrack.type() == 13 || aTrack.type() == -13) )continue;
    if(previous && ROOT::Math::VectorUtil::DeltaR(aTrack.momentum(),previous->momentum())<0.07) continue;
    if ( !result || aTrack.momentum().pt() > result->momentum().pt()) result = &aTrack;
  }
  return result;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////  
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(OMTFHitAnalyzer);
