/*
 * ProcessorBase.cpp
 *
 *  Created on: Jul 28, 2017
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlap/interface/ProcessorBase.h"

#include "CondFormats/L1TObjects/interface/L1TMuonOverlapParams.h"
#include "SimDataFormats/Track/interface/SimTrack.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "L1Trigger/L1TMuonOverlap/interface/GoldenPattern.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternParametrised.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternWithStat.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternPdfGen.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternPdf4D.h"

/*ProcessorBase::ProcessorBase(): myOmtfConfig(0) {
  // TODO Auto-generated constructor stub
}*/


template<class GoldenPatternType>
void ProcessorBase<GoldenPatternType>::resetConfiguration() {
  //myResults.clear();
  //for(auto it: theGPs) delete it;
  theGPs.clear();
}

///////////////////////////////////////////////
///////////////////////////////////////////////
template<class GoldenPatternType>
bool ProcessorBase<GoldenPatternType>::configure(const OMTFConfiguration* omtfConfig,
    const L1TMuonOverlapParams* omtfPatterns){
  resetConfiguration();

  myOmtfConfig = omtfConfig;

  //myResults.assign(myOmtfConfig->nTestRefHits(),ProcessorBase::resultsMap());

  const l1t::LUT* chargeLUT =  omtfPatterns->chargeLUT();
  const l1t::LUT* etaLUT =  omtfPatterns->etaLUT();
  const l1t::LUT* ptLUT =  omtfPatterns->ptLUT();
  const l1t::LUT* pdfLUT =  omtfPatterns->pdfLUT();
  const l1t::LUT* meanDistPhiLUT =  omtfPatterns->meanDistPhiLUT();

  unsigned int nGPs = myOmtfConfig->nGoldenPatterns();
  edm::LogInfo("MTFProcessor::configure")<<"myOmtfConfig->nGoldenPatterns() "<<nGPs<<std::endl;
  unsigned int address = 0;
  unsigned int iEta, iPt;
  int iCharge;
  int meanDistPhiSize = myOmtfConfig->nLayers() * myOmtfConfig->nRefLayers() * myOmtfConfig->nGoldenPatterns();
  for(unsigned int iGP=0;iGP<nGPs;++iGP){
    address = iGP;
    iEta = etaLUT->data(address);
    iCharge = chargeLUT->data(address)==0? -1:1;
    iPt = ptLUT->data(address);

    Key aKey(iEta,iPt,iCharge,iGP);
    edm::LogInfo("ProcessorBase::configure")<<"adding pattern "<<aKey<<" "<<std::endl; //<<myOmtfConfig->getPatternPtRange(iGP).ptFrom<<" - "<<myOmtfConfig->getPatternPtRange(iGP).ptTo<<" GeV"<<std::endl; PatternPtRange is not initialized here yet!!!!
    GoldenPatternType* aGP = new GoldenPatternType(aKey, myOmtfConfig);

    ///Mean dist phi data
    for(unsigned int iLayer=0;iLayer<myOmtfConfig->nLayers();++iLayer){
      for(unsigned int iRefLayer=0;iRefLayer<myOmtfConfig->nRefLayers();++iRefLayer){
        address = iRefLayer + iLayer*myOmtfConfig->nRefLayers() + iGP*(myOmtfConfig->nRefLayers()*myOmtfConfig->nLayers());
        int value = meanDistPhiLUT->data(address) - (1<<(meanDistPhiLUT->nrBitsData() -1));
        aGP->setMeanDistPhiValue(value, iLayer, iRefLayer, 0);
        if(meanDistPhiLUT->nrBitsAddress() == 15) {//for the new version of the meanDistPhi which have two values for each gp,iLayer,iRefLayer, FIXME: do it a better way
          value = meanDistPhiLUT->data(address + meanDistPhiSize) - (1<<(meanDistPhiLUT->nrBitsData() -1));
          //the second meanDistPhi is in the LUT at the position (address+meanDistPhiSize)
          aGP->setMeanDistPhiValue(value, iLayer, iRefLayer, 1);
        }

        //TODO add handling of selDistPhiShift - requires changing in the L1TMuonOverlapParams
      }
      ///Pdf data
      for(unsigned int iRefLayer=0;iRefLayer<myOmtfConfig->nRefLayers();++iRefLayer){
        for(unsigned int iPdf=0;iPdf<(unsigned int)(1<<myOmtfConfig->nPdfAddrBits());++iPdf){
          address = iPdf + iRefLayer*(1<<myOmtfConfig->nPdfAddrBits()) +
              iLayer*myOmtfConfig->nRefLayers()*(1<<myOmtfConfig->nPdfAddrBits()) +
              iGP*myOmtfConfig->nLayers()*myOmtfConfig->nRefLayers()*(1<<myOmtfConfig->nPdfAddrBits());
          int value = pdfLUT->data(address);//here only int is possible
          aGP->setPdfValue(value, iLayer, iRefLayer, iPdf);
        }
      }
    }
    addGP(aGP);
  }

  initPatternPtRange();

  return true;
}

///////////////////////////////////////////////
///////////////////////////////////////////////
template<class GoldenPatternType>
void ProcessorBase<GoldenPatternType>::addGP(GoldenPatternType* aGP) {
  theGPs.emplace_back(std::unique_ptr<GoldenPatternType>(aGP));
}

////////////////////////////////////////////
////////////////////////////////////////////
template<class GoldenPatternType>
OMTFinput::vector1D ProcessorBase<GoldenPatternType>::restrictInput(unsigned int iProcessor,
    unsigned int iRegion,
    unsigned int iLayer,
    const OMTFinput::vector1D & layerHits) {

  OMTFinput::vector1D myHits = layerHits;

  unsigned int iStart = myOmtfConfig->getConnections()[iProcessor][iRegion][iLayer].first;
  unsigned int iEnd = iStart + myOmtfConfig->getConnections()[iProcessor][iRegion][iLayer].second -1;

  for(unsigned int iInput=0;iInput<myHits.size();++iInput){
    if(iInput<iStart || iInput>iEnd) myHits[iInput] = myOmtfConfig->nPhiBins();
  }
  return myHits;
}

////////////////////////////////////////////
////////////////////////////////////////////
template <class GoldenPatternType>
void ProcessorBase<GoldenPatternType>::initPatternPtRange() {
  patternPts.clear();

  for(unsigned int iPat = 0; iPat < theGPs.size(); iPat++) {
    OMTFConfiguration::PatternPt patternPt;
    int charge = theGPs[iPat]->key().theCharge;
    if(theGPs[iPat] ==  0 || theGPs[iPat]->key().thePt == 0) {
      patternPts.push_back(patternPt);
      continue;
    }

    patternPt.ptFrom = OMTFConfiguration::hwPtToGev(theGPs[iPat]->key().thePt);

    unsigned int iPat1 = iPat;
    while(true) { //to skip the empty patterns with pt=0 and patterns with opposite charge
      iPat1++;
      if(iPat1 == theGPs.size())
        break;
      if(theGPs[iPat1]->key().thePt != 0 && theGPs[iPat1]->key().theCharge == charge)
        break;
    }

    if(iPat1 == theGPs.size() )
      patternPt.ptTo = 10000; //inf
    else
      patternPt.ptTo = OMTFConfiguration::hwPtToGev(theGPs[iPat1]->key().thePt );

    patternPt.charge = charge;
    patternPts.push_back(patternPt);
  }

/*  for(unsigned int iPat = 0; iPat < theGPs.size(); iPat++) {
    std::cout<<theGPs[iPat]->key()<<" ptFrom "<<patternPts[iPat].ptFrom<<" ptFrom "<<patternPts[iPat].ptTo<<std::endl;
  }*/
}

//to force compiler to compile the above methods with needed GoldenPatterns types
template class ProcessorBase<GoldenPattern>;
template class ProcessorBase<GoldenPatternParametrised>;
template class ProcessorBase<GoldenPatternWithStat>;
template class ProcessorBase<GoldenPatternPdfGen>;
template class ProcessorBase<GoldenPatternPdf4D>;
