#include <cassert>
#include <iostream>
#include <strstream>
#include <algorithm>
#include <bitset>


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFSorterWithThreshold.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternWithStat.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
template <class GoldenPatternType>
AlgoMuon OMTFSorterWithThreshold<GoldenPatternType>::sortRefHitResults(unsigned int procIndx, unsigned int iRefHit, const std::vector< std::shared_ptr<GoldenPatternType> >& gPatterns,
					  int charge){
  GoldenPatternWithThresh* bestGP = 0; //the GoldenPattern with the best result for this iRefHit
//  std::cout <<" ====== sortRefHitResults: " << std::endl;

  double p_deltaPhis1 = 0;
  //Calculating p_deltaPhis = P(delta_phis) = Sum_pt P(delta_phis | ptMu == gpPt) * P(ptMu == gpPt) - denominator of the Bayes formula
  unsigned int maxFiredLayerCnt = 0;
  for(auto& itGP: gPatterns) {
    GoldenPatternResult& result = itGP->getResults()[procIndx][iRefHit];
    if(!result.isValid())
      continue;

    if(charge!=0 && itGP->key().theCharge != charge)
      continue; //charge==0 means ignore charge

    ///Accept only candidates with >2 hits
    if(result.getFiredLayerCnt() < 3) //TODO - move 3 to the configuration??
      continue;

    //the class probability is stored in the pdf of refLayer, in the central bin, and is included in the getPdfWeigtSum()
    //result.getFiredLayerCnt() = P(delta_phis | ptMu == gpPt) * P(ptMu == gpPt)
    if(result.getFiredLayerCnt() > maxFiredLayerCnt ) {
      maxFiredLayerCnt = result.getFiredLayerCnt();
      p_deltaPhis1 = result.getPdfSum(); //cleaning p_deltaPhis acquired with the smaller number of fired layers
    }
    else if(result.getFiredLayerCnt() == maxFiredLayerCnt ) {
      p_deltaPhis1 += result.getPdfSum();
    }
  }

  double maxGpProbability1 = 0;
  //std::cout<<__FUNCTION__<<" line "<<__LINE__<<" procIndx "<<procIndx<<" iRefHit "<<iRefHit<<" charge "<<charge<<std::endl;
  double p_deltaPhis2 = 0;
  for(auto& itGP: gPatterns) {
    GoldenPatternResult& result = itGP->getResults()[procIndx][iRefHit];
    if(!result.isValid())
      continue;

    if(charge!=0 && itGP->key().theCharge != charge)
      continue; //charge==0 means ignore charge

    ///Accept only candidates with >2 hits
    if(result.getFiredLayerCnt() < 3) //TODO - move 3 to the configuration??
      continue;

    if(result.getFiredLayerCnt() < maxFiredLayerCnt)
      continue;

    //calculating P(ptMu = gpPt | delta_phis)
    if(p_deltaPhis1 > 0) {
      double gpProbability1 = result.getPdfSum() / p_deltaPhis1;
      result.setGpProbability1(gpProbability1);
    }

    //calculating P(ptMu >= gpPt | delta_phis)
    p_deltaPhis2 += result.getPdfSum();
    if(p_deltaPhis2 > 0) {
      double gpProbability2 = result.getPdfSum() / p_deltaPhis2;
      result.setGpProbability2(gpProbability2);
    }

/*

    int refLayerLogicNumber = myOmtfConfig->getRefToLogicNumber()[result.getRefLayer()];
    double ptProbability = itGP->pdfValue(refLayerLogicNumber, result.getRefLayer(), myOmtfConfig->nPdfBins()/2 );
    std::cout<<__FUNCTION__<<" line "<<__LINE__ << itGP->key()<<" refLayerLogicNumber "<<refLayerLogicNumber
        <<" FiredLayerCnt "<<result.getFiredLayerCnt()
        <<" PdfWeigtSum "<< result.getPdfWeigtSum()
        <<" p_deltaPhis "<<p_deltaPhis
        <<" ptProbability "<<ptProbability
        <<" gpProbability1 "<<gpProbability1
        << std::endl;
*/

    //assuming the patterns are sorted from the smallest pt
    /*if(result.getGpProbability1() > itGP->getTreshold(result.getRefLayer() )) {
      bestGP = itGP.get();
      //we take the one with the highest pattern number (i.e. pt) among these with the same FiredLayerCnt
      //std::cout<<__FUNCTION__<<" line "<<__LINE__ << " passed threshold "<<itGP->getTreshold(result.getRefLayer() )<<std::endl;
    }*/

    if(result.getGpProbability1() > maxGpProbability1) { //max likelihood option, TODO - remove
      maxGpProbability1 = result.getGpProbability1();
      bestGP = itGP.get();
      //std::cout<<__FUNCTION__<<" line "<<__LINE__ << " passed threshold "<<itGP->getTreshold(result.getRefLayer() )<<std::endl;
    }

  }


  if(bestGP) {
     AlgoMuon candidate(bestGP->getResults()[procIndx][iRefHit], bestGP, iRefHit);
     //std::cout<<__FUNCTION__<<" line "<<__LINE__ <<" return: " << candidate << std::endl;
     return candidate;
  }
  else {
    AlgoMuon candidate;
    candidate.setRefHitNumber(iRefHit);
    return candidate;
  }
}

template class OMTFSorterWithThreshold<GoldenPatternWithStat>;
template class OMTFSorterWithThreshold<GoldenPatternWithThresh>;
