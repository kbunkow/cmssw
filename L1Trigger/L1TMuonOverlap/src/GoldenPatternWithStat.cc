#include <iostream>
#include <iomanip>
#include <cmath>

#include <TF1.h>
#include <TLinearFitter.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternWithStat.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"


////////////////////////////////////////////////////
////////////////////////////////////////////////////
/*GoldenPatternWithStat::GoldenPatternWithStat(const Key& aKey, const OMTFConfiguration * omtfConfig): GoldenPattern(aKey, omtfConfig) {
  //GoldenPattern::reset();
  //reset();
}*/


void GoldenPatternWithStat::updateStat(unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, unsigned int what, double value) {
  statisitics[iLayer][iRefLayer][iBin][what] += value;
  //std::cout<<__FUNCTION__<<":"<<__LINE__<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" iBin "<<iBin<<" what "<<what<<" value "<<value<<std::endl;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
/*void GoldenPatternWithStat::updatePdfs(double learingRate) {
  //double f = 1;
  for(unsigned int iLayer = 0; iLayer < getPdf().size(); ++iLayer) {
    for(unsigned int iRefLayer=0; iRefLayer < getPdf()[iLayer].size(); ++iRefLayer) {
      for(unsigned int iPdf = 1; iPdf < getPdf()[iLayer][iRefLayer].size(); iPdf++) {
        double d = 0;
        if(statisitics[iLayer][iRefLayer][iPdf][whatSimNorm] != 0)
          d -= statisitics[iLayer][iRefLayer][iPdf][whatSimVal]/(double)statisitics[iLayer][iRefLayer][iPdf][whatSimNorm];

        if(statisitics[iLayer][iRefLayer][iPdf][whatOmtfNorm] != 0)
          d += statisitics[iLayer][iRefLayer][iPdf][whatOmtfVal]/(double)statisitics[iLayer][iRefLayer][iPdf][whatOmtfNorm] ;

        d = d * learingRate;
        pdfAllRef[iLayer][iRefLayer][iPdf] += d;
        if(d != 0) {
          std::cout<<__FUNCTION__<<":"<<__LINE__<<" "<< key()<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" iBin "<<iPdf<<" pdfVal "<<pdfAllRef[iLayer][iRefLayer][iPdf]<<" d "<<d<<std::endl;
        }
      }
    }
  }
}*/

std::ostream & operator << (std::ostream &out, const GoldenPatternWithStat & aPattern){

/*  out <<"GoldenPatternWithStat "<< aPattern.theKey <<std::endl;
  out <<"Number of reference layers: "<<aPattern.meanDistPhi[0].size()
          <<", number of measurement layers: "<<aPattern.pdfAllRef.size()
          <<std::endl;

  if(!aPattern.meanDistPhi.size()) return out;
  if(!aPattern.pdfAllRef.size()) return out;

  out<<"Mean dist phi per layer:"<<std::endl;
  for (unsigned int iRefLayer=0;iRefLayer<aPattern.meanDistPhi[0].size();++iRefLayer){
    out<<"Ref layer: "<<iRefLayer<<" (";
    for (unsigned int iLayer=0;iLayer<aPattern.meanDistPhi.size();++iLayer){   
      out<<std::setw(3)<<aPattern.meanDistPhi[iLayer][iRefLayer]<<"\t";
    }
    out<<")"<<std::endl;
  }

  if(aPattern.meanDistPhiCounts.size()){
    out<<"Counts number per layer:"<<std::endl;
    for (unsigned int iRefLayer=0;iRefLayer<aPattern.meanDistPhi[0].size();++iRefLayer){
      out<<"Ref layer: "<<iRefLayer<<" (";
      for (unsigned int iLayer=0;iLayer<aPattern.meanDistPhi.size();++iLayer){   
        out<<aPattern.meanDistPhiCounts[iLayer][iRefLayer]<<"\t";
      }
      out<<")"<<std::endl;
    }
  }

  unsigned int nPdfAddrBits = 7;
  out<<"PDF per layer:"<<std::endl;
  for (unsigned int iRefLayer=0;iRefLayer<aPattern.pdfAllRef[0].size();++iRefLayer){
    out<<"Ref layer: "<<iRefLayer;
    for (unsigned int iLayer=0;iLayer<aPattern.pdfAllRef.size();++iLayer){   
      out<<", measurement layer: "<<iLayer<<std::endl;
      for (unsigned int iRefPhiB=0; iRefPhiB < aPattern.pdfAllRef[iLayer][iRefLayer].size(); ++iRefPhiB) {
        for(unsigned int iPdf=0;iPdf<exp2(nPdfAddrBits);++iPdf){
          out<<std::setw(2)<<aPattern.pdfAllRef[iLayer][iRefLayer][iRefPhiB][iPdf]<<" ";
        }
      }
      out<<std::endl;
    }
  }*/

  return out;
}

