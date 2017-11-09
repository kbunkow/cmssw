/*
 * ProcessorBase.h
 *
 *  Created on: Jul 28, 2017
 *      Author: kbunkow
 */

#ifndef OMTF_PROCESSORBASE_H_
#define OMTF_PROCESSORBASE_H_

#include <L1Trigger/L1TMuonOverlap/interface/GoldenPatternBase.h>
#include <memory>

#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"

class L1TMuonOverlapParams;
class SimTrack;

template <class GoldenPatternType>
class ProcessorBase {
public:
  typedef std::vector< std::shared_ptr<GoldenPatternType> > GoldenPatternVec;

  ProcessorBase(const OMTFConfiguration* omtfConfig): myOmtfConfig(omtfConfig)  {
  };

  virtual ~ProcessorBase() {
  }

  ///Return vector of GoldenPatterns
  virtual const GoldenPatternVec& getPatterns() const  {
    return theGPs;
  };

  ///Fill GP vec with patterns from CondFormats object
  virtual bool configure(const OMTFConfiguration * omtfParams, const L1TMuonOverlapParams* omtfPatterns);

  ///Add GoldenPattern to pattern vec.
  virtual void addGP(GoldenPatternType *aGP);

  virtual void setGPs(const GoldenPatternVec& gps) {
    theGPs = gps;
    for(auto& gp : theGPs) {
      gp->setConfig(myOmtfConfig);
    }

    initPatternPtRange();
  }

  ///Reset all configuration parameters
  virtual void resetConfiguration();

  virtual void initPatternPtRange();

  const std::vector<OMTFConfiguration::PatternPt>& getPatternPtRange() const {
    return patternPts;
  }

protected:
  ///vector holding Golden Patterns
  GoldenPatternVec theGPs;

  ///Configuration of the algorithm. This object
  ///does not contain the patterns data.
  const OMTFConfiguration* myOmtfConfig;

  ///Remove hits whis are outside input range
  ///for given processor and cone
  virtual OMTFinput::vector1D restrictInput(unsigned int iProcessor,
            unsigned int iCone,
            unsigned int iLayer,
            const OMTFinput::vector1D & layerHits);

  std::vector<OMTFConfiguration::PatternPt> patternPts;
};

#endif /* OMTF_PROCESSORBASE_H_ */
