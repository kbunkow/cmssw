#ifndef PHASE2GMT_PRETRACKMATCHEDMUON
#define PHASE2GMT_PRETRACKMATCHEDMUON

#include <L1Trigger/Phase2L1GMT/interface/ApSignAbsInt.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"


namespace Phase2L1GMT {

  class PreTrackMatchedMuon {
  public:
    PreTrackMatchedMuon(const int& curvature,
                        const uint& charge,
                        const uint& pt,
                        const int& eta,
                        const int& phi,
                        const int& z0,
                        const int& d0,
                        const uint& beta = 15)
        : curvature_(curvature),
          charge_(charge),
          pt_(pt),
          eta_(eta),
          phi_(phi),
          z0_(z0),
          d0_(d0),
          beta_(beta),
          isGlobal_(false),
          quality_(0),
          stubID0_(511),
          stubID1_(511),
          stubID2_(511),
          stubID3_(511),
          stubID4_(511),
          valid_(false),
          deltaCoords1(5),
          deltaCoords2(5),
          deltaEtas1(5),
          deltaEtas2(5)
    {}

    const uint curvature() const { return curvature_; }
    const uint charge() const { return charge_; }
    const uint pt() const { return pt_; }
    const int eta() const { return eta_; }
    const int phi() const { return phi_; }
    const int z0() const { return z0_; }
    const int d0() const { return d0_; }
    const uint beta() const { return beta_; }

    bool isGlobalMuon() const { return isGlobal_; }
    const int quality() const { return quality_; }
    const int offline_pt() const { return offline_pt_; }
    const float offline_eta() const { return offline_eta_; }
    const float offline_phi() const { return offline_phi_; }

    const uint stubID0() const { return stubID0_; }
    const uint stubID1() const { return stubID1_; }
    const uint stubID2() const { return stubID2_; }
    const uint stubID3() const { return stubID3_; }
    const uint stubID4() const { return stubID4_; }

    const uint stubID(unsigned int layer) const {
      if(layer == 0)
        return stubID0_;
      if(layer == 1)
        return stubID1_;
      if(layer == 2)
        return stubID2_;
      if(layer == 3)
        return stubID3_;
      if(layer == 4)
        return stubID4_;

      return 0;
    }

    bool valid() const { return valid_; }

    void setQuality(uint quality) { quality_ = quality; }
    void setValid(bool v) { valid_ = v; }

    void setOfflineQuantities(float pt, float eta, float phi) {
      offline_pt_ = pt;
      offline_eta_ = eta;
      offline_phi_ = phi;
    }

    void addMuonRef(const l1t::RegionalMuonCandRef& ref) {
      muRef_.push_back(ref);
      isGlobal_ = true;
    }

    void resetGlobal() { isGlobal_ = false; }

    const std::vector<l1t::RegionalMuonCandRef>& muonRef() const { return muRef_; }
    void addStub(const l1t::MuonStubRef& stub) {
      stubs_.push_back(stub);
      if (stub->tfLayer() == 0)
        stubID0_ = stub->id();
      if (stub->tfLayer() == 1)
        stubID1_ = stub->id();
      if (stub->tfLayer() == 2)
        stubID2_ = stub->id();
      if (stub->tfLayer() == 3)
        stubID3_ = stub->id();
      if (stub->tfLayer() == 4)
        stubID4_ = stub->id();
    }

    const l1t::MuonStubRefVector& stubs() const { return stubs_; }

    void setTrkPtr(const edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> >& trkPtr) { trkPtr_ = trkPtr; }

    const edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> > trkPtr() const { return trkPtr_; }

    auto& getDeltaCoords1() {
      return deltaCoords1;
    }

    auto& getDeltaCoords2() {
      return deltaCoords2;
    }

    std::vector<ApSignAbsInt<BITSSIGMAETA> >& getDeltaEtas1() {
      return deltaEtas1;
    }

    std::vector<ApSignAbsInt<BITSSIGMAETA> >& getDeltaEtas2() {
      return deltaEtas2;
    }

    //each bit denotes one phi coordinate in the combined stub, i.e. bit 0 is deltaCoords1[layer=0], bit 1 deltaCoords2[layer=0], etc.
    unsigned short getHitsValid() const {
      return hitsValid;
    }

    void setHitsValid(unsigned int coordIdx, bool valid) {
      if(valid)
        hitsValid = hitsValid | (1 << coordIdx);
      else
        hitsValid = hitsValid & ~(1 << coordIdx);
    }

    friend std::ostream& operator<<(std::ostream& out, const PreTrackMatchedMuon& mu) {
      out<< "reconstructed muon : charge=" << mu.charge_ << " pt=" << mu.offline_pt_ << ","
                                           << mu.pt_ << " eta=" << mu.offline_eta_ << "," << mu.eta_ << " phi=" << mu.offline_phi_ << ","
                                           << mu.phi_ << " z0=" << mu.z0_ << " d0=" << mu.d0_ << " quality=" << mu.quality_
                                           << " isGlobal=" << mu.isGlobal_ << " valid=" << mu.valid_ << " stubs: " << mu.stubID0_
                                           << " " << mu.stubID1_ << " " << mu.stubID2_ << " " << mu.stubID3_ << " " << mu.stubID4_;

      return out;
    }

    uint64_t lsb() const {
      uint64_t w = charge_ & 0x1;
      w = w | (twos_complement(pt_, BITSPT) << 1);
      w = w | (twos_complement(phi_, BITSPHI) << (BITSPT + 1));
      w = w | (twos_complement(eta_, BITSETA) << (BITSPHI + BITSPT + 1));
      w = w | (twos_complement(z0_, BITSZ0) << (BITSETA + BITSPHI + BITSPT + 1));
      w = w | (twos_complement(d0_, BITSD0) << (BITSZ0 + BITSETA + BITSPHI + BITSPT + 1));
      return w;
    }

    uint64_t msb() const {
      uint64_t w2 = 0;
      w2 = twos_complement(stubID0_, BITSSTUBID);
      w2 = w2 | (twos_complement(stubID1_, BITSSTUBID) << BITSSTUBID);
      w2 = w2 | (twos_complement(stubID2_, BITSSTUBID) << (2 * BITSSTUBID));
      w2 = w2 | (twos_complement(stubID3_, BITSSTUBID) << (3 * BITSSTUBID));
      w2 = w2 | (twos_complement(stubID4_, BITSSTUBID) << (4 * BITSSTUBID));
      w2 = w2 | (twos_complement(isGlobal_, 1) << (5 * BITSSTUBID));
      w2 = w2 | (twos_complement(beta_, BITSMUONBETA) << (5 * BITSSTUBID + 1));
      w2 = w2 | (twos_complement(quality_, BITSMATCHQUALITY) << (BITSMUONBETA + 5 * BITSSTUBID + 1));
      w2 = w2 | (twos_complement(valid_, 1) << (BITSMATCHQUALITY + BITSMUONBETA + 5 * BITSSTUBID + 1));
      return w2;
    }

    void printWord() const {
      LogTrace("PreTrackMatchedMuon") << "PreTrackMatchedMuon : word=" << std::setfill('0') << std::setw(16) << std::hex
                                      << (long long unsigned int)(msb() >> 2) << std::setfill('0') << std::setw(16)
                                      << std::hex
                                      << (long long unsigned int)((lsb() | (msb() << 62)) & 0xffffffffffffffff);
    }

  private:
    int curvature_ = 0;
    uint charge_;
    uint pt_;
    int eta_;
    int phi_;
    int z0_;
    int d0_;
    uint beta_;
    bool isGlobal_;
    uint quality_;
    float offline_pt_ = 0;
    float offline_eta_ = 0;
    float offline_phi_ = 0;
    uint stubID0_;
    uint stubID1_;
    uint stubID2_;
    uint stubID3_;
    uint stubID4_;
    bool valid_;
    l1t::MuonStubRefVector stubs_;
    std::vector<l1t::RegionalMuonCandRef> muRef_;
    edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> > trkPtr_;

    //index is layer, 0 - to 4
    std::vector<ApSignAbsInt<BITSSIGMACOORD> > deltaCoords1;
    std::vector<ApSignAbsInt<BITSSIGMACOORD> > deltaCoords2;

    std::vector<ApSignAbsInt<BITSSIGMAETA> > deltaEtas1;
    std::vector<ApSignAbsInt<BITSSIGMAETA> > deltaEtas2;

    unsigned short hitsValid = 0;
  };
}  // namespace Phase2L1GMT

#endif
