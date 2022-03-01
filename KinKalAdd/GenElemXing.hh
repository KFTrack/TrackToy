#ifndef KinKal_GenElemXing_hh
#define KinKal_GenElemXing_hh
//
//  Describe the material effects of a kinematic trajectory crossing a straw
//  Used in the kinematic Kalman fit
//
#include "KinKal/Detector/ElementXing.hh"
//#include "KinKal/Detector/StrawMaterial.hh"
//#include "KinKal/Detector/StrawXingConfig.hh"
//#include "KinKal/Detector/WireHit.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"

namespace KinKal {
  template <class KTRAJ> class GenElemXing : public ElementXing<KTRAJ> {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using EXING = ElementXing<KTRAJ>;
      using PTCA = PiecewiseClosestApproach<KTRAJ,Line>;
      // construct from PTCA
      GenElemXing(PTCA const& tpoca, const MatEnv::DetMaterial *emat, double epath) : tpdata_(tpoca.tpData()), tprec_(tpoca.precision()), epath_(epath),
      emat_(emat),
//      sxconfig_(0.05*smat.strawRadius(),1.0),
      axis_(tpoca.sensorTraj()) {
        update(tpoca);
      }
      virtual ~GenElemXing() {}
      // ElementXing interface
      void update(PKTRAJ const& pktraj,MetaIterConfig const& miconfig) override;
      double crossingTime() const override { return tpdata_.particleToca(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // specific interface: this xing is based on PTCA
      void update(PTCA const& tpoca);
      // accessors
      ClosestApproachData const& closestApproach() const { return tpdata_; }
      double tpocaPrecision() const { return tprec_; }
//      StrawMaterial const& strawMaterial() const { return smat_; }
//      StrawXingConfig const& config() const { return sxconfig_; }
    private:
      ClosestApproachData tpdata_; // result of most recent TPOCA
      double tprec_; // precision of TPOCA
      double epath_;
//      StrawMaterial const& smat_;
      const MatEnv::DetMaterial* emat_; // material of the element
//      StrawXingConfig sxconfig_;
      Line axis_; // straw axis, expressed as a timeline

      // should add state for displace wire TODO
  };

  template <class KTRAJ> void GenElemXing<KTRAJ>::update(PTCA const& tpoca) {
    if(tpoca.usable()){
      EXING::matXings().clear();
      tpdata_ = tpoca.tpData();
      EXING::matXings().push_back(MaterialXing(*emat_,epath_));
    } else
         throw std::runtime_error("CA failure");

  }

  template <class KTRAJ> void GenElemXing<KTRAJ>::update(PKTRAJ const& pktraj,MetaIterConfig const& miconfig) {
    // search for an update to the xing configuration among this meta-iteration payload
//    auto sxconfig = miconfig.findUpdater<StrawXingConfig>();
//    if(sxconfig != 0) sxconfig_ = *sxconfig;
    // use current xing time create a hint to the CA calculation: this speeds it up
    CAHint tphint(this->crossingTime(), this->crossingTime());
    PTCA tpoca(pktraj,axis_,tphint,tprec_);
    update(tpoca);
  }

  template <class KTRAJ> void GenElemXing<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost <<"Generic Material Element Xing time " << this->crossingTime();
    if(detail > 0){
      for(auto const& mxing : this->matXings()){
        ost << " " << mxing.dmat_.name() << " pathLen " << mxing.plen_;
      }
    }
    if(detail > 1){
      ost << " Axis ";
      axis_.print(ost,0);
    }
    ost << std::endl;
  }

}
#endif
