//
// Detector element modeled as a 2-dimensional cylindrial shell
// Structured text file constructor should have contain: radius, rhalf, zpos, halfzlength
//
#ifndef TrackToy_Detector_DiscShell_hh
#define TrackToy_Detector_DiscShell_hh
#include "TrackToy/Detector/CylindricalShell.hh"

namespace TrackToy {
  class DiscShell : public CylindricalShell {
    public:
      DiscShell(): rmin_(-1.0), rmax_(-1.0) {}
      DiscShell(double radius, double rhalf, double zpos, double zhalf) : CylindricalShell( radius, rhalf, zpos, zhalf ) {
        rmin_=radius-rhalf;
        rmax_=radius+rhalf;
      }
      double rmin() const { return rmin_; }
      double rmax() const { return rmax_; }
      // find the 1st intersection of the trajectory with this cylinder, starting from the given time
      template<class PKTRAJ> KinKal::TimeRange intersect(PKTRAJ const& pktraj, double tstart, double tstep) const;
    private:
      double rmin_, rmax_;
  };

  template<class PKTRAJ> KinKal::TimeRange DiscShell::intersect(PKTRAJ const& pktraj, double tstart, double tstep) const {
    using KinKal::TimeRange;
    double ttest = tstart;
    TimeRange trange(ttest,ttest);

//scan version to search for intersection
//    auto pos = pktraj.position3(ttest);
//    double tin=-1, tout=-1;
//    double dz=0.0, dr = 0.0;
//    int nstep=0;
//    while(ttest < pktraj.range().end() && nstep<100){
//      auto vel = pktraj.velocity(ttest);
//      dz = fabs(tstep*vel.Z());
//      dr = fabs(tstep*vel.Rho());
//
//      bool inZ = pos.Z()>zmin()-dz && pos.Z()<zmax()+dz;
//      bool inR = pos.Rho()>rmin_-dr && pos.Rho()<rmax_+dr;
//      if ( inZ && inR ) {
//        if ( tin<0 ) {
//          tin=ttest;
//          double tmpStep = 2.0*zhalf()/vel.Z();
//          if (tstep>tmpStep) { tstep=tmpStep*0.1; }
//        } else {
//          tout=ttest;
//        }
//      } else if (tin>0 && tout>0 ) {
//        trange = TimeRange(tin,tout);
//        break;
//      }
//      ttest+= tstep;
//      pos = pktraj.position3(ttest);
//      ++nstep;
//
//    }
    auto pos = pktraj.position3(ttest);
    double tin=-1;
    double dz=0.0, dr = 0.0;
    bool inR=false;
    int niter=0;
    auto vel = pktraj.velocity(ttest);
    if (vel.Z()>0) {
      dz = tstep*vel.Z();
      niter=0;
      while( pos.Z()<(zmin()-dz) && ttest < pktraj.range().end() && niter<1000 ) {
        ttest+=tstep;
        pos = pktraj.position3(ttest);
        ++niter;
      }
      if (ttest < pktraj.range().end() && niter<1000) {
        vel = pktraj.velocity(ttest);
        tstep = 2.0*zhalf()/vel.Z()*0.05;
        dz = tstep*vel.Z();
        dr = tstep*vel.Rho();
        niter=0;
        while( pos.Z()<zmax() && ttest < pktraj.range().end() && niter<1000 ) {
          inR = pos.Rho()>rmin_-dr && pos.Rho()<rmax_+dr;
          if (tin<0 && pos.Z()>(zmin()-dz) && inR ) {
            tin=ttest;
          }
          if (tin>0 && !inR ) { break; }
          ttest+=tstep;
          pos = pktraj.position3(ttest);
          ++niter;
        }
      }
    } else {
      dz = -tstep*vel.Z();
      niter=0;
      while( pos.Z()>(zmax()+dz) && ttest < pktraj.range().end() && niter<1000 ) {
        ttest+=tstep;
        pos = pktraj.position3(ttest);
        ++niter;
      }
      if (ttest < pktraj.range().end() && niter<1000) {
        vel = pktraj.velocity(ttest);
        tstep = -2.0*zhalf()/vel.Z()*0.05;
        dz = tstep*vel.Z();
        dr = tstep*vel.Rho();
        niter=0;
        while( pos.Z()>zmin() && niter<1000 ) {
          inR = pos.Rho()>rmin_-dr && pos.Rho()<rmax_+dr;
          if (tin<0 && pos.Z()<(zmax()+dz) && inR ) {
            tin=ttest;
          }
          if (tin>0 && !inR ) { break; }
          ttest+=tstep;
          pos = pktraj.position3(ttest);
          ++niter;
        }
      }
    }
    if ( tin>0 && ttest < pktraj.range().end() && niter<1000 ) { trange = TimeRange(tin,ttest); }

//end of scan version

////analytic version of the intersection evaluation, note, currently it assume only track that are traveling in the forward direction
//    auto vel = pktraj.velocity(ttest);
//    double dt = 2.0*zhalf()/vel.Z();
//    double texit = ttest+dt; //ztime(pktraj,pktraj.back().range().begin(),zmax());
//    auto pstart = pktraj.position3(ttest);
//    auto pexit = pktraj.position3(texit);
////    bool reducedInt=false;
//    bool inRstart = pstart.Rho()>rmin_ && pstart.Rho()<rmax_;
//    bool inRexit = pexit.Rho()>rmin_ && pexit.Rho()<rmax_;
//
//    if ( inRstart && inRexit ) {
//      trange = TimeRange(ttest,texit);
//    } else if ( inRstart ) {
//      //      double dt=texit-ttest;
//      if (pexit.Rho()<rmin_) {
////        reducedInt=true;
//        dt*=(pstart.Rho()-rmin_)/(pstart.Rho()-pexit.Rho());
//        texit = ttest+dt;
//        trange = TimeRange(ttest,texit);
////        pexit = pktraj.position3(texit);
//      } else if ( pexit.Rho()>rmax_) {
////        reducedInt=true;
//        dt*=(rmax_-pstart.Rho())/(pexit.Rho()-pstart.Rho());
//        texit = ttest+dt;
//        trange = TimeRange(ttest,texit);
////        pexit = pktraj.position3(texit);
//      }
//    } else if ( inRexit ) {
//      //      double dt=texit-ttest;
//      if (pstart.Rho()<rmin_) {
////        reducedInt=true;
//        dt*=(pexit.Rho()-rmin_)/(pexit.Rho()-pstart.Rho());
//        ttest = texit-dt;
//        trange = TimeRange(ttest,texit);
////        pstart = pktraj.position3(ttest);
//      } else if ( pstart.Rho()>rmax_) {
////        reducedInt=true;
//        dt*=(rmax_-pexit.Rho())/(pstart.Rho()-pexit.Rho());
//        ttest = texit-dt;
//        trange = TimeRange(ttest,texit);
////        pstart = pktraj.position3(ttest);
//      }
//    }
////end of analytic version

    return trange;
  }

}
#endif

