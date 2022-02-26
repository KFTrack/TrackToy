//
// Cylindrical tracker
//
#ifndef TrackToy_Detector_Tracker_hh
#define TrackToy_Detector_Tracker_hh
#include "TrackToy/Detector/HollowCylinder.hh"
#include "TrackToy/Detector/TrckrWall.hh"
#include "TrackToy/General/TrajUtilities.hh"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "TrackToy/General/ELossDistributions.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Detector/StrawXing.hh"
#include "KinKal/Detector/StrawMaterial.hh"
#include "KinKal/Examples/SimpleWireHit.hh"
#include "TRandom3.h"
#include "Math/VectorUtil.h"
#include <string>
#include <vector>
#include <iostream>
namespace TrackToy {
  class Tracker {
    public:
      enum CellOrientation{azimuthal=0,axial};
      Tracker(MatEnv::MatDBInfo const& matdbinfo,std::string const& tfile,TRandom& tr);
      ~Tracker();
      auto const& cylinder() const { return cyl_; }
      unsigned nCells() const { return ncells_; }
      double density() const { return density_;} // gm/mm^3
      double driftVelocity() const { return vdrift_; } // mm/ns
      double propagationVelocity() const { return vprop_; } // mm/ns
      double driftResolution() const { return sigt_; } // ns
      double propagationResolution() const { return sigl_; } // ns
      double cellDensity() const { return cellDensity_;} // N/mm
      double cellRadius() const { return smat_->strawRadius(); }
      double rMin() const { return cyl_.rmin(); }
      double rMax() const { return cyl_.rmax(); }
      double zMin() const { return cyl_.zmin(); }
      double zMax() const { return cyl_.zmax(); }
      double zMid() const { return cyl_.zpos(); }
      bool InWallIsPresent() const { return inwall_!=nullptr; }
      bool EndCapWallUpStIsPresent() const { return ecpwallups_!=nullptr; }
      bool EndCapWallDwStIsPresent() const { return ecpwalldws_!=nullptr; }
      auto const& InWall() const { return *inwall_; }
      auto const& EndCapWallUpSt() const { return *ecpwallups_; }
      auto const& EndCapWallDwSt() const { return *ecpwalldws_; }
      void print(std::ostream& os) const;
      const KinKal::StrawMaterial* strawMaterial() const { return smat_; }
      // simulate hit and straw crossings along a particle trajectory.  This also updates the trajectory for BField and material effects
      template<class KTRAJ> void simulateHits(KinKal::BFieldMap const& bfield,
          KinKal::ParticleTrajectory<KTRAJ>& mctraj,
          std::vector<std::shared_ptr<KinKal::Hit<KTRAJ>>>& hits,
          std::vector<std::shared_ptr<KinKal::ElementXing<KTRAJ>>>& xings,
          std::vector<KinKal::TimeRange>& tinters,
          double tol) const;
//      template<class KTRAJ> void inWallcrossing(KinKal::BFieldMap const& bfield,
//          KinKal::ParticleTrajectory<KTRAJ>& mctraj,
//          std::vector<std::shared_ptr<KinKal::ElementXing<KTRAJ>>>& xings,
//          std::vector<KinKal::TimeRange>& tinters,
//          double tol) const;

    private:
      // helper functions for simulating hits
      // create a line representing the wire for a time on a particle trajector.  This embeds the timing information
      template <class KTRAJ> KinKal::Line wireLine(KinKal::ParticleTrajectory<KTRAJ> const& mctraj, KinKal::VEC3 wdir, double htime) const;
      // simulate hit and xing for a particular time on a particle trajectory and add them to the lists
      template <class KTRAJ> bool simulateHit(KinKal::BFieldMap const& bfield,
          KinKal::ParticleTrajectory<KTRAJ> const& mctraj,
          double htime,
          std::vector<std::shared_ptr<KinKal::Hit<KTRAJ>>>& hits,
          std::vector<std::shared_ptr<KinKal::ElementXing<KTRAJ>>>& xings ) const;
      // udpdate the trajectory for material effects
      template <class KTRAJ> void updateTraj(KinKal::BFieldMap const& bfield,
          KinKal::ParticleTrajectory<KTRAJ>& mctraj, const KinKal::ElementXing<KTRAJ>* sxing) const;
    private:
      HollowCylinder cyl_; // geometric form of the tracker
      CellOrientation orientation_; // orientation of the cells
      unsigned ncells_; // number of cells
      double cellDensity_; // effective linear cell density mm^-1
      double minpath_; // minimum path to score hits
      double density_; // total average density gm/mm^3
      KinKal::StrawMaterial* smat_; // straw material
      double vdrift_; // drift velocity
      double vprop_; // signal propagation velocity
      double sigt_; // transverse measurement time resolution sigma
      double sigl_; // longitudinal measurement time resolution sigma
      double lrdoca_; // minimum doca to resolve LR ambiguity
      double hiteff_; // hit efficiency
      double lowscat_, hiscat_, corefrac_; // factors for non-Gaussian scattering
      TRandom& tr_; // random number generator
      TrckrWall *inwall_;
      TrckrWall *ecpwallups_;
      TrckrWall *ecpwalldws_;
  };

  template<class KTRAJ> void Tracker::simulateHits(KinKal::BFieldMap const& bfield,
      KinKal::ParticleTrajectory<KTRAJ>& mctraj,
      std::vector<std::shared_ptr<KinKal::Hit<KTRAJ>>>& hits,
      std::vector<std::shared_ptr<KinKal::ElementXing<KTRAJ>>>& xings,
      std::vector<KinKal::TimeRange>& tinters, double tol) const {
    using KinKal::TimeRange;
    double tstart = mctraj.back().range().begin();
    double speed = mctraj.speed(tstart);
    double tstep = cellRadius()/speed;
    // find intersections with tracker
    TimeRange trange = cylinder().intersect(mctraj,tstart,tstep);
    while( (!trange.null()) && trange.end() < mctraj.range().end()) {
      double clen = trange.range()*speed;
      if(clen > minpath_){
        tinters.push_back(trange);
        unsigned ncells = (unsigned)floor(clen*cellDensity_);
        double hstep = trange.range()/(ncells+1);
        double htime = trange.begin()+0.5*hstep;
        for(unsigned icell=0;icell<ncells;++icell){
          // extend the trajectory to this time
          extendTraj(bfield,mctraj,htime,tol);
          // create hits and xings for this time
          bool hashit=simulateHit(bfield,mctraj,htime,hits,xings);
          // update the trajector for the effect of this material
          if(hashit){
            updateTraj(bfield, mctraj,xings.back().get());
          }
          // update to the next
          htime += hstep;
        }
      }
      trange = cylinder().intersect(mctraj,trange.end(),tstep);
    }
  }

  template <class KTRAJ> bool Tracker::simulateHit(KinKal::BFieldMap const& bfield, KinKal::ParticleTrajectory<KTRAJ> const& mctraj,
      double htime,
      std::vector<std::shared_ptr<KinKal::Hit<KTRAJ>>>& hits,
      std::vector<std::shared_ptr<KinKal::ElementXing<KTRAJ>>>& xings ) const {
    using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,KinKal::Line>;
    using WIREHIT = KinKal::SimpleWireHit<KTRAJ>;
    using STRAWXING = KinKal::StrawXing<KTRAJ>;
    using STRAWXINGPTR = std::shared_ptr<STRAWXING>;
    using ROOT::Math::VectorUtil::PerpVector;
    // define the drift and wire directions
    static const KinKal::VEC3 zdir(0.0,0.0,1.0);
    KinKal::VEC3 wdir = zdir;
    if(orientation_ == azimuthal){
      auto pos = mctraj.position3(htime);
      static double pi_over_3 = M_PI/3.0;
      // generate azimuthal angle WRT the hit
      double wphi = tr_.Uniform(0.0,pi_over_3);
      // acceptance test, using a triangular inner hole
      double rmax = 0.5*cyl_.rmax()/cos(wphi);
      if(pos.Rho() < rmax)return false;
      auto rdir = PerpVector(pos,zdir).Unit(); // radial direction
      auto fdir = zdir.Cross(rdir).Unit(); // azimuthal direction (increasing)
      wdir = fdir*cos(wphi) + rdir*sin(wphi);
      //      std::cout << "wire rmin " << pos.Rho()*cos(phimax) << std::endl;
    }
    // create the line representing this hit's wire.  The line embeds the timing information
    KinKal::Line const& wline = wireLine(mctraj,wdir,htime);
    // find the TPOCA between the particle trajectory and the wire line
    KinKal::CAHint tphint(htime,htime);
    static double tprec(1e-8); // TPOCA precision
    PTCA tp(mctraj,wline,tphint,tprec);
    // create the straw xing (regardless of inefficiency)
    auto xing = std::make_shared<STRAWXING>(tp,*smat_);
    xings.push_back(xing);
    // test for inefficiency
    double eff = tr_.Uniform(0.0,1.0);
    if(eff < hiteff_){
      // check
      //    std::cout << "doca " << tp.doca() << " sensor TOCA " << tp.sensorToca() - fabs(tp.doca())/vdrift_ << " particle TOCA " << tp.particleToca() << " hit time " << htime << std::endl;
      // define the initial ambiguity; it is the MC true value by default
      // preset to a null hit
      KinKal::WireHitState::State ambig(KinKal::WireHitState::null);
      if(fabs(tp.doca())> lrdoca_){
        ambig = tp.doca() < 0 ? KinKal::WireHitState::left : KinKal::WireHitState::right;
      }
      KinKal::WireHitState whstate(ambig);
      double mindoca = std::min(lrdoca_,cellRadius());
      // create the hit
      hits.push_back(std::make_shared<WIREHIT>(bfield, tp, whstate, mindoca, vdrift_, sigt_*sigt_, cellRadius()));
    }
    return true;
  }

  template <class KTRAJ> KinKal::Line Tracker::wireLine(KinKal::ParticleTrajectory<KTRAJ> const& mctraj, KinKal::VEC3 wdir, double htime) const {
    // find the position and direction of the particle at this time
    auto pos = mctraj.position3(htime);
    auto pdir = mctraj.direction(htime);
    // drift direction is perp to wire and particle
    auto ddir = pdir.Cross(wdir).Unit();
    // uniform drift distance = uniform impact parameter (random sign)
    double rdrift = tr_.Uniform(-cellRadius(),cellRadius());
    auto wpos = pos + rdrift*ddir;
    // find the wire ends; this is where the wire crosses the outer envelope
    double dprop, wlen;
    if(orientation_ == azimuthal){
      double rwire = wpos.Rho(); // radius of the wire position
      // find crossing of outer cylinder
      double rmax = std::max(rwire,rMax());
      double wdot = -wpos.Dot(wdir);
      double term = sqrt(wdot*wdot + (rmax*rmax - rwire*rwire));
      double d1 = wdot+term;
      double d2 = wdot-term;
      wlen = fabs(d1)+fabs(d2);
      // choose the shortest propagation distance to define the measurement (earliest time)
      if(fabs(d1) < fabs(d2)){
        dprop = fabs(d1);
      } else {
        wdir *= -1.0; // points from source to measurement location
        dprop = fabs(d2);
      }
    } else {
      // wire ends are at zmin and zmax
      double zmax = zMax() - wpos.Z();
      double zmin = zMin() - wpos.Z();
      wlen = zmax - zmin;
      if(fabs(zmax) < fabs(zmin)){
        dprop = fabs(zmax);
      } else {
        wdir *= -1.0; // points from source to measurement location
        dprop = fabs(zmin);
      }
    }
    // measurement time includes propagation and drift
    double mtime = htime + dprop/vprop_ + fabs(rdrift)/vdrift_;
    // smear measurement time by the resolution
    mtime = tr_.Gaus(mtime,sigt_);
    // construct the trajectory for this hit.  Note this embeds the timing and location information together
    auto mpos = wpos + dprop*wdir;
//    std::cout << "measurement radius " << mpos.Rho() << std::endl;
    return KinKal::Line(mpos,mtime,wdir*vprop_,wlen);
  }

  template <class KTRAJ> void Tracker::updateTraj(KinKal::BFieldMap const& bfield,
      KinKal::ParticleTrajectory<KTRAJ>& mctraj, const KinKal::ElementXing<KTRAJ>* sxing) const {
    static unsigned moyalterms_(20); // number of terms in Moyal expansion
    // simulate energy loss and multiple scattering from this xing
    auto txing = sxing->crossingTime();
    auto const& endpiece = mctraj.nearestPiece(txing);
    auto endmom = endpiece.momentum4(txing);
    KinKal::VEC3 momvec = endmom.Vect();
    double mom = momvec.R();
    auto endpos = endpiece.position4(txing);
    std::array<double,3> dmom {0.0,0.0,0.0}, momvar {0.0,0.0,0.0};
    sxing->materialEffects(mctraj,KinKal::TimeDir::forwards, dmom, momvar);
//   std::cout << "DMOM " << dmom[0] << std::endl;
   // note materialEffects returns the normalized energy change
    MoyalDist edist(MoyalDist::MeanRMS(fabs(dmom[KinKal::MomBasis::momdir_])*mom,sqrt(momvar[KinKal::MomBasis::momdir_])),moyalterms_);
    // radiation energy loss model
    BremssLoss bLoss;
    for(int idir=0;idir<=KinKal::MomBasis::phidir_; idir++) {
      auto mdir = static_cast<KinKal::MomBasis::Direction>(idir);
      double momsig = sqrt(momvar[idir]);
      double ionloss, bremloss, dloss;
      double dm;
      double radFrac = sxing->radiationFraction()/10; // convert to cm
      // only include wall material for now; gas can be added later TODO
      DeltaRayLoss dLoss(&(sxing->matXings()[0].dmat_), mom,sxing->matXings()[0].plen_/10.0, endpiece.mass());
      // generate a random effect given this variance and mean.  Note momEffect is scaled to momentum
      switch( mdir ) {
        case KinKal::MomBasis::perpdir_: case KinKal::MomBasis::phidir_ :
          // non-Gaussian scattering tails
          dm = tr_.Uniform(0.0,1.0) < corefrac_ ?  tr_.Gaus(dmom[idir],lowscat_*momsig) : tr_.Gaus(dmom[idir],hiscat_*momsig);
          dm *= mom; // scale to momentum
          break;
        case KinKal::MomBasis::momdir_ :
          // calculate a smeared energy loss using a Moyal distribution
          ionloss = edist.sample(tr_.Uniform(0.0,1.0));
          bremloss = bLoss.sampleSSPGamma(mom,radFrac);
          dloss = dLoss.sampleDRL();
//          std::cout << "Tracker Ionization eloss = " << ionloss << " Delta eloss " << dloss << " rad eloss "  << bremloss << std::endl;
          dm = std::max(-(ionloss+bremloss+dloss),-mom); // must be positive
          break;
        default:
          throw std::invalid_argument("Invalid direction");
      }
      momvec += endpiece.direction(txing,mdir)*dm;
    }
//    std::cout << "startmom " << mom << " endmom " << momvec.R() << std::endl;
    endmom.SetCoordinates(momvec.X(), momvec.Y(), momvec.Z(),endmom.M());
    // generate a new piece and append
    auto bnom = bfield.fieldVect(endpos.Vect());
//    auto bnom = endpiece.bnom();
    KTRAJ newend(endpos,endmom,endpiece.charge(),bnom,KinKal::TimeRange(txing,mctraj.range().end()));
//    mctraj.append(newend);
    mctraj.append(newend,true); // allow truncation if needed
  }

//  template<class KTRAJ> void Tracker::inWallcrossing(KinKal::BFieldMap const& bfield,
//      KinKal::ParticleTrajectory<KTRAJ>& mctraj,
//      std::vector<std::shared_ptr<KinKal::ElementXing<KTRAJ>>>& xings,
//      std::vector<KinKal::TimeRange>& tinters, double tol) const {
//    if (inwall_!=nullptr) {
//      double tstart = mctraj.back().range().begin();
//      double speed = mctraj.speed(tstart);
//      double tstep = inwall_->cylinder().rhalf()/speed;
//      // find intersections with tracker
//      inwall_->cylinder().intersect(mctraj,tinters,tstart,tstep);
//      //    std::cout << "ninters " << tinters.size() << std::endl;
////      for(auto const& tinter : tinters) {
////        double clen = tinter.range()*speed;
////        //
////        if(clen > minpath_){
////          unsigned ncells = (unsigned)floor(clen*cellDensity_);
////          double hstep = tinter.range()/(ncells+1);
////          double htime = tinter.begin()+0.5*hstep;
////          for(unsigned icell=0;icell<ncells;++icell){
////            // extend the trajectory to this time
////            extendTraj(bfield,mctraj,htime,tol);
////            // create hits and xings for this time
////            bool hashit=simulateHit(bfield,mctraj,htime,hits,xings);
////            // update the trajector for the effect of this material
////            if(hashit){
////              updateTraj(bfield, mctraj,xings.back().get());
////            }
////            // update to the next
////            htime += hstep;
////          }
////        }
////      }
//    }
//  }

}
#endif

