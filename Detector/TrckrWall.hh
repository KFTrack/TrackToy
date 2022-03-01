//
//  Inner proton absorber.  Currently a thin cylinder
//
#ifndef TrackToy_Detector_TrckrWall_hh
#define TrackToy_Detector_TrckrWall_hh
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "TrackToy/General/ELossDistributions.hh"
#include "KinKal/General/TimeRange.hh"
#include "TrackToy/General/TrajUtilities.hh"
//#include "TrackToy/Detector/CylindricalShell.hh"
#include "TrackToy/Detector/DiscShell.hh"
#include "TRandom3.h"
namespace TrackToy {
  class TrckrWall {
    public:
      enum TrckrWallType { unknown=-1, cylindrical=1, disc=2 };
      TrckrWall(MatEnv::MatDBInfo const& matdbinfo,std::string const& file,TRandom& tr);
      ~TrckrWall() { if (cyl_!=nullptr) { delete cyl_; } }
      auto const& cylinder() const { return *cyl_; }
      auto const& material() const { return *mat_; }
      auto const* materialPtr() const { return mat_; }
      auto type() const { return type_; }
      // extend a trajectory through the TrckrWall.  Return value specifies if the particle continues downsream (true) or stops in the target or exits the field upstream (false)
      template<class PKTRAJ> bool extendTrajectory(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj,TimeRanges& intersections, double tol=1e-4) const;
      template<class PKTRAJ> bool intersect(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj, KinKal::TimeRange& intersection, double tol=1e-4) const;
      template<class PKTRAJ> void intersects(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj, TimeRanges& intersections, double tol=1e-4) const;
      void print(std::ostream& os) const;
    private:
      TrckrWallType type_;
      CylindricalShell *cyl_;
      const MatEnv::DetMaterial* mat_;
      TRandom& tr_; // random number generator
  };

//  template<class PKTRAJ> bool TrckrWall::extendTrajectoryCyl(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj, TimeRanges& intersections,double tol) const {
  template<class PKTRAJ> bool TrckrWall::extendTrajectory(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj, TimeRanges& intersections, double tol) const {
	    using KinKal::VEC3;
	    using KinKal::TimeRange;
	    static unsigned moyalterms_(20); // number of terms in Moyal expansion
	    intersections.clear();
	    // record the end of the previous extension; this is where new extensions start
	    double tstart = pktraj.back().range().begin();
	    // extend through the TrckrWall or exiting the BField (backwards)
	    bool retval = extendZ(pktraj,bfield, cyl_->zmax(), tol);
	    if(retval){
	//      auto pstart = pktraj.position3(tstart);
	//      std::cout << "zstart " << pstart.Z() << " Z extend " << zipa << std::endl;
	      // first find the intersections.
	      static double tstep(0.01);
//	      TimeRange trange = cyl_->intersect(pktraj,tstart,tstep);
        TimeRange trange(0.0,0.0);
        if (type_==disc) {
          trange = ((DiscShell*)cyl_)->intersect(pktraj,tstart,tstep);
        } else {
          trange = cyl_->intersect(pktraj,tstart,tstep);
        }
	      /*while*/ if( (!trange.null()) && trange.end() < pktraj.range().end()) {
	        intersections.push_back(trange);
	        double energy = pktraj.energy(trange.mid());
	        auto momvec = pktraj.momentum3(trange.mid());
	        double mom = momvec.R();
	        double plen = pktraj.speed(trange.mid())*trange.range();
	        // Moyal dist. models ionization loss
	        double demean = mat_->energyLoss(mom,plen,pktraj.mass());
	        double derms = mat_->energyLossRMS(mom,plen,pktraj.mass());
	        // model ionization energy loss using a Moyal distribution
	//        std::cout << "TrckrWall demean " << demean << " derms " << derms << std::endl;
	        MoyalDist edist(MoyalDist::MeanRMS(fabs(demean), derms),moyalterms_);
	        double ionloss = edist.sample(tr_.Uniform(0.0,1.0));
	        // add radiative energy loss.  note we have to convert to cm!!!
	        double radFrac = mat_->radiationFraction(trange.range())/10;
	        BremssLoss bLoss;
	        double bremloss = bLoss.sampleSSPGamma(energy,radFrac);
	        // delta energy loss; note unit change mm->cm!
	        DeltaRayLoss dLoss(mat_, energy,plen/10, pktraj.mass());
	        double dloss = dLoss.sampleDRL();
	        double totloss = ionloss + bremloss + dloss;
	        //          std::cout << "TrckrWall Ionization eloss = " << ionloss << " Delta eloss " << dloss << " rad eloss "  << bremloss << " tot " << totloss << std::endl;
	        //          double oldenergy = energy;
	        energy -= totloss;
	        //          std::cout << "old energy " << oldenergy << " new energy " << energy << std::endl;
	        // scattering
	        double scatterRMS = mat_->scatterAngleRMS(mom,plen,pktraj.mass());
	        // generate random momentum scatter
	        VEC3 phidir = VEC3(momvec.Y(),-momvec.X(),0.0).Unit();
	        double momtan = tan(momvec.Theta());
	        VEC3 thedir = VEC3(momvec.X()/momtan,momvec.Y()/momtan,-momvec.Z()*momtan).Unit();
	        auto dmom = tr_.Gaus(0.0,scatterRMS)*mom*phidir + tr_.Gaus(0.0,scatterRMS)*mom*thedir;
	        retval = updateEnergy(pktraj,trange.mid(),energy,dmom);
//	        if(!retval)break;
//	        trange = cyl_->intersect(pktraj,trange.end(),tstep);
	      }
	    }
	    return retval;
	  }

  template<class PKTRAJ> bool TrckrWall::intersect(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj, KinKal::TimeRange& intersection, double tol) const {

    // record the end of the previous extension; this is where new extensions start
    double tstart = intersection.end();//ztime(pktraj,pktraj.back().range().begin(),cyl_->zmin());
    // extend through the TrckrWall or exiting the BField (backwards)
    bool retval = extendZ(pktraj,bfield, cyl_->zmax(), tol);
    if(retval){
        static double tstep(0.01);
        if (type_==disc) {
          intersection = ((DiscShell*)cyl_)->intersect(pktraj,tstart,tstep);
        } else {
          intersection = cyl_->intersect(pktraj,tstart,tstep);
        }
//        std::cout<<"trange: "<<intersection.begin()<<" "<<intersection.end()<<std::endl;
      if ( (!intersection.null()) && intersection.end() < pktraj.range().end()) {
        retval = true;
      } else {
        retval = false;
      }
    }
    return retval;
  }

  template<class PKTRAJ> void TrckrWall::intersects(KinKal::BFieldMap const& bfield, PKTRAJ& pktraj, TimeRanges& intersections, double tol) const {
    intersections.clear();
    KinKal::TimeRange intersection(0.0,0.0);
    while (intersect(bfield, pktraj, intersection, tol)) {
      intersections.push_back(intersection);
    }
  }

}
#endif

