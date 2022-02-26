//
//  Implementation of tracker
//
#include "TrackToy/Detector/Tracker.hh"
#include "TrackToy/General/FileFinder.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <stdexcept>

namespace TrackToy {
  Tracker::Tracker(MatEnv::MatDBInfo const& matdbinfo,std::string const& trkfile,TRandom& tr):
    ncells_(0), cellDensity_(-1.0), density_(-1.0), smat_(0),
    vdrift_(-1.0), vprop_(-1.0), sigt_(-1.0), sigl_(-1.0), lrdoca_(-1.0), hiteff_(-1.0), hiscat_(10.0), tr_(tr)
	, inwall_(nullptr), ecpwallups_(nullptr), ecpwalldws_(nullptr)
  {
    FileFinder filefinder;
    std::string fullfile = filefinder.fullFile(trkfile);
    std::string line;
    static std::string comment("#");
    std::ifstream tracker_stream(fullfile,std::ios_base::in);
    if(tracker_stream.fail()){
      std::string errmsg = std::string("File doesn't exist" )+ fullfile;
      throw std::invalid_argument(errmsg.c_str());
    }
    double rmin(-1.0), rmax(-1.0), zpos(-1.0), zhalf(-1.0);
    double rcell(-1.0), lcell(-1.0), rwire(-1.0), wthick(-1.0);
    unsigned orient(0);
    std::string strawwall, strawgas, strawwire;
    std::string inwallfile, ecupwallfile, ecdwwallfile;
    while (std::getline(tracker_stream, line)) {
      // skip comments and blank lines
      if (line.compare(0,1,comment) != 0 && line.size() > 0 ) {
        // strip leading whitespace
        line = std::regex_replace(line, std::regex("^ +"), "");
        std::istringstream iss(line);
        // first, overall geometry
        if(rmin < 0.0){
          iss >> rmin >> rmax >> zpos >> zhalf;
          cyl_ = HollowCylinder(rmin,rmax,zpos,zhalf);
        } else if(rcell < 0.0) {
          // then cell description
          iss >> orient >> ncells_ >> rcell >> lcell >> wthick >> rwire >> minpath_;
          // std::cout << "ncells " << ncells_  << " rcell " << rcell << " lcell " << lcell << " wthick " << wthick << " rwire " << rwire << std::endl;
          if(rcell<0.0)throw std::invalid_argument("Invalid cell parameters");
        } else if(strawwall=="") {
          iss >> strawwall >> strawgas >> strawwire;
          std::cout << "Straw materials: wall " << strawwall << " gas " << strawgas << " wire " << strawwire << std::endl;
          smat_ = new KinKal::StrawMaterial(matdbinfo, rcell, wthick, rwire, strawwall.c_str(), strawgas.c_str(), strawwire.c_str());
          // calculate the scatter fraction; these are used to keep the RMS unchanged
          corefrac_ = smat_->wallMaterial().scatterFraction();
          lowscat_ = sqrt( (1.0 +hiscat_*hiscat_*(corefrac_-1.0))/corefrac_);
          orientation_ = (CellOrientation)orient;
          // cell density (transverse to the cell)
          cellDensity_ = ncells_*2.0*rcell*lcell/cyl_.volume();
          // compute the total mass of 1 cell
          double cmass = smat_->wallMaterial().density()*2.0*M_PI*rcell*lcell*wthick
            + smat_->gasMaterial().density()*M_PI*rcell*rcell*lcell
            + smat_->wireMaterial().density()*M_PI*rwire*rwire*lcell;
//          std::cout << "cell mass = " << cmass << std::endl;
          density_ = cmass*ncells_/cyl_.volume();
        } else if(vdrift_ < 0.0){
          // hit properties
          iss >> vdrift_ >> vprop_ >> sigt_ >> sigl_ >> lrdoca_ >> hiteff_;
        } else if(inwallfile=="") {
        	iss >> inwallfile;
        	std::cout << "Inner wall present, description in: " << inwallfile << std::endl;
          inwall_ = new TrckrWall(matdbinfo, inwallfile, tr_);
        } else if(ecupwallfile=="") {
        	iss >> ecupwallfile;
            std::cout << "Upstream EndCap wall present, description in: " << ecupwallfile << std::endl;
            ecpwallups_ = new TrckrWall(matdbinfo, ecupwallfile, tr_);
        } else if(ecdwwallfile=="") {
          iss >> ecdwwallfile;
            std::cout << "Downstream EndCap wall present, description in: " << ecdwwallfile << std::endl;
            ecpwalldws_ = new TrckrWall(matdbinfo, ecdwwallfile, tr_);
        }
      }
    }
    std::cout << "Read Tracker from file " << fullfile << std::endl;
  }

  void Tracker::print(std::ostream& os ) const {
    std::cout << "Tracker Z between " << zMin() << " and " << zMax() << " Rho between " << rMin() << " and " << rMax()
      << " average material density (gm/mm^3)" << density_ << std::endl;
    std::cout << "scatter core fraction " << corefrac_ << " low factor " << lowscat_ << " high factor " << hiscat_ << std::endl;
    std::cout << "Cell density (cells/mm) " << cellDensity_ << " with " << ncells_ << " cells oriented ";
    if(orientation_ == azimuthal)
      std::cout << " azimuthally " << " and miminum path " << minpath_ << std::endl;
    else
      std::cout << " axially " << std::endl;
    std::cout << "Cell radius " << cellRadius() << " wall thickness " << smat_->wallThickness() << std::endl;
    std::cout << "Vdrift " << vdrift_ << " Vprop " << vprop_ << " transverse resolution " << sigt_
    << " longitudinal resolution " << sigl_ << " minimum LR doca " << lrdoca_
    << " hit efficiency " << hiteff_ << std::endl;
    if ( inwall_ ) { inwall_->print(os); }
    if ( ecpwallups_ ) { ecpwallups_->print(os); }
    if ( ecpwalldws_ ) { ecpwalldws_->print(os); }
  }

  Tracker::~Tracker() {
	  if ( inwall_ ) { delete inwall_; }
	  if ( ecpwallups_ ) { delete ecpwallups_; }
	  if ( ecpwalldws_ ) { delete ecpwalldws_; }
  }
}


