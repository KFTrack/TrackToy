//
//  Compute the range and other properties of a muon of a givem momentum in a given material description file from the PDG
//  https://pdg.lbl.gov/2021/AtomicNuclearProperties/
//
#include "TrackToy/General/MuonRange.hh"
#include "TrackToy/General/FileFinder.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <regex>
#include <stdexcept>

namespace TrackToy {
  MuonRange::MuonRange(const char* rangefile, double density) : density_(density) {
    FileFinder filefinder;
    std::string fullfile = filefinder.fullFile(rangefile);
    //    std::cout << "Reading file " << fullfile << std::endl;
    // read the file
    std::ifstream range_stream(fullfile);
    std::string line;
    static std::string comment("#");
    double energy, momentum, eion, ebrems, epair, ephoto, erad, dEdx, range, delta, beta, dEdx_R;
    while (std::getline(range_stream, line)) {
      // skip comments and blank lines
      if (line.compare(0,1,comment) != 0 && line.size() > 0 ) {
        // strip leading whitespace
        line = std::regex_replace(line, std::regex("^ +"), "");
        std::istringstream iss(line);
        while (iss >> energy >> momentum >> eion >> ebrems >> epair >> ephoto >> erad >> dEdx >> range >> delta >> beta >> dEdx_R) {
          range_.push_back(RangeData(energy,momentum,range));
        }
      }
    }
    if(range_.size() ==0)throw std::invalid_argument("No range data read");

  }

  double  MuonRange::rangeEnergy(double energy) const {
    double retval(0.0);
    // find the energy range
    size_t imax;
    if(energy <= range_.front().energy_) {
      imax = 1;
    } else if (energy >= range_.back().energy_) {
      imax = range_.size()-1;
    } else {
      // loop
      imax =1;
      while(energy > range_[imax].energy_ )++imax;
    }
    size_t imin = imax-1;
    double dRdE = (range_[imax].range_-range_[imin].range_)/(range_[imax].energy_-range_[imin].energy_);
    retval = std::max(0.0,(range_[imin].range_ + dRdE*(energy-range_[imin].energy_)));
    return 10.0*retval/density_; // convert to mm
  }

  double  MuonRange::rangeMomentum(double momentum) const {
    double retval(0.0);
    // find the momentum range
    size_t imax;
    if(momentum <= range_.front().momentum_) {
      imax = 1;
    } else if (momentum >= range_.back().momentum_) {
      imax = range_.size()-1;
    } else {
      // loop
      imax =1;
      while(momentum > range_[imax].momentum_ )++imax;

    }
    size_t imin = imax-1;
    double dRdM = (range_[imax].range_-range_[imin].range_)/(range_[imax].momentum_-range_[imin].momentum_);
    retval = std::max(0.0,(range_[imin].range_ + dRdM*(momentum-range_[imin].momentum_)));
    return 10.0*retval/density_; // convert to mm
  }

}
