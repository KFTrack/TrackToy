//
//  Implementation of TrckrWall
//
#include "KinKal/Detector/MaterialXing.hh"
#include "TrackToy/Detector/TrckrWall.hh"
#include "TrackToy/General/FileFinder.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <stdexcept>
#include <string>
namespace TrackToy {
  TrckrWall::TrckrWall(MatEnv::MatDBInfo const& matdbinfo,std::string const& tgtfile,TRandom& tr) : type_(unknown), tr_(tr) {
    FileFinder filefinder;
    std::string fullfile = filefinder.fullFile(tgtfile);
    std::string line;
    static std::string comment("#");
    std::ifstream tgt_stream(fullfile);
    while (std::getline(tgt_stream, line)) {
      // skip comments and blank lines
      if (line.compare(0,1,comment) != 0 && line.size() > 0 ) {
        // strip leading whitespace
        line = std::regex_replace(line, std::regex("^ +"), "");
        std::istringstream iss(line);
        // first get type and material
        if(type_ == unknown){
          int type;
          std::string material;
          iss >> type >> material;
          type_ = (TrckrWallType)type;
          cyl_=nullptr;
          // lookup material
          mat_ = matdbinfo.findDetMaterial(material);
          if(mat_ == 0){
            std::string errmsg = std::string("Invalid Material ")+material;
            throw std::invalid_argument(errmsg);
          }
          // then geometry
        } else if (type_ == cylindrical) {
          double radius, rhalf, zpos, zhalf;
          iss >> radius >> rhalf >> zpos >> zhalf;
          if(radius < 0.0 || rhalf < 0.0 || zhalf < 0.0)throw std::invalid_argument("Invalid CylindricalShell parameters");
            cyl_ = new CylindricalShell(radius,rhalf,zpos,zhalf);
        } else if (type_ == disc) {
          double radius, rhalf, zpos, zhalf;
          iss >> radius >> rhalf >> zpos >> zhalf;
          if(radius < 0.0 || rhalf < 0.0 || zhalf < 0.0)throw std::invalid_argument("Invalid CylindricalShell parameters");
            cyl_ = new DiscShell(radius,rhalf,zpos,zhalf);
        }

      }
    }
    std::cout << "Read TrckrWall from file " << fullfile << std::endl;
  }

  void TrckrWall::print(std::ostream& os ) const {
    switch (type_ ) {
      default:
      case TrckrWall::cylindrical:
        std::cout << "Cylindrical TrckrWall with radius " << cyl_->radius() << " thickness " << 2*cyl_->rhalf() << " Zmid " << cyl_->zpos() << " Length " << 2*cyl_->zhalf();
        break;
      case TrckrWall::disc:
        std::cout << "Disc TrckrWall with radius " << cyl_->radius() << " thickness " << 2*cyl_->rhalf() << " Zmid " << cyl_->zpos() << " Length " << 2*cyl_->zhalf();
        break;
    }
    std::cout << " out of material " << mat_->name()  << std::endl;
  }

}

