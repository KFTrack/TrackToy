//
//  Test Ce tracks from stopped muons
//
#include "KinKal/General/AxialBFieldMap.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/MatEnv/SimpleFileFinder.hh"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/Fit/Track.hh"
#include "KinKal/Examples/ScintHit.hh"
#include "TrackToy/General/FileFinder.hh"
#include "TrackToy/General/TrajUtilities.hh"
#include "TrackToy/Detector/HollowCylinder.hh"
#include "TrackToy/Detector/Target.hh"
#include "TrackToy/Detector/IPA.hh"
#include "TrackToy/Detector/Tracker.hh"
#include "TrackToy/Detector/Calorimeter.hh"
#include "TrackToy/Spectra/CeMinusSpectrum.hh"
#include "TrackToy/Spectra/DIOSpectrum.hh"
#include "TrackToy/Test/TrkInfo.hh"
#include "TFile.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TAxis3D.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TView.h"
#include "TTUBE.h"
#include "TBRIK.h"
#include "TNode.h"
#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <chrono>
#include <ctime>

//using namespace TrackToy;
using namespace std;
using namespace TrackToy;
using namespace KinKal;

void print_usage() {
  printf("Usage: CeTrackTest --mustopsfile s --rmue f --bfield s --trkfield i --targetfile s --trackerfile s --ipafile s --calofile s --fitschedule s--extschedule s --process s --endpoint f --endrange f --lifetime f --tol f  --seed i --npts i --ntrks i --draw i --ttree i --tfile s --minnhits i --npot i --saveall i --faildetail i \n");
}

int makeConfig(string const& cfile, KinKal::Config& config) {
  string fullfile;
  if(strncmp(cfile.c_str(),"/",1) == 0) {
    fullfile = string(cfile);
  } else {
    if(const char* source = std::getenv("TRACKTOY_SOURCE_DIR")){
      fullfile = string(source) + ("/Data/") + string(cfile);
    } else {
      cout << "TRACKTOY_SOURCE_DIR not defined" << endl;
      return -1;
    }
  }
  std::ifstream ifs (fullfile, std::ifstream::in);
  if ( (ifs.rdstate() & std::ifstream::failbit ) != 0 ){
    std::cerr << "Error opening " << fullfile << std::endl;
    return -1;
  }
  string line;
  int plevel(-1);
  while (getline(ifs,line)){
    if(line.size() > 0 && strncmp(line.c_str(),"#",1)!=0){
      istringstream ss(line);
      if(plevel < 0) {
        ss >> config.maxniter_ >> config.dwt_ >> config.convdchisq_ >> config.divdchisq_ >>
        config.pdchi2_ >> config.tbuff_ >> config.tol_ >> config.minndof_ >> config.bfcorr_ >>
        plevel;
        config.plevel_ = Config::printLevel(plevel);
      } else {
        double temp, mindoca(-1.0),maxdoca(-1.0), minprob(-1.0);
        ss >> temp >> mindoca >> maxdoca >> minprob;
        MetaIterConfig mconfig(temp);
        if(mindoca >0.0 || maxdoca > 0.0){
          // setup and insert the updater
          cout << "SimpleWireHitUpdater with mindoca " << mindoca << " maxdoca " << maxdoca << " minprob " << minprob << endl;
          SimpleWireHitUpdater updater(mindoca,maxdoca,minprob);
          mconfig.addUpdater(std::any(updater));
        }
        config.schedule_.push_back(mconfig);
      }
    }
  }
  return 0;
}

int main(int argc, char **argv) {
  using KTRAJ=LoopHelix;
  using PKTRAJ = ParticleTrajectory<KTRAJ>;
  using KKTRK = Track<KTRAJ>;
  using KKHIT = Measurement<KTRAJ>;
  using KKMAT = Material<KTRAJ>;
  using KKBF = BField<KTRAJ>;
  using Clock = std::chrono::high_resolution_clock;
  int ntrks(100);
  string bfile("Data/DSMapDump.dat"), mstops("MDC2020n_10pc"), targetfile("Data/Mu2eTarget.dat"), trackerfile("Data/Mu2eTracker.dat");
  string calofile("Data/Mu2eCalo.dat");
  string ipafile("Data/Mu2e_IPA.dat");
  string efile_my("Data/EStar_Mylar.dat"); // should come from tracker FIXME
  string sfile("Schedule_driftfit.txt"), extfile; // fit schedule
  string process("CeMinus");
  string diofile("Data/DIOAl_fine.dat"); // this should be a parameter FIXME
  double mustopeff; // mu stops/POT.  This comes from running MuBeam, which also produces the Mustops.root file.  These MUST be consistent to get physically correct results!!
  double npot(3.6e20); // total number of POTs
  double rmue(1e-16); // set as needed
  double endpoint(105.0), lifetime(864.0); // these should be specified by target material FIXME
  double decayfrac(0.391); // target material specific FIXME!
  double endrange(5.0); // spectrum generation range
  double tol(1e-4);
  double emass(0.511); //electron
  size_t npts(5000);
  bool draw(false), ttree(true), saveall(false);
  string tfile("");
  int trkfieldtype(0);
  unsigned minnhits(15); // minimum # of hits
  //  double mine(90.0); // minimum energy to simulate
  // ttree variables
  TTree* trktree_(0);
  float targetde_, ipade_, trackerde_;
  float origine_, targete_, ipae_, trackere_;
  float trkrECUPSde_, trkrECDWSde_, trkrIWentrde_, trkrIWotrde_, trkrWentrde_;
  float trkrECUPSe_, trkrECDWSe_, trkrIWentre_, /*trkrIWotre_,*/ trkrWentre_;
  float bfrtrcke_;
//  double bfrtrktime_;
  VEC3 originpos_, originmom_;
  float origintime_;
  int ntarget_, nipa_;
  int ntrkrECUPS_, ntrkrECDWS_, ntrkrIWentr_, ntrkrIWotr_, ntrkrWentr_;
  int itrk_;
  double weight_; // event weight
  TrkInfo tinfo_;
  VEC3 kkentmom_, kkmidmom_, kkextmom_;
  float kkentmomerr_, kkmidmomerr_, kkextmomerr_;
  float kkentt0_, kkmidt0_, kkextt0_;
  VEC3 mcentmom_, mcmidmom_, mcextmom_;
  VEC3 kkentpos_, kkmidpos_, kkextpos_;
  VEC3 mcentpos_, mcmidpos_, mcextpos_;
  // fit parameters
  KinKal::DVEC sigmas(0.5, 0.5, 0.5, 0.5, 0.002, 0.5); // expected parameter sigmas for loop helix
  double seedsmear(1.0);
  unsigned seed(0);
  Config::printLevel faildetail(Config::none);

  static struct option long_options[] = {
    {"mustopsfile",     required_argument, 0, 'm' },
    {"bfield",     required_argument, 0, 'F' },
    {"trkfield",     required_argument, 0, 'k' },
    {"targetfile",     required_argument, 0, 't' },
    {"trackerfile",     required_argument, 0, 'T' },
    {"ipafile",     required_argument, 0, 'i' },
    {"calofile",     required_argument, 0, 'c' },
    {"fitschedule",     required_argument, 0, 'X' },
    {"extschedule",     required_argument, 0, 'z' },
    {"process",     required_argument, 0, 's'  },
    {"endpoint",     required_argument, 0, 'e' },
    {"endrange",     required_argument, 0, 'E' },
    {"tol",     required_argument, 0, 'x' },
    {"seed",     required_argument, 0, 'I'  },
    {"ntrks",     required_argument, 0, 'n'  },
    {"npts",     required_argument, 0, 'q'  },
    {"draw",     required_argument, 0, 'd'  },
    {"saveall",     required_argument, 0, 'S'  },
    {"ttree",     required_argument, 0, 'r'  },
    {"tfile",     required_argument, 0, 'R'  },
    {"minnhits",     required_argument, 0, 'N' },
    {"faildetail",     required_argument, 0, 'f' },
    {"npot",     required_argument, 0, 'P' },
    {"rmue",     required_argument, 0, 'M' },
    {NULL, 0,0,0}
  };
  int opt;
  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'F' : bfile = string(optarg);
                 break;
      case 'k' : trkfieldtype = atoi(optarg);
                 break;
      case 'm' : mstops = string(optarg);
                 break;
      case 't' : targetfile = string(optarg);
                 break;
      case 'T' : trackerfile = string(optarg);
                 break;
      case 'i' : ipafile = string(optarg);
                 break;
      case 'c' : calofile = string(optarg);
                 break;
      case 'X' : sfile = string(optarg);
                 break;
      case 'z' : extfile = string(optarg);
                 break;
      case 's' : process = string(optarg);
                 break;
      case 'e' : endpoint = atof(optarg);
                 break;
      case 'E' : endrange = atof(optarg);
                 break;
      case 'x' : tol = atof(optarg);
                 break;
      case 'I' : seed = atoi(optarg);
                 break;
      case 'q' : npts = atoi(optarg);
                 break;
      case 'n' : ntrks = atoi(optarg);
                 break;
      case 'd' : draw = atoi(optarg);
                 break;
      case 'S' : saveall = atoi(optarg);
                 break;
      case 'r' : ttree = atoi(optarg);
                 break;
      case 'R' : tfile  = string(optarg);
                 break;
      case 'N' : minnhits = atoi(optarg);
                 break;
      case 'f' : faildetail = (Config::printLevel)atoi(optarg);
                 break;
      case 'P' : npot = atof(optarg);
                 break;
      case 'M' : rmue = atof(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }
  auto now = Clock::now();
  // open the input muonstops file
  string mfile =mstops + string("MuStops.root");
  TFile* mustopsfile = TFile::Open(mfile.c_str(),"READ");
  if(!mustopsfile){
    cout << "MuStop file " << mfile << " not found: did you forget to run MuStops_test?   terminating" << endl;
    return 1;
  }
  // find the TTree in the pfile
  TTreeReader reader("MuStops",mustopsfile);
  TTreeReaderValue<VEC4> mustoppos(reader, "Pos");
  TTree* mtree = (TTree*)mustopsfile->Get("MuStops");
  if(ntrks<0)ntrks = mtree->GetEntries();
// get stop efficiency (per pot)
  mfile = mstops + string("MuonStopEff.txt");
  ifstream mstopeffstream(mfile,std::ios_base::in);
  if(mstopeffstream.fail()){
    std::string errmsg = std::string("File doesn't exist " )+ mfile;
    throw std::invalid_argument(errmsg.c_str());
  }
  mstopeffstream >> mustopeff;
  cout << "MuStops TTree from file " << mfile << " has " << mtree->GetEntries() << " Entries, Stops/POT " << mustopeff << endl;
  // build the materials database
  MatEnv::SimpleFileFinder ttfinder(std::string("TRACKTOY_SOURCE_DIR"),std::string("/Data/"));
  cout << "Using Materials file " << ttfinder.matMtrDictionaryFileName() << endl;
  MatEnv::MatDBInfo matdb_(ttfinder,MatEnv::DetMaterial::moyalmean);

  // random; set the seed from the clock
  if(seed == 0){
    auto tdiff = Clock::now() - now;
    double nsecs = std::chrono::duration_cast<std::chrono::nanoseconds>(tdiff).count();
    seed = static_cast<unsigned>(rint(nsecs));
  }
  cout << "Random seed " << seed << endl;
  TRandom3 tr_(seed); // random number generator
  // setup target
  Target target(matdb_,targetfile,tr_);
//  Target target(targetfile,tr_);
  target.print(cout);
  // setup ipa
  IPA ipa(matdb_,ipafile,tr_);
  ipa.print(cout);
  // setup tracker
  Tracker tracker(matdb_,trackerfile,tr_);
  tracker.print(cout);
  EStar trackerEStar(efile_my); // this is no longer needed, except for comparisons FIXME
  // setup calo
  Calorimeter calo(calofile,tr_);
  calo.print(cout);
  // setup BField
  FileFinder filefinder;
  std::string fullfile = filefinder.fullFile(bfile);
  cout << "Building BField from file " << fullfile << endl;
  auto bfield = new AxialBFieldMap(fullfile);
  KinKal::BFieldMap* trkfield = bfield;
  auto bent = bfield->fieldVect(VEC3(0.0,0.0,tracker.zMin()));
  if(trkfieldtype == 1){
    cout << "Using fixed BField for tracker"<< endl;
    auto bent = bfield->fieldVect(VEC3(0.0,0.0,tracker.zMin()));
    trkfield = new UniformBFieldMap(bent.Z());
  } else if(trkfieldtype == 2) {
    cout << "Using gradient for tracker"<< endl;
    trkfield = new GradientBFieldMap(bent.Z(),0.96*bent.Z(),tracker.zMin(),tracker.zMax());
  } else {
    cout << "Using full field for tracker"<< endl;
  }
  bfield->print(cout);
//  auto trkfield = bfield;
  // setup fit configuration
  Config config;
  makeConfig(sfile,config);
  cout << "Fit configuration " << config << endl;
  // setup extension (if provided).
  Config extconfig;
  if(extfile.size() > 0){
    makeConfig(extfile,extconfig);
    cout << "Ext configuration " << extconfig << endl;
  }
  // setup spectrum
  Spectrum* spectrum(0);
  bool flat(false);
  double wfac = npot*mustopeff/(double)ntrks;
  if(process == "CeMinus"){
    // endpoint energy should come from target material TODO
    CeMinusSpectrumParams ceparams(endpoint);
    spectrum = new CeMinusSpectrum(ceparams);
    weight_ = (1.0-decayfrac)*rmue*wfac;
  } else if (process == "DIO") {
    spectrum = new DIOSpectrum(diofile.c_str(),endpoint-endrange,endpoint);
    weight_ = decayfrac*wfac/spectrum->normalization();
  } else if (process == "FlatDIO") {
    spectrum = new DIOSpectrum(diofile.c_str(),endpoint-endrange,endpoint);
    flat = true;
  } else {
    cout << "Unknown process " << process << ": aborting" << endl;
    return -2;
  }
  // timing

  // output file
  string ofile= tfile+process+string("Tracks.root");
  TFile trkfile(ofile.c_str(),"RECREATE");
  TH1F* ipade = new TH1F("ipde","IPA dE",100,-3.0,0.0);
  TH1F* tarde = new TH1F("tarde","Target dE;dE (MeV)",100,-3.0,0.0);
  TH1F* trkde = new TH1F("trkde","Tracker dE;dE (MeV)",100,-3.0,0.0);
  TH1F* trknc = new TH1F("trknc","Tracker N Cells;N Cells",100,0.001,100.0);
  if(ttree){
    trktree_ = new TTree("trks","trks");
    trktree_->Branch("itrk",&itrk_,"itrk/I");
    trktree_->Branch("weight",&weight_,"weight/D");
    trktree_->Branch("tinfo",&tinfo_);
    trktree_->Branch("origine",&origine_,"origine/F");
    trktree_->Branch("originpos",&originpos_);
    trktree_->Branch("originmom",&originmom_);
    trktree_->Branch("origintime",&origintime_,"origintime/F");
    trktree_->Branch("ntarget",&ntarget_,"ntarget/I");
    trktree_->Branch("targetde",&targetde_,"targetde/F");
    trktree_->Branch("targete",&targete_,"targete/F");
    trktree_->Branch("ipae",&ipae_,"ipae/F");
    trktree_->Branch("nipa",&nipa_,"nipa/I");
    trktree_->Branch("ipade",&ipade_,"ipade/F");
    trktree_->Branch("trkrECUPSe",&trkrECUPSe_,"trkrECUPSe/F");
    trktree_->Branch("ntrkrECUPS",&ntrkrECUPS_,"ntrkrECUPS/I");
    trktree_->Branch("trkrECUPSde",&trkrECUPSde_,"trkrECUPSde/F");
    trktree_->Branch("trkrECDWSe",&trkrECDWSe_,"trkrECDWSe/F");
    trktree_->Branch("ntrkrECDWS",&ntrkrECDWS_,"ntrkrECDWS/I");
    trktree_->Branch("trkrECDWSde",&trkrECDWSde_,"trkrECDWSde/F");
    trktree_->Branch("trkrIWentre",&trkrIWentre_,"trkrIWentre/F");
    trktree_->Branch("ntrkrIWentr",&ntrkrIWentr_,"ntrkrIWentr/I");
    trktree_->Branch("trkrIWentrde",&trkrIWentrde_,"trkrIWentrde/F");
//    trktree_->Branch("trkrIWotre",&trkrIWotre_,"trkrIWotre/F");
    trktree_->Branch("ntrkrIWotr",&ntrkrIWotr_,"ntrkrIWotr/I");
    trktree_->Branch("trkrIWotrde",&trkrIWotrde_,"trkrIWotrde/F");
    trktree_->Branch("trkrWentre",&trkrWentre_,"trkrWentre/F");
    trktree_->Branch("ntrkrWentr",&ntrkrWentr_,"ntrkrWentr/I");
    trktree_->Branch("trkrWentrde",&trkrWentrde_,"trkrWentrde/F");
    trktree_->Branch("trackere",&trackere_,"trackere/F");
    trktree_->Branch("trackerde",&trackerde_,"trackerde/F");
    trktree_->Branch("mcentmom",&mcentmom_);
    trktree_->Branch("mcmidmom",&mcmidmom_);
    trktree_->Branch("mcextmom",&mcextmom_);
    trktree_->Branch("mcentpos",&mcentpos_);
    trktree_->Branch("mcmidpos",&mcmidpos_);
    trktree_->Branch("mcextpos",&mcextpos_);
    trktree_->Branch("kkentmom",&kkentmom_);
    trktree_->Branch("kkmidmom",&kkmidmom_);
    trktree_->Branch("kkextmom",&kkextmom_);
    trktree_->Branch("kkentmomerr",&kkentmomerr_,"kkentmomerr/F");
    trktree_->Branch("kkmidmomerr",&kkmidmomerr_,"kkmidmomerr/F");
    trktree_->Branch("kkextmomerr",&kkextmomerr_,"kkextmomerr/F");
    trktree_->Branch("kkentt0",&kkentt0_,"kkentt0/F");
    trktree_->Branch("kkmidt0",&kkmidt0_,"kkmidt0/F");
    trktree_->Branch("kkextt0",&kkextt0_,"kkextt0/F");
    trktree_->Branch("kkentpos",&kkentpos_);
    trktree_->Branch("kkmidpos",&kkmidpos_);
    trktree_->Branch("kkextpos",&kkextpos_);
  }

  std::vector<TPolyLine3D*> plhel;
  // loop over stops
  int icolor(kBlue);

  unsigned nfit(0), nfail(0), ndiv(0), npdiv(0), nlow(0), nconv(0), nuconv(0);
  itrk_ = 0;
  auto start = Clock::now();
  while (itrk_ < ntrks) {
    ++itrk_;
    if(fmod(itrk_,ntrks/10) == 0)cout << "Processing track " << itrk_ << endl;
    tinfo_.reset();
    if(!reader.Next()){
      reader.Restart();
      if(!reader.Next()){
        cout << "Unable to rewind file" << mustopsfile << endl;
        return -2;
      }
    }
    //    cout << "Track " << itrk_ << endl;
    // reset tree variables
    targetde_ = ipade_ = trackerde_ = 0.0;
    targete_ = ipae_ = 0.0;
    ntarget_ = nipa_ = -1;
    trkrECUPSde_ = 0.0;
    trkrECUPSe_ = 0.0;
    ntrkrECUPS_ = -1;
    kkentmom_ = kkmidmom_ = kkextmom_ = VEC3();
    kkentmomerr_ = kkmidmomerr_ = kkextmomerr_ = -1.0;
    kkentt0_ = kkmidt0_ = kkextt0_ = 0.0;
    mcentmom_ = mcmidmom_ = mcextmom_ = VEC3();
    kkentpos_ = kkmidpos_ = kkextpos_ = VEC3();
    mcentpos_ = mcmidpos_ = mcextpos_ = VEC3();

    // generate a random energy
    if(flat){
      origine_ = tr_.Uniform(endpoint-endrange,endpoint);
      weight_ = wfac*spectrum->rate(origine_);
    } else {
      double prob = tr_.Uniform(0.0,1.0);
      origine_ = spectrum->sample(prob);
    }
    // generate a random muon decay time; this should come from the target FIXME
    double tdecay = tr_.Exp(lifetime);
    // generate random phi and cos(theta)
    double phi = tr_.Uniform(-M_PI,M_PI);
    double cost = tr_.Uniform(-1.0,1.0);
    double sint = sqrt(1.0-cost*cost);
    VEC4 const& pos4 = *mustoppos;
    originpos_ = pos4.Vect();
    origintime_ = tdecay+pos4.T(); // add decay time to stopping time
    double mom = sqrt(origine_*origine_ - emass*emass);
    originmom_ = VEC3(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost);
    ParticleState cestate(originpos_,originmom_,origintime_,emass,-1);
    TimeRange range(cestate.time(),cestate.time()+1000.0); // much longer than physical: is truncated later
    auto bstart = bfield->fieldVect(cestate.position3());
    KTRAJ lhelix(cestate,bstart,range);
    // sim tolerance is smaller than fit
    double mctol = tol/3;
    //    cout << "Initial trajectory " << lhelix << endl;
    // initialize piecetraj
    PKTRAJ mctraj(lhelix);
    TimeRanges targetinters, ipainters, trackerinters;
    TimeRanges trckrECUPSinters, trckrECDWSinters, trckrIWinters;
    // extend through the target
    if(target.extendTrajectory(*bfield,mctraj,targetinters,mctol)){
      //      cout << "Extended to target " << targetinters.size() << endl;
      //      mctraj.print(cout,2);
      targete_ = mctraj.energy(mctraj.range().end());
      targetde_ = targete_ - origine_;
      ntarget_ = targetinters.size();
      // extend through the IPA
      if(ipa.extendTrajectory(*bfield,mctraj,ipainters,mctol)){
        //        cout << "Extended to ipa " << ipainters.size() << endl;
        //        mctraj.print(cout,2);
        nipa_ = ipainters.size();
        ipae_ = mctraj.energy(mctraj.range().end());
        ipade_ = ipae_ - targete_;

        bfrtrcke_ = ipae_;

        // extend  to the tracker entrance
        double trckZmin = tracker.zMin();
        if(tracker.EndCapWallUpStIsPresent()){
          trckZmin = tracker.EndCapWallUpSt().cylinder().zmin();
         }
        extendZ(mctraj,*bfield, trckZmin, mctol);
//        bfrtrktime_ = mctraj.back().range().begin();

        // if the tracker field is different from the general field, change the trajector bnom
        if(trkfield != bfield){
          double tent = ztime(mctraj,mctraj.back().range().begin(),trckZmin);
          auto pstate = mctraj.back().state(tent);
          auto pos = pstate.position3();
          auto bend = trkfield->fieldVect(pos);
          KTRAJ endtraj(pstate,bend,TimeRange(tent,mctraj.range().end()));
          mctraj.append(endtraj);
        }

        ntrkrECUPS_=ntrkrECDWS_=ntrkrIWentr_=ntrkrIWotr_=ntrkrWentr_=-1;
        trkrECUPSe_=trkrECDWSe_=trkrIWentre_=/*trkrIWotre_=*/trkrWentre_=0.0;
        trkrECUPSde_=trkrECDWSde_=trkrIWentrde_=trkrIWotrde_=trkrWentrde_=0.0;
        if (tracker.EndCapWallUpStIsPresent() && tracker.EndCapWallUpSt().extendTrajectory(*bfield,mctraj,trckrECUPSinters,mctol)){
          ntrkrECUPS_ = trckrECUPSinters.size();
          trkrECUPSe_ = mctraj.energy(mctraj.range().end());
          trkrECUPSde_ = trkrECUPSe_ - ipae_;
          if (ntrkrECUPS_>0) {
            trkrWentre_ = trkrECUPSe_;
            trkrWentrde_ = trkrECUPSde_;
            bfrtrcke_ = trkrECUPSe_;
//            bfrtrktime_ = trckrECUPSinters.at(0).end();
          }
//          std::cout<<"EndCap Up intersection found, n. "<<trckrECUPSinters.size()<<std::endl;
        }
        if (tracker.InWallIsPresent() && ntrkrECUPS_<1 && tracker.InWall().extendTrajectory(*bfield,mctraj,trckrIWinters,mctol)) {
          if (trckrIWinters.size()>0) {
            ntrkrIWentr_=1;
            auto iwlinter = trckrIWinters.at(0);
            trkrIWentre_ = mctraj.energy(iwlinter.end());
            trkrIWentrde_ = trkrIWentre_ - ipae_;
            trkrWentre_ = trkrIWentre_;
            trkrWentrde_ = trkrIWentrde_;
            bfrtrcke_ = trkrIWentre_;
//            bfrtrktime_ = iwlinter.end();
//            auto posIn = mctraj.position3(iwlinter.begin());
//            auto posEx = mctraj.position3(iwlinter.end());
//            std::cout<<"Inner Wall intersection at:"<<std::endl;
//            std::cout<<"entering pos "<<posIn<<" rad "<<posIn.Rho()<<" exiting "<<posEx<<" rad "<<posEx.Rho()<<std::endl;
          }

        }

        // now create hits and straw intersections
        std::vector<std::shared_ptr<Hit<KTRAJ>>> hits;
        std::vector<std::shared_ptr<ElementXing<KTRAJ>>> xings;
        double speed = mctraj.speed(mctraj.range().end());
//        // if the tracker field is different from the general field, change the trajector bnom
//        if(trkfield != bfield){
//          double tent = ztime(mctraj,mctraj.back().range().begin(),tracker.zMin());
//          auto pstate = mctraj.back().state(tent);
//          auto pos = pstate.position3();
//          auto bend = trkfield->fieldVect(pos);
//          KTRAJ endtraj(pstate,bend,TimeRange(tent,mctraj.range().end()));
//          mctraj.append(endtraj);
//        }
        tracker.simulateHits(*trkfield,mctraj,hits,xings,trackerinters,trckrIWinters,trckrECDWSinters,mctol);
        tinfo_.ncells = xings.size();
        tinfo_.ntrkhits = hits.size();
        tinfo_.narcs = trackerinters.size();
        if(hits.size() >= minnhits){
// calcluate the estart energy loss

          if (trckrIWinters.size()>0) {
            int iwintotfst=0;
            if (ntrkrECUPS_<1) { iwintotfst=1; }
            ntrkrIWotr_ = trckrIWinters.size()-iwintotfst;
            for (int iwint=iwintotfst; iwint<trckrIWinters.size(); ++iwint) {
              auto iwlinter = trckrIWinters.at(iwint);
              trkrIWotrde_ += mctraj.energy(iwlinter.end())-mctraj.energy(iwlinter.begin());
            }
          }
          if (trckrECDWSinters.size()>0) {
            auto ewlinter = trckrIWinters.at(0);
            ntrkrECDWS_ = trckrECDWSinters.size();
            trkrECDWSe_ = mctraj.energy(ewlinter.end());
            trkrECDWSde_ = trkrECDWSe_ - mctraj.energy(ewlinter.begin());
          }

          double trackerpath(0.0);
          for(auto const& inter : trackerinters) { trackerpath += speed*inter.range(); }
//          double ke = cestate.energy() - cestate.mass();
//          trackerde_ = -100*trackerEStar.dEIonization(ke)*tracker.density()*trackerpath; // unit conversion
//          trackere_ = mctraj.energy(mctraj.range().end());
          trackere_ = mctraj.energy(trackerinters.back().end());
          trackerde_ = trackere_ - bfrtrcke_;//ipae_;
          tarde->Fill(targetde_);
          ipade->Fill(ipade_);
          trkde->Fill(trackerde_);
          trknc->Fill(tinfo_.ncells);
          // add calo hit
          tinfo_.ncalohits = calo.simulateHits(*trkfield,mctraj,hits,mctol);
          // truncate the true trajectory
          mctraj.setRange(TimeRange(mctraj.range().begin(),mctraj.back().range().begin()+0.1));
          // get the true times at entrance and exit
          double mctent = ztime(mctraj,trackerinters.front().begin(),tracker.zMin());
          double mctext = ztime(mctraj,trackerinters.back().end(),tracker.zMax());
          double mctmid = ztime(mctraj,0.5*(mctent+mctext),tracker.zMid());
          // record true momentum at tracker entranc, mid, and exit
          auto entstate = mctraj.stateEstimate(mctent);
          auto midstate = mctraj.stateEstimate(mctmid);
          auto extstate = mctraj.stateEstimate(mctext);
          mcentmom_ = entstate.momentum3();
          mcmidmom_ = midstate.momentum3();
          mcextmom_ = extstate.momentum3();
          mcentpos_ = entstate.position3();
          mcmidpos_ = midstate.position3();
          mcextpos_ = extstate.position3();

          // reconstruct track from hits and material info
          auto seedtraj = mctraj.nearestPiece(mctmid);
          seedtraj.setRange(mctraj.range());
          for(size_t ipar=0; ipar < NParams(); ipar++){
            double perr = sigmas[ipar]*seedsmear;
            seedtraj.params().covariance()[ipar][ipar] = perr*perr;
            seedtraj.params().parameters()[ipar] += tr_.Gaus(0.0,perr);
          }
          KKTRK kktrk(config,*trkfield,seedtraj,hits,xings);
          if(extconfig.schedule().size()>0){
            std::vector<std::shared_ptr<Hit<KTRAJ>>> exthits;
            std::vector<std::shared_ptr<ElementXing<KTRAJ>>> extxings;
            kktrk.extend(extconfig,exthits,extxings);
          }
          nfit++;
          // fill hit counting
          for(auto const& hit : hits) {
            if(hit->active()){
              ++tinfo_.kknactive;
              auto swh = dynamic_cast<KinKal::SimpleWireHit<KTRAJ>*>(hit.get());
              if(swh != 0 && !swh->hitState().useDrift())++tinfo_.kknnull;
            }
          }

          // fill fit information
          auto const& fstat = kktrk.fitStatus();
          if(fstat.status_ == Status::failed)nfail++;
          if(fstat.status_ == Status::converged)nconv++;
          if(fstat.status_ == Status::unconverged)nuconv++;
          if(fstat.status_ == Status::lowNDOF)nlow++;
          if(fstat.status_ == Status::diverged)ndiv++;
          if(fstat.status_ == Status::paramsdiverged)npdiv++;
          tinfo_.kkstatus = fstat.status_;
          tinfo_.kkchisq = fstat.chisq_.chisq();
          tinfo_.kkprob = fstat.chisq_.probability();
          tinfo_.kkndof = fstat.chisq_.nDOF();
          for(auto const& eff: kktrk.effects()) {
            const KKHIT* kkhit = dynamic_cast<const KKHIT*>(eff.get());
            const KKBF* kkbf = dynamic_cast<const KKBF*>(eff.get());
            const KKMAT* kkmat = dynamic_cast<const KKMAT*>(eff.get());
            if(kkhit != 0){
              tinfo_.kknhit++;
            } else if(kkmat != 0){
              tinfo_.kknmat++;
            } else if(kkbf != 0){
              tinfo_.kknbf++;
            }
          }
          if(fstat.usable()){
            auto const& kktraj = kktrk.fitTraj();
            double kktent = ztime(kktraj,kktraj.range().begin(),tracker.zMin());
            double kktmid = ztime(kktraj,kktraj.range().mid(),tracker.zMid());
            double kktext = ztime(kktraj,kktraj.range().end(),tracker.zMax());
            // t0 comes from mid section
            auto const& kkmidtraj = kktraj.nearestPiece(kktmid);
            tinfo_.kkt0err = sqrt(kkmidtraj.paramVar(kkmidtraj.t0Index()));
            auto entstate = kktraj.stateEstimate(kktent);
            auto midstate = kktraj.stateEstimate(kktmid);
            auto extstate = kktraj.stateEstimate(kktext);
            kkentmom_ = entstate.momentum3();
            kkmidmom_ = midstate.momentum3();
            kkextmom_ = extstate.momentum3();
            kkentmomerr_ = sqrt(entstate.momentumVariance());
            kkmidmomerr_ = sqrt(midstate.momentumVariance());
            kkextmomerr_ = sqrt(extstate.momentumVariance());
// find time a crossing the entrance, mid, and exit z
            kkentt0_ = ztime(kktraj,kktraj.range().begin(),tracker.zMin());
            kkmidt0_ = ztime(kktraj,kktraj.range().mid(),tracker.zMid());
            kkextt0_ = ztime(kktraj,kktraj.range().end(),tracker.zMax());
            kkentpos_ = entstate.position3();
            kkmidpos_ = midstate.position3();
            kkextpos_ = extstate.position3();
          } else if(faildetail > Config::none) {
            cout << "Bad Fit event " << itrk_ << " status " << kktrk.fitStatus() << endl;
            cout << "True Traj " << mctraj << endl;
            cout << "Seed Traj " << seedtraj << endl;
            kktrk.print(cout,faildetail);
          }

          if(!saveall)trktree_->Fill();
          //
          if(draw){
            plhel.push_back(new TPolyLine3D(npts));
            plhel.back()->SetLineColor(icolor++%10);
            double tstart = mctraj.range().begin();
            double ts = mctraj.range().range()/(npts-1);
            VEC3 ppos;
            for(unsigned ipt=0;ipt<npts;ipt++){
              double t = tstart + ipt*ts;
              ppos = mctraj.position3(t);
              plhel.back()->SetPoint(ipt,ppos.X(),ppos.Y(),ppos.Z());
            }
          }
        } // particle hits the tracker
      } // particle stops in IPA
    } // particle exits the target going downstream
    if(saveall)trktree_->Fill();
  }

  auto stop = Clock::now();
  double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
  cout
    << nfit << " KinKal fits "
    << nconv << " Converged fits "
    << nuconv << " Unconverged fits "
    << nfail << " Failed fits "
    << nlow << " low NDOF fits "
    << ndiv << " Diverged fits "
    << npdiv << " ParameterDiverged fits "
    << duration/double(ntrks) << " nanoseconds/track" << endl;
  // Canvas of basic parameters
  TCanvas* trkcan = new TCanvas("Tracks");
  trkcan->Divide(2,2);
  trkcan->cd(1);
  tarde->Draw();
  trkcan->cd(2);
  ipade->Draw();
  trkcan->cd(3);
  trkde->Draw();
  trkcan->cd(4);
  trknc->Draw();

  trkcan->Write();

  if(draw){
    TCanvas* trkcan = new TCanvas("Tracks");
    TTUBE* ttarget= new TTUBE("ttarget","ttarget","void",target.cylinder().rmin(),target.cylinder().rmax(),target.cylinder().zhalf());
    ttarget->SetLineColor(kBlack);
    ttarget->SetLineWidth(4);
    ttarget->SetFillColorAlpha(kBlack, 0.5);
    ttarget->Draw();
    TTUBE* ttracker= new TTUBE("ttracker","ttracker","void",tracker.cylinder().rmin(),tracker.cylinder().rmax(),tracker.cylinder().zhalf());
    ttracker->SetLineColor(kRed);
    ttracker->SetLineWidth(4);
    ttracker->SetFillColorAlpha(kRed, 0.5);
    ttracker->Draw();
    TTUBE *DS  = new TTUBE("DS","DS","void",0.0,800,4000);
    DS->SetVisibility(0);
    DS->Draw();
    TNode *node1 = new TNode("NODE1","NODE1","DS");
    node1->cd();
    TNode* targetnode = new TNode("targetnode","targetnode","ttarget",0.0,0.0,-4300);
    TNode* trackernode = new TNode("trackernode","trackernode","ttracker",0.0,0.0,0);
    targetnode->Draw();
    trackernode->Draw();
    node1->cd();
    node1->Draw();
    // draw tracks
    for(auto const& ph : plhel) {
      ph->Draw();
    }
    // draw the origin and axes
    TAxis3D* rulers = new TAxis3D();
    rulers->GetXaxis()->SetAxisColor(kBlue);
    rulers->GetXaxis()->SetLabelColor(kBlue);
    rulers->GetYaxis()->SetAxisColor(kCyan);
    rulers->GetYaxis()->SetLabelColor(kCyan);
    rulers->GetZaxis()->SetAxisColor(kOrange);
    rulers->GetZaxis()->SetLabelColor(kOrange);
    rulers->SetAxisRange(target.cylinder().zmin(),tracker.cylinder().zmax(),"Z");
    rulers->Draw();
    trkcan->Write();
  }

  trkfile.Write();
  trkfile.Close();

  // now draw
  mustopsfile->Close();
  return 0;
}

