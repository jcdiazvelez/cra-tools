#include <SimpleDST.h>
#include <SimpleTrigger.h>

#include "cuts.h"
#include <TChain.h>
#include <TH1D.h>
#include <TMath.h>
#include <TRandom.h>
#include <TStopwatch.h>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <healpix_cxx/fitshandle.h>
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/healpix_map_fitsio.h>
#include <healpix_cxx/pointing.h>

#include <photospline/splinetable.h>
#include <photospline/bspline.h>

#include <boost/program_options.hpp>

#include <astro/astro.h>
#include <astro/time.h>
#include <Direction.h>
#include <solardipole.h>

using namespace std;
namespace po = boost::program_options;

const double deg2rad = TMath::DegToRad();
const Double_t hour = 1 / 24.;
const Double_t minute = hour / 60.;
const Double_t second = minute / 60.;
const Double_t millisecond = 1e-3 * second;
const Double_t microsecond = 1e-6 * second;

void tScramble(po::variables_map vm, vector<string> inFiles);
//void tScramble(po::variables_map vm, const char* inFilesStr);


int main(int argc, char* argv[]) {

  po::options_description desc("Allowed options");
  desc.add_options()
      // Options used for all configurations
      ("help", "Produce help message")
      //("inFiles", po::value< vector<string> >()->multitoken(),"")
      ("batchFile", po::value<string>(), "Text file with input root filenames")
      ("batch_idx", po::value<string>(), "Line number to read in txt file")
      ("outBase", po::value<string>(), "Base name for outfile")
      ("yyyymmdd", po::value<string>(), "Desired date")
      ("config", po::value<string>(), "Detector configuration")
      ("nInt", po::value<string>(), "Integration time in hours")
      ("method", po::value<string>(), "Sidereal, Anti, Solar, Extended")
      ("sd", po::value<bool>(), "Correct for solar dipole for each event")
      ("sd2", po::value<bool>(), "apply 2nd-order dipole correction")
      // IceCube specific options
      ("spline", po::value<string>(), "File containing spline tables")
      // Icetop specific options
      ("filter", po::value<string>(), "Filter for IceTop data")
      ("comp", po::value< vector<string> >()->multitoken(), "Comp bins")
      ("sbins", po::value< vector<string> >()->multitoken(), "S125 bins")
      ("emin", po::value<string>(), "Minimum reconstructed energy")
      // Options for either detector
      ("ebins", po::value< vector<string> >()->multitoken(), "Energy bins")
  ;

  po::variables_map vm;
  // Disable short flags ('-') to allow use of negative signs in input
  po::store(po::parse_command_line(argc, argv, desc,
      po::command_line_style::unix_style ^ po::command_line_style::allow_short),
      vm);

  if (vm.count("help")) {
    cout << desc << "\n";
    return 1;
  }

  // Check for all necessary parameters
  string arr[] = {"outBase", "config", "nInt", "method", "yyyymmdd", "sd", "sd2"};
  int nKeys = 6;
  vector<string> keyParams(arr, arr+nKeys);
  for (unsigned i = 0; i < keyParams.size(); ++i) {
    if (not vm.count(keyParams[i])) {
      cerr << "\nUsage: " << keyParams[i] << " parameter not defined\n";
      return 1;
    }
  }

  // Read in filelist from batch_idx element of batchfile
  ifstream batchFile(vm["batchFile"].as<string>().c_str());
  string fileListStr;
  int batch_idx = atoi(vm["batch_idx"].as<string>().c_str());
  for (int i = 0; i < batch_idx; ++i) {
    getline(batchFile, fileListStr); 
  }

  getline(batchFile, fileListStr);
  batchFile.close();

  // Convert fileList string to vector
  istringstream iss(fileListStr);
  vector<string> fileList;
  copy(istream_iterator<string>(iss),
       istream_iterator<string>(),
       back_inserter(fileList));

  cout << "Input files are: " << endl;
  for (unsigned i = 0; i < fileList.size(); ++i) {
    cout << fileList[i] << endl;
  }


  if (vm.count("ebins")) {
    vector<string> test = vm["ebins"].as< vector<string> >();
    cout << "Ebin values:" << endl;
    for (unsigned i = 0; i < test.size(); ++i) {
      cout << " " << test[i];
    }
    cout << endl;
  }

  if (vm.count("sbins")) {
    vector<string> test = vm["sbins"].as< vector<string> >();
    cout << "S125 bin values:" << endl;
    for (unsigned i = 0; i < test.size(); ++i) {
      cout << " " << test[i];
    }
    cout << endl;
  }

  tScramble(vm, fileList);

  return 0;
}


void tScramble(po::variables_map vm, vector<string> inFiles_) {

  TStopwatch timer;
  timer.Start();

  // Read input parameters
  string outBase = vm["outBase"].as<string>();
  string config = vm["config"].as<string>();
  string filter = vm["filter"].as<string>();
  string method = vm["method"].as<string>();
  string yyyymmdd = vm["yyyymmdd"].as<string>();
  const Int_t nInt = atoi(vm["nInt"].as<string>().c_str());
  string detector = config.substr(0,2);
  bool sd = vm["sd"].as<bool>();
  bool sd2 = vm["sd2"].as<bool>();

  // Split infiles by space

  // Read in spline tables if provided
  //struct splinetable table;
  photospline::splinetable<> spline;
  if (vm.count("spline")) {
    string splineFile = vm["spline"].as<string>();
    spline.read_fits(splineFile.c_str());
  }

  // Energy binning setup
  unsigned nMaps = 1;
  vector<float> ebins, sbins;
  vector<string> ebinstr, sbinstr;
  if (vm.count("ebins")) {
    ebinstr = vm["ebins"].as<vector <string> >();
    for (unsigned n=0; n<ebinstr.size(); ++n) {
      ebins.push_back(atof(ebinstr[n].c_str()));
    }
    nMaps = ebins.size() - 1;
  }
  if (vm.count("sbins")) {
    sbinstr = vm["sbins"].as<vector <string> >();
    for (unsigned n=0; n<sbinstr.size(); ++n) {
      sbins.push_back(atof(sbinstr[n].c_str()));
    }
    nMaps = sbins.size() - 1;
  }


  int NSide = 64;

  // Allow for a map for each energy bin
  cout << "Number of maps: " << nMaps << endl;


  vector< Healpix_Map<float> > LocalMapInt(nMaps);
  vector< Healpix_Map<float> > LocalMap(nMaps);
  vector< Healpix_Map<float> > DataMapInt(nMaps);
  vector< Healpix_Map<float> > DataMap(nMaps);
  vector< Healpix_Map<float> > BGMap(nMaps);

  for (unsigned i = 0; i < nMaps; ++i) {
    LocalMapInt[i].SetNside(NSide, RING);
    LocalMapInt[i].fill(0.);
    LocalMap[i].SetNside(NSide, RING);
    LocalMap[i].fill(0.);
    DataMapInt[i].SetNside(NSide, RING);
    DataMapInt[i].fill(0.);
    DataMap[i].SetNside(NSide, RING);
    DataMap[i].fill(0.);
    BGMap[i].SetNside(NSide, RING);
    BGMap[i].fill(0.);
  }

  pointing sphereDir;
  int pixelID;

  stringstream sstr;
  sstr.str("");

  // Get info from previous- and next-day files
  int yy = atoi(yyyymmdd.substr(0, 4).c_str());
  int mm = atoi(yyyymmdd.substr(5, 2).c_str());
  int dd = atoi(yyyymmdd.substr(8, 2).c_str());
  astro::Time t(yy, mm, dd, 0, 0, 0);
  double MJD0 = t.GetMJD();

  const char* masterTree;
  const char* triggerTree;
  bool sundp2 = true;

  Config cfg(vm);

  if (cfg.detector == Config::IceCube) { 
        masterTree = "CutDST"; 
        triggerTree = "TDSTTriggers"; 
  } 
  if (cfg.detector == Config::IceTop) { 
        triggerTree = "";   // Unused? Will probably break IT functionality...  
        sundp2 = false; 
        if (cfg.cfg == Config::ITv3) { 
            masterTree = "MasterTree"; 
        } else { 
            masterTree = "master_tree"; 
        }
  }


  // Initialize the chain and read data
  TChain *cutDST = new TChain(masterTree);
  for (unsigned i = 0; i < inFiles_.size(); ++i) {
    cutDST->Add(inFiles_[i].c_str());
  }
  SimpleDST dst(cutDST, config);

  // Need to also initialize triggers if IC86-2016 or newer
  TChain *trigDST = new TChain(triggerTree);

  if (cfg.newConfig()) { 
    for (unsigned i = 0; i < inFiles_.size(); ++i) {
      trigDST->Add(inFiles_[i].c_str());
    }
  }
  SimpleTrigger dst_trig(trigDST);


  cout << "Number of chained files: " << cutDST->GetNtrees() << endl;

  Long64_t nEntries = cutDST->GetEntries();
  vector<Long64_t> nEvents(nMaps, 0);
  vector<Long64_t> nUsedEvents(nMaps, 0);

  const int nBGResample = 20;
  const double alpha = 1. / nBGResample;
  const double pi = TMath::Pi();

  Double_t startMJD, mjd2, mjd1=0;
  double zenith, azimuth, theta, phi, rndMJD;
  bool isGood;

  // Integration time
  const Double_t dt = nInt * hour;
  cout << "Integration time = " << nInt << " (hours) " << dt << " (day)\n";
  cout << "Reading " << nEntries << " entries...\n";

  mjd1 = MJD0;
  mjd2 = mjd1 + 1;
  startMJD = mjd1;

  // Setup histograms for storing time information
  vector<TH1D*> histMJD(nMaps);
  const char* histName;
  for (unsigned i = 0; i < nMaps; ++i) {
    sstr.str("");
    sstr << "histMJD_" << i;
    histName = sstr.str().c_str();
    histMJD[i] = new TH1D(histName, ";modified julian day;events",
    Int_t((mjd2 - mjd1) / (10. * second)), mjd1, mjd2);
  }

  // Track the local coordinates
  //vector< vector<Double_t> > LocCoord_theta(nMaps);
  //vector< vector<Double_t> > LocCoord_phi(nMaps);
  vector< vector<Float_t> > LocCoord_theta(nMaps);
  vector< vector<Float_t> > LocCoord_phi(nMaps);

  // Timers to figure out what's taking so long...
  TStopwatch timer1, timer2, timer3, timer4;
  timer1.Start();

  //*********************************************************************//
  // Begin iterating through events
  //*********************************************************************//
  int dayCounter = 0;
  int validCounter = 0;
  int mapIdx;
  bool temp;

  for (Long64_t jentry=0; jentry<nEntries; ++jentry) {

    cutDST->GetEntry(jentry);
    if (cfg.newConfig()) { 
      trigDST->GetEntry(jentry);
    }

    // Basic time check
    if (dst.ModJulDay < mjd1) {
      if (jentry % 10000000 == 0)
        cout << "Processed " << jentry << " entries before starting..." << endl;
      continue;
    }

    dayCounter += 1;
    if (dayCounter == 1) {
      cout << "First entry: " << jentry << endl;
      timer1.Stop();
      printf("Time to first entry: %7.3fs\n", timer.RealTime());
    }

    if (dayCounter % 1000000 == 0) {
      cout << "Processed " << dayCounter << " entries in the right day of " << jentry+1 << " total entries..." << endl;
      //printf("Time 2: %7.3fs\n", timer2.RealTime());
      //printf("Time 3: %7.3fs\n", timer3.RealTime());
      //printf("Time 4: %7.3fs\n", timer4.RealTime());
    }

    // Timer for all setup
    //timer2.Start(dayCounter == 1);

    // Additional checks on data
    isGood = true;
    mapIdx = 0;

    if (detector == "IC") {
      zenith = dst.LLHZenithDeg * deg2rad;
      azimuth = dst.LLHAzimuthDeg * deg2rad;
    }
    if (detector == "IT") {
      zenith = dst.Zenith;
      azimuth = dst.Azimuth;
    }

    // SimpleDST cuts (automatically included in Segev-processed files)
    // First: throw away reconstructions too close to poles
    //const Double_t zLo = 0.002;                  // 0.11 degrees
    //const Double_t zHi = TMath::Pi() - 0.002;    // 179.89 degrees
    const Float_t zLo = 0.002;                  // 0.11 degrees
    const Float_t zHi = TMath::Pi() - 0.002;    // 179.89 degrees
    if (zenith < zLo || zenith > zHi) {
      isGood = false;
    }
    // Second: require SMT08 trigger
    if (cfg.newConfig()) { 
      if (dst_trig.TriggID_1006 == 0) {
        isGood = false;
      }
    }

    // Reconstruction cuts
    if (!dst.isReco || zenith != zenith || azimuth != azimuth)
      isGood = false;

    // IceTop filter cut
    if (vm.count("filter")) {
      temp = filterCut(cfg,dst);
      if (not temp)
        isGood = false;
    }

    // Energy cuts for IceTop and IceCube
    if (detector == "IT" && vm.count("ebins"))
      mapIdx = ITenergyCut(dst, ebins);
    if (vm.count("spline"))
      mapIdx = ICenergyCut(dst, spline, zenith, ebins);
    if (vm.count("sbins"))
      mapIdx = ITs125Cut(dst, sbins);

    if (mapIdx == -1)
      isGood = false;

    //timer2.Stop();

    if (isGood && dst.ModJulDay <= (startMJD+dt)) {

      validCounter += 1;
      //timer3.Start(validCounter == 1);

      // Store local coordinates
      ++nEvents[mapIdx];
      LocCoord_theta[mapIdx].push_back(zenith);
      LocCoord_phi[mapIdx].push_back(azimuth);
      sphereDir.theta = zenith;
      sphereDir.phi = azimuth;
      pixelID = LocalMapInt[mapIdx].ang2pix(sphereDir);

      // Calculate solar dipole weighting
      double mjd = dst.ModJulDay;
      Direction dir(zenith,azimuth);
      Equatorial eq = GetEquatorialFromDirection(dir, mjd);
      double mjdTime = dst.ModJulDay;


      double eventweight = 1.0;
      if (sd) { 
        // if sd2 is true, apply 2nd-order correction
        eventweight = solar_dipole(dst.ModJulDay, eq.ra, eq.dec,sd2);
      }
      LocalMapInt[mapIdx][pixelID] += eventweight;

      // Calculate equatorial coordinates in other time frame
      double lst = GetGMST(mjd);
      double ra = eq.ra; 
      double dec = eq.dec;

      if (cfg.method == Config::antisid) { 
          double localAntiS = GetGMAST(mjd); 
          ra = fmod( eq.ra - (lst + localAntiS)*pi/12,2*pi);
      } else if (cfg.method == Config::extsid) { 
          double localExtS = GetGMEST(mjd); 
          ra = fmod( eq.ra - (lst + localExtS)*pi/12,2*pi);
      } else if (cfg.method == Config::solar) { 
          double tod = ( mjd - int(mjd) )* 24.; 
          ra = fmod(eq.ra - (lst + tod)*pi/12.,2*pi);
      }
      //timer3.Stop();
      //timer4.Start(validCounter == 1);

      // Write to map
      sphereDir.theta = pi/2. - dec;
      sphereDir.phi = ra;
      // Solar coordinates need a 180 deg flip in phi (definition difference)
      if (cfg.method == Config::solar) 
        sphereDir.phi -= pi;
      while (sphereDir.phi < 0)
        sphereDir.phi += 2.*pi;
      pixelID = DataMapInt[mapIdx].ang2pix(sphereDir);
      DataMapInt[mapIdx][pixelID] += eventweight;

      // Store time
      histMJD[mapIdx]->Fill(dst.ModJulDay);
      //timer4.Stop();
    }
  }

  // Can't stop early, because some files have events out of time order
  //if ((dst.ModJulDay > (startMJD + dt)) || (dst.ModJulDay > mjd2) ||
  //    (jentry + 1 == nEntries)) {
  for (unsigned mEntry = 0; mEntry<nMaps; mEntry++) {

    if (vm.count("ebins")) {
      cout << "Working on energy bin " << ebinstr[mEntry] << "-"
           << ebinstr[mEntry+1] << "GeV..." << endl;
    }
    if (vm.count("sbins")) {
      cout << "Working on s125 bin " << sbinstr[mEntry] << " to "
           << sbinstr[mEntry+1] << "s125..." << endl;
    }
    nUsedEvents[mEntry] += (nEvents[mEntry]);

    // Scramble the time
    cout << "  Scrambling time for (" << nBGResample << "x "
         << nEvents[mEntry] << " events)..." << endl;
    gRandom->SetSeed(0);

    for (Long64_t iEntry=0; iEntry<(Long64_t)(nEvents[mEntry]); iEntry++) {

      // Get local coordinates
      theta = LocCoord_theta[mEntry][iEntry];
      phi = LocCoord_phi[mEntry][iEntry];
      Direction dir(theta,phi);

      for (int k=0; k<nBGResample; ++k) {

        // Generate new equatorial coordinates
        rndMJD = histMJD[mEntry]->GetRandom();
        Equatorial eq = GetEquatorialFromDirection(dir, rndMJD);
        double new_ra = eq.ra; 
        double new_dec = eq.dec;

        double lst = GetGMST(rndMJD);
        if (cfg.method == Config::antisid) { 
            double localAntiS = GetGMAST(rndMJD);
            new_ra = fmod( eq.ra - (lst + localAntiS)*pi/12,2*pi);
        } else if (cfg.method == Config::extsid) { 
            double localExtS = GetGMEST(rndMJD);
            new_ra = fmod( eq.ra - (lst + localExtS)*pi/12,2*pi);
        } else if (cfg.method == Config::solar) {
            double tod = (rndMJD - int(rndMJD) )* 24.;
            new_ra = fmod(eq.ra - (lst + tod)*pi/12.,2*pi);
        }


        // Write to map
        sphereDir.theta = (pi/2. - new_dec);
        sphereDir.phi = new_ra;
        if (cfg.method == Config::solar) 
          sphereDir.phi -= pi;
        while (sphereDir.phi < 0)
          sphereDir.phi += 2.*pi;
        pixelID = DataMapInt[mEntry].ang2pix(sphereDir);

        BGMap[mEntry][pixelID] += 1.0;
      }
    }

    // Update the data map for this time interval
    for (int i=0; i<DataMap[mEntry].Npix(); ++i) {
      DataMap[mEntry][i] += DataMapInt[mEntry][i];
      LocalMap[mEntry][i] += LocalMapInt[mEntry][i];
    }
  }

  //cout << "jentry : " << jentry << endl;
  //jentry = jentry - 1;
  //startMJD += dt;

  //if (startMJD + second >= mjd2)
  //  break;
  //else {
  //  cout << "new startMJD :" << setprecision(12) << startMJD << endl;
  //  for (unsigned kEntry = 0; kEntry < nMaps; ++kEntry) {
  //    LocCoord_phi[kEntry].erase(LocCoord_phi[kEntry].begin(), 
  //        LocCoord_phi[kEntry].end());
  //    LocCoord_theta[kEntry].erase(LocCoord_theta[kEntry].begin(), 
  //        LocCoord_theta[kEntry].end());
  //    histMJD[kEntry]->Reset();
  //    DataMapInt[kEntry].fill(0);
  //    LocalMapInt[kEntry].fill(0);
  //    nEvents[kEntry] = 0;
  //  }
  //}
  //}

  for (unsigned m=0; m<nMaps; ++m) {

    // Scale background map
    for (int i=0; i<BGMap[m].Npix(); ++i)
      BGMap[m][i] *= alpha;

    cout << "Read " << nEntries << " events" << "\n"
         << "Used " << nUsedEvents[m] << " events" << endl;

    // Save BG, Data, and Local maps in one file
    arr<std::string> colname(3);
    colname[0] = "data map";
    colname[1] = "background map";
    colname[2] = "local map";

    sstr.str("");
    sstr << outBase;
    if (vm.count("ebins"))
      sstr << "_" << ebinstr[m] << "-" << ebinstr[m+1] << "GeV";
    if (vm.count("sbins"))
      sstr << "_" << sbinstr[m] << "to" << sbinstr[m+1] << "s125";
    sstr << "_" << yyyymmdd << ".fits";
    fitshandle fitsOut;
    fitsOut.create(sstr.str().c_str());

    fitsOut.add_comment("Maps: data, bg, local");
    //prepare_Healpix_fitsmap(fitsOut, DataMap[m], 
    //    FITSUTIL<float>::DTYPE, colname);
    //    FITSUTIL<double>::DTYPE, colname);
    // Temporary workaround - Planck Data Type for double is "9"
    prepare_Healpix_fitsmap(fitsOut, DataMap[m], PLANCK_FLOAT64, colname);
    fitsOut.write_column(1, DataMap[m].Map());
    fitsOut.write_column(2, BGMap[m].Map());
    fitsOut.write_column(3, LocalMap[m].Map());
    fitsOut.close();
  }

  // Clean up
  delete cutDST;
  if (cfg.newConfig()) { 
    delete trigDST;
  }
  for (unsigned m=0; m<nMaps; ++m)
    delete histMJD[m];

  timer.Stop();
  printf("RT=%7.3f s, Cpu=%7.3f s\n",timer.RealTime(),timer.CpuTime());

}




