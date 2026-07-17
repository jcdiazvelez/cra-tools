#include <SimpleDST.h>
#include <config.h>

#include <TChain.h>
#include <TH1D.h>
#include <TRandom.h>

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include <algorithm>

#include <healpix_cxx/fitshandle.h>
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/healpix_map_fitsio.h>
#include <healpix_cxx/pointing.h>

#include <photospline/splinetable.h>
#include <photospline/bspline.h>

#include <boost/program_options.hpp>
#include <boost/make_shared.hpp>
#include <boost/filesystem.hpp>

#include <astro/astro.h>
#include <astro/time.h>
#include <Direction.h>
#include <solardipole.h>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

const double hour = 1 / 24.;
const double second = hour / 3600.;

typedef Healpix_Map<float> SkyMap;
typedef boost::shared_ptr<SkyMap> SkyMapPtr;


// Default values for inputs
struct ProgramOptions {

  // Input/output
  std::vector< std::string > input;
  std::string outdir = "./sample";
  std::string outfile = "CR_TS_sample";

  // Time-scrambling inputs
  std::string config = "IC86";
  std::string yyyymmdd;
  int nInt;
  std::string method = "sid";
  int resample = 20;

  // Energy reconstruction
  std::string spline;
  std::vector<float> ebins;
  std::vector<float> sbins;

  // Solar dipole flags
  bool sd = false;
  bool sd2 = false;

};


// NOTES:
// - suggestions from chatgpt/Claude include:
//    - avoid reconstructing Direction if it's not lightweight
//    - 

int IT_s125_bin(const SimpleDST& dst, const std::vector<float>& sbins);
int IC_energy_bin(const SimpleDST& dst, photospline::splinetable<> &table, double zenith, const std::vector<float>& ebins);
double time_standard_ra(double mjd, const Equatorial& eq, const Config& cfg);
SkyMapPtr MakeSkyMap(int nside);


int main(int argc, char* argv[]) {

  //=====================================================================//
  // Input options
  //=====================================================================//

  ProgramOptions opts;

  po::options_description desc("Allowed options");
  desc.add_options()

      // Options used for all configurations
      ("help", "Produce help message")
      ("input", po::value<std::vector< std::string > >(&opts.input)
          ->multitoken()
          ->required(),
          "Input files")
      ("outdir", po::value<std::string>(&opts.outdir)
          ->default_value(opts.outdir),
          "Directory of output")
      ("outfile", po::value<std::string>(&opts.outfile)
          ->default_value(opts.outfile),
          "Base name for output file")
      ("config", po::value<std::string>(&opts.config)
          ->default_value(opts.config),
          "Detector configuration. Options: IC86 | IT81 (pre-11-year) | "
          "ITpass2 (11-year+)")
      ("method", po::value<std::string>(&opts.method)
          ->default_value(opts.method),
          "Time standard (sid|anti|solar|ext)")
      ("yyyymmdd", po::value<std::string>(&opts.yyyymmdd),
          "Restrict scrambling to times within target date")

      // Less-common options
      ("nInt", po::value<int>(&opts.nInt),
          "Integration time in hours. Defaults to length of input files OR "
          "24h if specifying a calendar date")
      ("resample", po::value<int>(&opts.resample)
          ->default_value(opts.resample),
          "Number of times to resample background values")

      // Solar dipole flags
      ("sd", po::value<bool>(&opts.sd)
          ->default_value(opts.sd),
          "Correct for solar dipole for each event")
      ("sd2", po::value<bool>(&opts.sd2)
          ->default_value(opts.sd2),
          "Apply 2nd-order dipole correction")

      // IceCube specific options
      ("spline", po::value<std::string>(&opts.spline),
          "File containing spline tables")
      ("ebins", po::value< std::vector<float> >(&opts.ebins)
          ->multitoken(),
          "Energy bins")

      // Icetop specific options
      ("sbins", po::value< std::vector<float> >(&opts.sbins)
          ->multitoken(),
          "S125 bins")
  ;

  po::variables_map vm;

  // Disable short flags ('-') to allow use of negative signs in input
  po::store(po::parse_command_line(argc, argv, desc,
      po::command_line_style::unix_style ^ po::command_line_style::allow_short),
      vm);

  po::notify(vm);

  // Help functionality
  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  // Print input filelist
  std::cout << "Input files are:\n";
  for (unsigned i = 0; i < opts.input.size(); ++i)
      std::cout << opts.input[i] << "\n";


  //=====================================================================//
  // Setup for binning in energy or S125
  //=====================================================================//

  // Default setup: one map
  unsigned nMaps = 1;

  // Energy and S125
  const bool useEnergyBins = ( vm.count("ebins") && vm.count("spline") );
  const bool useS125Bins = vm.count("sbins");

  if (useEnergyBins && useS125Bins)
      throw std::runtime_error("Critical error: Can't have S125 and Ebins!");

  // Additional behavior if binning in energy
  photospline::splinetable<> spline;
  if (useEnergyBins) {

    // Minimum energy bin size
    if (opts.ebins.size() < 2)
      throw std::runtime_error("Critical error: Must be >=2 energy bins");

    std::cout << "Ebin values:" << "\n";
    for (float val : opts.ebins)
        std::cout << val << " ";
    std::cout << "\n";

    // Update number of maps
    nMaps = opts.ebins.size() - 1;

    // Read in spline file for energy lookup
    std::cout << "Reading spline file: " << opts.spline << "\n";
    spline.read_fits(opts.spline.c_str());
  }

  // Additional behavior if binning in S125
  if (useS125Bins) {

    // Minimum energy bin size
    if (opts.sbins.size() < 2)
      throw std::runtime_error("Critical error: Must be >=2 S125 bins");

    std::cout << "Sbin values:" << "\n";
    for (float val : opts.sbins)
        std::cout << val << " ";
    std::cout << "\n";

    // Update number of maps
    nMaps = opts.sbins.size() - 1;
  }

  std::cout << "Number of maps: " << nMaps << "\n";


  //=====================================================================//
  // Load input files into chain
  //=====================================================================//

  Config cfg(vm);
  const char* masterTree;

  // Get naming for tree from detector configuration
  if (cfg.detector == Config::IceCube)
      masterTree = "CutDST";
  // IceTop tree differently named starting with eleven-year analysis
  else if (cfg.cfg == Config::ITv3)
      masterTree = "MasterTree";
  // Name for earlier IceTop trees
  else
      masterTree = "master_tree";

  // Initialize the chain and read data
  auto cutDST = std::make_unique<TChain>(masterTree);
  for (unsigned i = 0; i < opts.input.size(); ++i) {
    cutDST->Add(opts.input[i].c_str());
  }

  std::string detector_string = (cfg.detector==Config::IceCube) ? "IC" : "IT";
  SimpleDST dst(cutDST.get(), detector_string);
  std::cout << "Number of chained files: " << cutDST->GetNtrees() << "\n";

  const auto nEntries = cutDST->GetEntries();
  std::cout << "Number of entries: " << nEntries << "\n";


  //=====================================================================//
  // Time-scrambling window
  //=====================================================================//

  if ( vm.count("nInt") && !vm.count("yyyymmdd") ) {
    throw std::runtime_error("Critical error: Specifying a time-scrambling "
                             "window means program will run for 24h "
                             "(need to use yyyymmdd argument too");
  }

  // Use fixed start and stop times if yyyymmdd option provided
  double start_mjd, stop_mjd;
  if (vm.count("yyyymmdd")) {
    int yy = std::stoi(opts.yyyymmdd.substr(0, 4).c_str());
    int mm = std::stoi(opts.yyyymmdd.substr(4, 2).c_str());
    int dd = std::stoi(opts.yyyymmdd.substr(6, 2).c_str());
    astro::Time t(yy, mm, dd, 0, 0, 0);
    start_mjd = t.GetMJD();
    stop_mjd = start_mjd + 1;
  }

  // Otherwise, set start and stop time based on data
  else {
    cutDST->GetEntry(0);
    start_mjd = dst.ModJulDay;
    cutDST->GetEntry(nEntries - 1);
    stop_mjd = dst.ModJulDay;
  }

  // MJD1 will update to a new start time if scrambling window is passed
  double mjd1 = start_mjd;

  // Time-scrambling window defaults to input value
  int dt_hrs;
  if (vm.count("nInt"))
      dt_hrs = opts.nInt;
  // Alt: 24h scrambling if a calendar day is specified
  else if (vm.count("yyyymmdd"))
      dt_hrs = 24;
  // Allow for >24h scrambling for everything else
  else
      dt_hrs = static_cast<int>(std::ceil((stop_mjd - start_mjd) / hour));

  // Integration time
  const double dt = dt_hrs * hour;
  std::cout << "Integration time = " << dt_hrs << " (hrs) " << dt << " (day)\n";
  std::cout << "Reading " << nEntries << " entries...\n";


  //=====================================================================//
  // Data storage and variable definitions
  //=====================================================================//

  // Skymaps
  std::vector<SkyMapPtr> LocalMap;
  std::vector<SkyMapPtr> DataMap;
  std::vector<SkyMapPtr> BGMap;

  const int NSide = 64;
  for (unsigned int i=0; i<nMaps; ++i) {
    LocalMap.push_back(MakeSkyMap(NSide));
    DataMap.push_back(MakeSkyMap(NSide));
    BGMap.push_back(MakeSkyMap(NSide));
  }

  // Histograms for storing time information
  gRandom->SetSeed(0);
  std::vector< std::unique_ptr<TH1D> > histMJD;
  for (unsigned i = 0; i < nMaps; ++i) {
    std::string histName = "histMJD_" + std::to_string(i);
    histMJD.emplace_back(std::make_unique<TH1D>(histName.c_str(),
        ";modified julian day;events",
        Int_t((stop_mjd - start_mjd) / (10. * second)), start_mjd, stop_mjd));
  }

  // Storage for local coordinates
  std::vector< std::vector<Float_t> > LocCoord_theta(nMaps);
  std::vector< std::vector<Float_t> > LocCoord_phi(nMaps);
  // Save some reallocation time if not binning maps
  if (nMaps == 1) {
    int n_windows = std::max(24/dt_hrs, 1);
    int expected_events = nEntries / n_windows;
    LocCoord_theta[0].reserve(expected_events);
    LocCoord_phi[0].reserve(expected_events);
  }

  // General counters for number of events
  std::vector<Long64_t> nEvents(nMaps, 0);
  std::vector<Long64_t> nUsedEvents(nMaps, 0);
  Long64_t dayCounter = 0;

  // Rescaling factor from number of resamples
  const double alpha = 1. / opts.resample;

  // Start timer
  auto t_start = std::chrono::steady_clock::now();

  //=====================================================================//
  // Begin iterating through events
  //=====================================================================//

  Long64_t nevent = 0;
  double mjd = 0;

  // Read all events (potentially up to stop time)
  while ( nevent < nEntries && mjd < stop_mjd ) {

    cutDST->GetEntry(nevent);

    //==================================================================//
    // Basic tracking output
    //==================================================================//

    // Events before start of time window
    mjd = dst.ModJulDay;
    if (mjd < start_mjd) {
      if (nevent % 10000000 == 0)
        std::cout << "Processed " << nevent << " entries before starting...\n";
      nevent++;
      continue;
    }

    // Time to hit first entry
    dayCounter += 1;
    if (dayCounter == 1) {
      std::cout << "First entry: " << nevent << "\n";
      auto t1 = std::chrono::steady_clock::now();
      std::chrono::duration<double> t_to_start = t1 - t_start;
      std::cout << "Time to first entry: " << std::fixed
                << std::setprecision(3) << t_to_start.count() << "s\n";
    }

    // Counted event tracker
    if (dayCounter % 1000000 == 0) {
      std::cout << "Processed " << dayCounter << " entries in the right day of " << nevent+1 << " total entries...\n";
    }

    //==================================================================//
    // Event cuts
    //==================================================================//

    bool event_passed = true;

    // Extract zenith, azimuth, and fit_status information
    double zenith, azimuth;
    bool fitPassed;
    if (cfg.detector == Config::IceCube) {
      zenith = dst.LLHZenithDeg * M_PI/180.;
      azimuth = dst.LLHAzimuthDeg * M_PI/180.;
      fitPassed = dst.isReco;
    }
    else {
      zenith = dst.SPZenith;
      azimuth = dst.SPAzimuth;
      fitPassed = (dst.SPFitStatus == 0);
    }

    // Throw away reconstructions too close to poles
    const float zLo = 0.002;           // 0.11 degrees
    const float zHi = M_PI - 0.002;    // 179.89 degrees
    if (zenith < zLo || zenith > zHi)
      event_passed = false;

    // Reconstruction cuts
    if (!fitPassed || std::isnan(zenith) || std::isnan(azimuth))
      event_passed = false;

    // Energy cuts for IceTop and IceCube
    int mapIdx = 0;
    if (useEnergyBins)
      mapIdx = IC_energy_bin(dst, spline, zenith, opts.ebins);
    if (useS125Bins)
      mapIdx = IT_s125_bin(dst, opts.sbins);
    // Energy bin of -1 is outside range
    if (mapIdx == -1)
      event_passed = false;

    //==================================================================//
    // Store local coords and event time
    //==================================================================//

    if ( event_passed && mjd <= (mjd1+dt) ) {

      // Store local coordinates
      ++nEvents[mapIdx];
      LocCoord_theta[mapIdx].push_back(zenith);
      LocCoord_phi[mapIdx].push_back(azimuth);

      // Calculate coordinates in equatoral time
      Direction dir(zenith,azimuth);
      Equatorial eq = GetEquatorialFromDirection(dir, mjd);

      // RA can change depending on time standard (anti, ext, solar, sid)
      double ra = time_standard_ra(mjd, eq, cfg);

      // Calculate solar dipole weighting
      double eventweight = 1.0;
      if (opts.sd)
        eventweight = solar_dipole(mjd, eq.ra, eq.dec, opts.sd2);
      pointing localDir(zenith, azimuth);
      SkyMap& tmp_local = *LocalMap[mapIdx];
      int pixelID = tmp_local.ang2pix(localDir);
      tmp_local[pixelID] += eventweight;

      // Write to map
      pointing eqDir(M_PI/2.-eq.dec, ra);
      SkyMap& tmp_data = *DataMap[mapIdx];
      pixelID = tmp_data.ang2pix(eqDir);
      tmp_data[pixelID] += eventweight;

      // Store time
      histMJD[mapIdx]->Fill(mjd);
    }

    //==================================================================//
    // Time Scrambling
    //==================================================================//

    if (mjd > (mjd1+dt) || mjd > stop_mjd || nevent+1 == nEntries) {

      for (unsigned mEntry = 0; mEntry<nMaps; ++mEntry) {

        if (useEnergyBins) {
          std::cout << "Working on energy bin " << opts.ebins[mEntry] << "-"
               << opts.ebins[mEntry+1] << "GeV...\n";
        }
        if (useS125Bins) {
          std::cout << "Working on s125 bin " << opts.sbins[mEntry] << " to "
               << opts.sbins[mEntry+1] << "s125...\n";
        }
        nUsedEvents[mEntry] += (nEvents[mEntry]);

        // Scramble the time
        std::cout << "  Scrambling time for (" << opts.resample << " x "
             << nEvents[mEntry] << " events)...\n";

        for (Long64_t i = 0; i < nEvents[mEntry]; ++i) {

          // Get local coordinates
          double theta = LocCoord_theta[mEntry][i];
          double phi = LocCoord_phi[mEntry][i];
          Direction dir(theta,phi);

          for (int k=0; k<opts.resample; ++k) {

            // Generate new equatorial coordinates
            double rndMJD = histMJD[mEntry]->GetRandom();
            Equatorial eq = GetEquatorialFromDirection(dir, rndMJD);

            // Calculate solar dipole weighting
            double eventweight = 1.0;
            if (opts.sd)
              eventweight = solar_dipole(rndMJD, eq.ra, eq.dec, opts.sd2);

            // RA can change depending on time standard (anti, ext, solar, sid)
            double tmp_ra = time_standard_ra(rndMJD, eq, cfg);

            // Write to map
            pointing eqDir(M_PI/2.-eq.dec, tmp_ra);
            SkyMap& tmp_bg = *BGMap[mEntry];
            int pixelID = tmp_bg.ang2pix(eqDir);
            tmp_bg[pixelID] += eventweight * alpha;
          }
        }
      }

      // If not on the last entry, scrambling triggered by event outside dt
      if (nevent + 1 != nEntries) {

        // Update beginning of time-scrambling window
        mjd1 += dt;
        std::cout << "new start_mjd :" << std::setprecision(12) << mjd1 << "\n";

        // Clear storage
        for (unsigned kEntry = 0; kEntry < nMaps; ++kEntry) {
          LocCoord_phi[kEntry].clear();
          LocCoord_theta[kEntry].clear();
          histMJD[kEntry]->Reset();
          nEvents[kEntry] = 0;
        }

        // Return to loop without incrementing counter to revisit event
        continue;

      }
    }

    nevent += 1;

  }


  //=====================================================================//
  // Finish up
  //=====================================================================//

  for (unsigned m=0; m<nMaps; ++m) {

    const SkyMap& data = *DataMap[m];
    const SkyMap& bg = *BGMap[m];
    const SkyMap& local = *LocalMap[m];

    std::cout << "Read " << nEntries << " events\n"
         << "Used " << nUsedEvents[m] << " events\n";

    // Save BG, Data, and Local maps in one file
    arr<std::string> colname(3);
    colname[0] = "data map";
    colname[1] = "background map";
    colname[2] = "local map";

    std::stringstream namefits;

    namefits << opts.outdir << "/" << opts.outfile;
    if (useEnergyBins)
      namefits << "_" << opts.ebins[m] << "-" << opts.ebins[m+1] << "GeV";
    if (useS125Bins)
      namefits << "_" << opts.sbins[m] << "to" << opts.sbins[m+1] << "s125";
    namefits << ".fits.gz";

    // Create output directory if it does not exist yet
    fs::path out_directory(opts.outdir);
    if (!(fs::exists(out_directory))) {
      std::cout << "Directory " << opts.outdir << " doesn't exist\n";
      if (fs::create_directory(out_directory))
          std::cout << "....successfully created!\n";
    }

    // Overwrite file if it already exists
    if (fs::exists(namefits.str()))
         fs::remove(namefits.str());

    std::cout << "Writing output to " << namefits.str() << "\n";
    fitshandle fitsOut;
    fitsOut.create(namefits.str().c_str());

    fitsOut.add_comment("Maps: data, bg, local");
    prepare_Healpix_fitsmap(fitsOut, data, PLANCK_FLOAT64, colname);
    fitsOut.write_column(1, data.Map());
    fitsOut.write_column(2, bg.Map());
    fitsOut.write_column(3, local.Map());
    fitsOut.close();
  }

  auto t_stop = std::chrono::steady_clock::now();
  std::chrono::duration<double> t_to_end = t_stop - t_start;
  std::cout << "Total duration: " << std::fixed << std::setprecision(3) 
            << t_to_end.count() << "s\n";

  return 0;
}



//=====================================================================//
// Helper functions
//=====================================================================//

// Automatic creation of empty skymap pointers
SkyMapPtr MakeSkyMap(int nside) {

  auto map = boost::make_shared<SkyMap>();

  map->SetNside(nside, RING);
  map->fill(0.0);

  return map;
}


// Adjust RA based on selected time standard (sid|solar|anti|ext)
double time_standard_ra(double mjd, 
                        const Equatorial& eq, 
                        const Config& cfg) {

  // Default behavior assumes sidereal time
  double ra = eq.ra;

  // Shift RA for chosen time standard
  double lst = GetGMST(mjd);
  if (cfg.method == Config::antisid) {
    double localAntiS = GetGMAST(mjd);
    ra = fmod( eq.ra - (lst + localAntiS) * M_PI/12, 2*M_PI);
  }

  else if (cfg.method == Config::extsid) {
    double localExtS = GetGMEST(mjd);
    ra = fmod( eq.ra - (lst + localExtS) * M_PI/12, 2*M_PI);
  }

  else if (cfg.method == Config::solar) {
    double tod = ( mjd - int(mjd) ) * 24.;
    ra = fmod(eq.ra - (lst + tod) * M_PI/12., 2*M_PI);
    // Solar coordinates need an additional 180-deg flip (definition)
    ra -= M_PI;
  }

  // Catch anything that went negative
  while (ra < 0)
    ra += 2.*M_PI;

  return ra;

}


// Bin event in energy
int IC_energy_bin(const SimpleDST& dst, 
                  photospline::splinetable<> &spline, 
                  double zenith, 
                  const std::vector<float>& ebins) {

  // Setup basic parameters
  double x = cos(zenith);
  double y = log10(dst.NChannels);

  // Boundary check (energy cut tables go to 0.3 in cos(zenith))
  if (x < 0.3)
    return -1;

  // Catch additional outliers
  double coords[2] = {x, y};
  int centers[spline.get_ndim()];
  if (!spline.searchcenters(coords, centers)) {
    std::cout << "Variables outside of table boundaries\n";
    std::cout << "x: " << x << " y: " << y << "\n";
    return -1;
  }

  // Calculate reconstructed energy
  double median = spline.ndsplineeval(coords, centers, 0);
  // Make sure we're in the energy bin range
  if ((median < ebins[0]) || (median > ebins.back()))
    return -1;

  // Get energy bin
  int ebin = 0;
  while (median > ebins[ebin+1])
    ebin += 1;

  return ebin;
}


// Bin event in S125
int IT_s125_bin(const SimpleDST& dst, const std::vector<float>& sbins) {

  // Get desired s125 value
  double s125 = (dst.nStations >= 5) ? dst.s125 : dst.ss125;
  double logS125 = log10(s125);

  // Make sure we're in the bin range
  if ((logS125 < sbins[0]) || (logS125 > sbins.back()))
    return -1;

  // Get s125 bin
  int sbin = 0;
  while (logS125 > sbins[sbin+1])
    sbin += 1;

  return sbin;
}

