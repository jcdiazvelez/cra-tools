#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <random>

#include <healpix_cxx/fitshandle.h>
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/healpix_map_fitsio.h>
#include <healpix_cxx/pointing.h>

#include <boost/program_options.hpp>
#include <boost/histogram.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/json.hpp>

#include <astro/astro.h>
#include <astro/time.h>
#include <Direction.h>
#include <units.h>

#include <omp.h>  // <-- for OpenMP

using namespace std;
namespace bh   = boost::histogram;
namespace io   = boost::iostreams;
namespace json = boost::json;
namespace po   = boost::program_options;

const double pi = constants::pi;


// Load json file into memory buffer
json::value load_json_file_gz(const string &filename) {

  ifstream file(filename, ios::binary | ios::ate);
  if (!file) throw runtime_error("Cannot open file: " + filename);

  const auto endpos = file.tellg();
  size_t compressed_size = 0;
  if (endpos != streampos(-1)) {
    compressed_size = static_cast<size_t>(static_cast<streamoff>(endpos));
  }
  file.seekg(0, ios::beg);

  io::filtering_istream in;
  if (filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz") {
    in.push(io::gzip_decompressor());
  }
  in.push(file);

  vector<char> buffer;
  if (compressed_size > 0) buffer.reserve(compressed_size);

  constexpr size_t chunk_size = 1 << 16;
  char temp[chunk_size];
  while (in.read(temp, chunk_size) || in.gcount() > 0) {
    buffer.insert(buffer.end(), temp, temp + in.gcount());
  }

  return json::parse({buffer.data(), buffer.size()});
}


// Event cuts
bool quality_cuts_passed(const json::object &event) {

  // NChannel cut
  int nchannel = event.at("nchannel").as_int64();
  if (nchannel < 10) {
      return false;
  }

  // Throw away reconstructions too close to the poles
  auto muE = event.at("PoleMuonLlhFitNanoDSTMuE").as_object();
  double zenith = muE.at("zenith").as_double();
  double azimuth = muE.at("azimuth").as_double();
  if (zenith < 0.002 || zenith > pi - 0.002) {
      return false;
  }

  // Reconstruction cuts
  // NOTE: this used to include an IsGoodLineFit or something similar. Discuss
  if (zenith != zenith || azimuth != azimuth) {
      return false;
  }
  return true;
}


// Fill time histogram
template <typename Histogram>
void fill_time_histogram(const json::array &data, Histogram& histMJD) {

    // Get the minimum and maximum time values from the histogram
    const auto& axis = histMJD.axis(0);
    double mjd0 = axis.bin(0).lower();
    double mjd1 = axis.bin(axis.size()-1).upper();

    // Iterate over events
    for (const auto &event_val : data) {

      // Access event information
      const json::object &event = event_val.as_object();

      // Apply time check -- could break when mjd > mjd1 if time-sorted
      double mjd = event.at("ModJulDate").as_double();
      if (mjd < mjd0 || mjd >= mjd1)
        continue;

      // Apply quality cuts
      bool qc_check = quality_cuts_passed(event);
      if (!qc_check)
        continue;

      // Fill time histogram
      histMJD(mjd);

    }
}


// Calculate equatorial coordinates and fill skymap
void fill_skymap(Healpix_Map<double>& skymap, double zenith, 
                 double azimuth, double mjd, string method) {

    // Calculate equatorial coordinates
    Direction dir(zenith,azimuth);
    Equatorial eq = GetEquatorialFromDirection(dir, mjd);
    double ra = eq.ra;
    double dec = eq.dec;

    // Additional calculations needed for other time frames
    double lst = GetGMST(mjd);
    if (method == "anti") {
      double localAntiS = GetGMAST(mjd);
      ra = fmod( eq.ra - (lst + localAntiS)*pi/12,2*pi);
    } 
    if (method == "ext") {
      double localExtS = GetGMEST(mjd);
      ra = fmod( eq.ra - (lst + localExtS)*pi/12,2*pi);
    } 
    if (method == "solar") {
      double tod = ( mjd - int(mjd) )* 24.;
      ra = fmod( eq.ra - (lst + tod)*pi/12.,2*pi);
    }

    // Finalize coordinates
    pointing sphere_dir(pi/2.-dec, ra);
    // Solar coordinates need a 180 deg flip in phi (definition difference)
    if (method == "solar")
      sphere_dir.phi -= pi;
    // Adjust negative values
    while (sphere_dir.phi < 0)
      sphere_dir.phi += 2.*pi;

    // Write to map
    int pixel_id = skymap.ang2pix(sphere_dir);
    skymap[pixel_id] += 1;

}


// -----------------------------------------------------------------------------
// This function fills local and data maps in parallel, using thread-local
// + reduce. It was written by Chat-GPT4 when asked to speed up the code,
// and it does so immensely (from ~30s to <1s/file). 
// -----------------------------------------------------------------------------

template <typename Histogram>
void time_scramble(const json::array &data,
                   const string &method,
                   Healpix_Map<double>& LocalMap,
                   Healpix_Map<double>& DataMap,
                   Healpix_Map<double>& RefMap,
                   const Histogram& histMJD) {

  // Time window from histogram
  const auto& axis = histMJD.axis(0);
  const double mjd0 = axis.bin(0).lower();
  const double mjd1 = axis.bin(axis.size() - 1).upper();

  // Bin weights and centers for "scrambled" times
  vector<double> weights;
  vector<double> bin_centers;
  for (auto x : bh::indexed(histMJD)) {
    weights.push_back(*x); // bin count
    bin_centers.push_back((x.bin().upper() + x.bin().lower())/2.);
  }

  // Random engine + discrete distribution
  mt19937 rng(random_device{}());
  discrete_distribution<> dist(weights.begin(), weights.end());

  // Thread-local maps (avoid writes to shared maps)
  const int nthreads = omp_get_max_threads();
  vector<Healpix_Map<double>> localsLocal(nthreads);
  vector<Healpix_Map<double>> localsData(nthreads);
  vector<Healpix_Map<double>> localsRef(nthreads);

  for (int i = 0; i < nthreads; ++i) {
    localsLocal[i].SetNside(LocalMap.Nside(), LocalMap.Scheme());
    localsLocal[i].fill(0.0);
    localsData[i].SetNside(DataMap.Nside(), DataMap.Scheme());
    localsData[i].fill(0.0);
    localsRef[i].SetNside(RefMap.Nside(), RefMap.Scheme());
    localsRef[i].fill(0.0);
  }

  // Initialize section to parallelize
  #pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    auto& locLocal = localsLocal[tid];
    auto& locData  = localsData[tid];
    auto& locRef  = localsRef[tid];

    // For loop to execute in parallel
    #pragma omp for schedule(static)
    for (size_t i = 0; i < data.size(); ++i) {
      const json::object &event = data[i].as_object();

      // Apply time check -- could adjust to break if events are time-sorted
      auto itMJD = event.find("ModJulDate");
      if (itMJD == event.end()) continue;
      const double mjd = itMJD->value().as_double();
      if (mjd < mjd0 || mjd >= mjd1)
        continue;

      // Apply quality cuts
      bool qc_check = quality_cuts_passed(event);
      if (!qc_check)
        continue;

      // Extract zenith and azimuth information
      auto itMuE = event.find("PoleMuonLlhFitNanoDSTMuE");
      if (itMuE == event.end() || !itMuE->value().is_object()) continue;
      const json::object &muE = itMuE->value().as_object();
      auto itZen = muE.find("zenith");
      auto itAzi = muE.find("azimuth");
      if (itZen == muE.end() || itAzi == muE.end()) continue;

      const double zenith  = itZen->value().as_double();
      const double azimuth = itAzi->value().as_double();

      // Fill local map (thread-local)
      pointing sphere_dir(zenith, azimuth);
      int pixel_id = locLocal.ang2pix(sphere_dir);
      locLocal[pixel_id] += 1.0;

      // Fill equatorial map (thread-local)
      fill_skymap(locData, zenith, azimuth, mjd, method);

      // Fill reference map (thread-local)
      for (int i = 0; i < 20; ++i) {
        int rnd_t_index = dist(rng);
        double rnd_mjd = bin_centers[rnd_t_index];
        fill_skymap(locRef, zenith, azimuth, rnd_mjd, method);
      }
    }
  }

  // Reduce: sum thread-local maps into the global maps
  for (int t = 0; t < nthreads; ++t) {
    for (int p = 0; p < LocalMap.Npix(); ++p) {
      LocalMap[p] += localsLocal[t][p];
      DataMap[p]  += localsData[t][p];
      RefMap[p]   += localsRef[t][p];
    }
  }
}




int main(int argc, char* argv[]) {

  po::options_description desc("Allowed options");
  desc.add_options()
      // Options used for all configurations
      ("help", "Produce help message")
      ("files", po::value<vector<string>>()->multitoken()->required(), "Files")
      ("yyyymmdd", po::value<string>(), "Desired date")
      ("method", po::value<string>(), "sid, solar, anti, ext")
      ("outBase", po::value<string>(), "Base name for outfile")
      //("spline", po::value<string>(), "File containing spline tables")
      //("ebins", po::value< vector<string> >()->multitoken(), "Energy bins")
  ;

  po::variables_map vm;

  // Disable short flags ('-') to allow use of negative signs in input
  po::store(po::parse_command_line(argc, argv, desc,
      po::command_line_style::unix_style ^ po::command_line_style::allow_short),
      vm);
  po::notify(vm);

  // Printout for help functionality
  if (vm.count("help")) {
    cout << desc << "\n";
    return 1;
  }

  // Establish desired starting time
  string yyyymmdd = vm["yyyymmdd"].as<string>();
  int yy = stoi(yyyymmdd.substr(0, 4));
  int mm = stoi(yyyymmdd.substr(4, 2));
  int dd = stoi(yyyymmdd.substr(6, 2));
  astro::Time t(yy, mm, dd, 0, 0, 0);
  double mjd0 = t.GetMJD();

  // Setup histogram for storing time information
  double mjd1 = mjd0 + 1.;
  int n_tbins = 8640;   // 10 second bins (86400 seconds/day)
  auto histMJD = bh::make_histogram(
      bh::axis::regular<>(n_tbins, mjd0, mjd1, "Times")
  );

  // Fill time histogram first - fast and avoids storing local coords in memory
  vector<string> files = vm["files"].as<vector<string>>();
  cout << "Now filling the time histogram..." << endl;

  // Time the process
  auto t0 = chrono::high_resolution_clock::now();

  // Note: is a try/catch statement here necessary?
  for (const auto &file : files) {
    try {
      // Load and parse JSON file
      json::value root = load_json_file_gz(file);
      const json::array &data = root.at("data").as_array();

      // Fill time histogram from data
      fill_time_histogram(data, histMJD);
    } 
    catch (exception &e) {
      cerr << "Error: " << e.what() << "\n";
    }
  }

  // Timing update
  auto t1 = chrono::high_resolution_clock::now();
  chrono::duration<double> dt = t1 - t0;
  cout << "Time to fill time histogram: " << dt.count() << " s\n";

  /*
  // Print histogram
  int count = 0;
  for (auto&& x : bh::indexed(histMJD)) {
    if (count++ == 10) break;
    cout << "[" << x.bin(0).lower() << ", " << x.bin(0).upper() << "): "
              << *x << "\n";
  }
  */

  // Initialize skymaps
  int nside = 64;
  Healpix_Map<double> LocalMap, DataMap, RefMap;
  LocalMap.SetNside(nside, RING);
  LocalMap.fill(0.);
  DataMap.SetNside(nside, RING);
  DataMap.fill(0.);
  RefMap.SetNside(nside, RING);
  RefMap.fill(0.);

  // Time the process
  t0 = chrono::high_resolution_clock::now();

  // Fill skymaps
  string method = vm["method"].as<string>();
  cout << "Calculating and storing equatorial coordinates..." << endl;
  for (const auto &file : files) {
    try {
      // Load and parse JSON file
      json::value root = load_json_file_gz(file);
      const json::array &data = root.at("data").as_array();

      // Calculate equatorial coordinates and fill skymaps
      time_scramble(data, method, LocalMap, DataMap, RefMap, histMJD);
    } 
    catch (exception &e) {
      cerr << "Error: " << e.what() << "\n";
    } 
  }

  // Timing update
  t1 = chrono::high_resolution_clock::now();
  dt = t1 - t0;
  cout << "Time to fill skymaps: " << dt.count() << " s\n";

  // Closing up shop: scale the background map
  for (int i=0; i<RefMap.Npix(); ++i)
    RefMap[i] *= 1/20.;

  // Save Reference, Data, and Local maps in one file
  arr<std::string> colname(3);
  colname[0] = "data map";
  colname[1] = "background map";
  colname[2] = "local map";
  string outfile = vm["outBase"].as<string>();
  fitshandle fits_out;
  fits_out.create(outfile.c_str());
  fits_out.add_comment("Maps: data, reference, local");
  prepare_Healpix_fitsmap(fits_out, DataMap, PLANCK_FLOAT64, colname);
  fits_out.write_column(1, DataMap.Map());
  fits_out.write_column(2, RefMap.Map());
  fits_out.write_column(3, LocalMap.Map());
  fits_out.close();

  return 0;

}
