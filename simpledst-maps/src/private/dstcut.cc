
#include <astro/astro.h>
#include <TDST.h>
#include <Direction.h>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TSystem.h>
#include <TTree.h>

#include <bitset>
#include <iomanip>
#include <memory>
#include <sstream>
#include <iostream>
#include <string>


// -----------------------------------------------------------------------------
// Cuts: zenith angle, trigger tag, NaN removal, etc.
// -----------------------------------------------------------------------------
Bool_t
cut_zenith_poles(const TDST* dst, Bool_t cutOnLineFit = kFALSE)
{
  // Rasha's cuts on poles
//  const Direction* dirLin = &(dst->reco1.direction);
  const Direction* dirLLH = &(dst->reco2.direction);
  const Double_t zLo = 0.002;                  // 0.11 degrees
  const Double_t zHi = TMath::Pi() - 0.002;    // 179.89 degrees

//  cout << dst->reco2.direction.GetZenith() << endl;

//  if (cutOnLineFit) 
//    return dirLin->GetZenith() > zLo && dirLin->GetZenith() < zHi &&
//           dirLLH->GetZenith() > zLo && dirLLH->GetZenith() < zHi; 
//  else
    return dirLLH->GetZenith() > zLo && dirLLH->GetZenith() < zHi; 
}

Bool_t
cut_nan(const TDST* dst, Bool_t cutOnLineFit = kFALSE)
{
  // Rasha's cuts NaNs and other garbage input
  const Direction* dirLin = &(dst->reco1.direction);
  Double_t a1 = dirLin->GetAzimuth();
  Double_t z1 = dirLin->GetZenith();

  const Direction* dirLLH = &(dst->reco2.direction);
  Double_t a2 = dirLLH->GetAzimuth();
  Double_t z2 = dirLLH->GetZenith();

  if (cutOnLineFit)
    return (z1 == z1 && a1 == a1 && z2 == z2 && a2 == a2 &&
            dst->time == dst->time &&
            dst->ndir == dst->ndir &&
            dst->nstring == dst->nstring &&
            dst->nchan == dst->nchan);
  else
    return (z2 == z2 && a2 == a2 &&
            dst->time == dst->time &&
            dst->ndir == dst->ndir &&
            dst->nstring == dst->nstring &&
            dst->nchan == dst->nchan);
}

Bool_t
cut_smt_trigger(const UInt_t TriggID_1006)
{
  // See wiki.icecube.wisc.edu/index.php/TDST#A:_Trigger_tag_Information
  // ID = 2: simple majority
  return TriggID_1006 != 0;
}


Bool_t
cut_reco_idtag(const TDST* dst)
{
  // Make sure the LLH reconstruction is used, not the TOI reconstruction
  return (dst->reco2.reco_id == 2);
}

Bool_t
accept_dst(const TDST* dst)
{
  Double_t zenith1 = dst->reco1.direction.GetZenith();
  Double_t azimuth1 = dst->reco1.direction.GetAzimuth();
  Double_t zenith2 = dst->reco2.direction.GetZenith();
  Double_t azimuth2 = dst->reco2.direction.GetAzimuth();
  UShort_t reco2id = dst->reco2.reco_id;
  UInt_t trigtag = dst->triggertag;
  Double_t time = dst->time;
  UInt_t ndir = dst->ndir;
  UInt_t nstring = dst->nstring;
  UShort_t ndoms = dst->nchan;

  return (zenith2==zenith2) && (zenith2>0.02) && (zenith2<TMath::Pi()-0.02) &&
         //(zenith1==zenith1) && (reco2id==2.0) && ((trigtag & 2) != 0) &&
         (zenith1==zenith1) && (reco2id==2.0) &&
         (azimuth1==azimuth1) && (azimuth2==azimuth2) && 
         (time==time) && (ndir==ndir) && (nstring==nstring) && (ndoms==ndoms);
}

// -----------------------------------------------------------------------------
// Calculate the opening angle between the line fit and log-likelihood fit
// -----------------------------------------------------------------------------
Double_t
linllh_opening_angle(const TDST* dst)
{
  const Direction* dirLin = &(dst->reco1.direction);
  const Direction* dirLLH = &(dst->reco2.direction);

  const Double_t a1 = dirLin->GetAzimuth();
  const Double_t z1 = dirLin->GetZenith();

  const Double_t a2 = dirLLH->GetAzimuth();
  const Double_t z2 = dirLLH->GetZenith();

  const Double_t dp = TMath::Sin(z1)*TMath::Sin(z2)*TMath::Cos(a1 - a2) +
                      TMath::Cos(z1)*TMath::Cos(z2);

  return TMath::ACos(dp);
}

// -----------------------------------------------------------------------------
// Extract the subrun Id from a DST file name
// -----------------------------------------------------------------------------
UShort_t
get_subrun_id(const Char_t* filename)
{
  string name(filename);
  size_t i = name.find_last_of("_") + 1;
  size_t j = name.find(".root");

  stringstream sstr;
  UShort_t id;
  sstr << name.substr(i, j - i) << endl;
  sstr >> id;

  return id;
}

int main(int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [output.root] [infiles...]\n\n";
    return 1;
  }

  // Set the maximum file size to be 256G
  Long64_t maxTreeSize = 1;
  maxTreeSize <<= 38;
  TTree::SetMaxTreeSize(maxTreeSize);

  // Set up the output tree
  TFile* output = new TFile(argv[1], "RECREATE");
  TTree* cutDST = new TTree("CutDST", "IC79 DST + zenith + NaN + trigger cuts");

  Double_t mjd, usec, mjdTime;
  cutDST->Branch("ModJulDay", &mjdTime, "ModJulDay/D");

  Double_t secsInDay, nsInDay, localMST;
  cutDST->Branch("LocalMST", &localMST, "LocalMST/D");

  Float_t llhAzimuth;
  cutDST->Branch("LLHAzimuthDeg", &llhAzimuth, "LLHAzimuthDeg/F");

  Float_t llhZenith;
  cutDST->Branch("LLHZenithDeg", &llhZenith, "LLHZenithDeg/F");

  Float_t linllhOpeningAngle;
  cutDST->Branch("LinLLHOpeningAngleDeg", &linllhOpeningAngle, "LinLLHOpeningAngleDeg/F");

  // RA, Dec in current epoch
  Float_t RA;
  cutDST->Branch("RADeg", &RA, "RADeg/F");

  Float_t Dec;
  cutDST->Branch("DecDeg", &Dec, "DecDeg/F");

  // RA, Dec calculated with anti-sidereal time
  Float_t RAAntiS;
  cutDST->Branch("RAAntiS", &RAAntiS, "RAAntiS/F");

  Float_t DecAntiS;
  cutDST->Branch("DecAntiS", &DecAntiS, "DecAntiS/F");

  // RA, Dec calculated with solar time
  Float_t RASolar;
  cutDST->Branch("RASolar", &RASolar, "RASolar/F");

  Float_t DecSolar;
  cutDST->Branch("DecSolar", &DecSolar, "DecSolar/F");

  // Position of the sun
  Float_t RASun;
  cutDST->Branch("RASun", &RASun, "RASun/F");

  Float_t DecSun;
  cutDST->Branch("DecSun", &DecSun, "DecSun/F");

  // Position of the moon
  Float_t RAMoon;
  cutDST->Branch("RAMoon", &RAMoon, "RAMoon/F");

  Float_t DecMoon;
  cutDST->Branch("DecMoon", &DecMoon, "DecMoon/F");

  // Muon energy
  Float_t logMuE;
  cutDST->Branch("LogMuE", &logMuE, "LogMuE/F");

  // Other geometrical quantities: impact parameter, center of gravity, etc.
  Float_t rlogl;
  cutDST->Branch("RLogL", &rlogl, "RLogL/F");

  Float_t sdcog;
  cutDST->Branch("SDCoG", &sdcog, "SDCoG/F");

  Float_t cogx;
  cutDST->Branch("CoG_X", &cogx, "CoG_X/F");

  Float_t cogy;
  cutDST->Branch("CoG_Y", &cogy, "CoG_Y/F");

  Float_t cogz;
  cutDST->Branch("CoG_Z", &cogz, "CoG_Z/F");

  UInt_t ldir;
  cutDST->Branch("LDir", &ldir, "LDir/i");

  // Run information
  UInt_t runId;
  cutDST->Branch("RunId", &runId, "RunId/i");

  // Number of hits and strings
  UInt_t ndir;
  cutDST->Branch("NDirHits", &ndir, "NDirHits/i");

  UShort_t nchan;
  cutDST->Branch("NChannels", &nchan, "NChannels/s");

  UShort_t nstring;
  cutDST->Branch("NStrings", &nstring, "NStrings/s");

  // Sub-run information
  UShort_t subrunId;
  cutDST->Branch("SubRunId", &subrunId, "SubRunId/s");

  Bool_t isGoodLineFit;
  cutDST->Branch("IsGoodLineFit", &isGoodLineFit, "IsGoodLineFit/O");

  // Loop through all input files
  cout << "Reading through " << argc - 2 << " DST files..." << endl;
  for (int i = 2; i < argc; ++i) {
    cout << "File " << i - 1 << ":" << endl;
    string currentFile = argv[i];
    string badFile = "Run00126158_Subrun00000156";
    size_t found = currentFile.find(badFile);
    if (found != string::npos) {
      cout << "Bad file skipped" << endl;
      continue;
    }

    // Set up the input chain
    auto_ptr<TChain> chain = auto_ptr<TChain>(new TChain("tree"));
    chain->Add(argv[i]);

    TDST* dst = NULL;
    chain->SetBranchAddress("dst", &dst);

    UInt_t TriggID_1006;
    chain->SetBranchAddress("TriggID_1006", &TriggID_1006);

    subrunId = get_subrun_id(argv[i]);

    // Loop through the input chain
    Int_t nEvents = 0;
    const Int_t step = chain->GetEntries() / 79;

    cout << "Entries: " << chain->GetEntries() << endl;

    for (UInt_t j = 0; j < chain->GetEntries(); ++j) {
      chain->GetEntry(j);
      if (cut_nan(dst) && 
          cut_zenith_poles(dst))// &&
          //(dst->nhit > 8)) //cut_smt_trigger(TriggID_1006) && //cut_reco_idtag(dst)*/)
      {
        // Extra cuts on the line fit
        isGoodLineFit = cut_nan(dst, kTRUE) && cut_zenith_poles(dst, kTRUE);

        // Get string hit statistics
        nstring = dst->nstring;
        nchan = dst->nchan;
        ndir = dst->ndir;

        // Get energy and direct length
        logMuE = dst->logE;
        rlogl = dst->rlogl;

        // Get COG information
        sdcog = dst->sdcog;
        cogx = dst->cog.GetX();
        cogy = dst->cog.GetY();
        cogz = dst->cog.GetZ();
        ldir = dst->ldir;

        // Get event time
        mjd = dst->mjd * units::day;

        // Change in IC79: DST now bins data in 10-us chunks, and TDST
        // converts 10-us time bins into microseconds.  So remove conversion:
        //usec = dst->time * 1e-2 * microsecond;

        usec = dst->time * units::microsecond;
        mjdTime = mjd + usec;
        secsInDay = (mjdTime - mjd) / units::second;
        nsInDay = (secsInDay - TMath::Floor(secsInDay)) / units::nanosecond;
        localMST = GetGMST(mjd);
        mjdTime /= units::day;

        // Extract geometry info
        llhAzimuth = dst->reco2.direction.GetAzimuth(); 
        llhZenith = dst->reco2.direction.GetZenith();
        linllhOpeningAngle = linllh_opening_angle(dst);

        // Convert to equatorial coordinates
        Direction dir(llhZenith,llhAzimuth);
        Equatorial eq = GetEquatorialFromDirection(dir, mjdTime);

        RA = eq.ra/units::degree;
        Dec = eq.dec/units::degree;

        // Convert using antisidereal time
        double localAntiS = GetGMAST(mjdTime);
        RAAntiS = fmod( (eq.ra - (localMST+ localAntiS)*constants::pi/12)/ units::degree,360) /units::degree;
        DecAntiS = Dec;

        // Convert using solar time
        double tod = ( mjdTime - int(mjdTime) )* 24.;
        RASolar = fmod(eq.ra - (localMST+ tod)*constants::pi/12.,2*constants::pi)/units::degree;
        DecSolar = Dec;

        // Get the position of the sun
        Equatorial sunEq = GetEquatorialFromDirection(GetSunDirection(mjdTime), mjdTime);
        RASun = fmod(sunEq.ra / units::degree,360);
        if (RASun < 0) 
            RASun += 360;
        DecSun = sunEq.dec/units::degree;

        // Get the position of the moon
        Equatorial moonEq = GetEquatorialFromDirection(GetMoonDirection(mjdTime), mjdTime);
        RAMoon = fmod(moonEq.ra / units::degree,360);
        if (RAMoon < 0) 
            RAMoon += 360; 
        DecMoon = moonEq.dec / units::degree;


        // Get run info
        runId = dst->runId;

        cutDST->Fill();

        if (!(++nEvents % step)) {
          cutDST->AutoSave("SaveSelf");
        }
      }
      if (!(j % step)) {
        cout << '.';
        cout.flush();
      }
    }
    cout << endl;
    cout << "Passed " << nEvents << " events" << endl;
  }
  cutDST->AutoSave("SaveSelf");

  return 0;
}

