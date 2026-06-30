
#include <SimpleDST.h>
#include <string>
#include <TBranch.h>
#include <TChain.h>

void
SimpleDST::SetupChain(TChain* chain, std::string config)
{
  std::string detector = config.substr(0,2);

  if (detector == "IC") {

    chain->SetBranchAddress("NChannels", &NChannels, &b_nchan);
    chain->SetBranchAddress("IsGoodLineFit", &isReco, &b_isReco);
    chain->SetBranchAddress("ModJulDay", &ModJulDay, &b_ModJulDay);
    chain->SetBranchAddress("LLHAzimuthDeg", &LLHAzimuthDeg, &b_LLHAzimuth);
    chain->SetBranchAddress("LLHZenithDeg", &LLHZenithDeg, &b_LLHZenith);
    chain->SetBranchAddress("RLogL", &RLogL, &b_RLogL);
    chain->SetBranchAddress("NDirHits", &NDirHits, &b_NDirHits);
    chain->SetBranchAddress("LDir", &LDir, &b_LDir);
  }

  if (detector == "IT") {

    // Generic parameters
    std::string reco;
    chain->SetBranchAddress("NStations.value", &nStations, &b_nStations);
    chain->SetBranchAddress("I3EventHeader.time_start_mjd", 
              &ModJulDay, &b_ModJulDay);
    chain->SetBranchAddress("ShowerPlane.azimuth", &SPAzimuth, &b_SPAzimuth);
    chain->SetBranchAddress("ShowerPlane.zenith", &SPZenith, &b_SPZenith);
    chain->SetBranchAddress("ShowerPlane.fit_status", &SPFitStatus, &b_SPFitStatus);

    //Laputop values
    chain->SetBranchAddress("Laputop.azimuth", &LapAzimuth, &b_LapAzimuth);
    chain->SetBranchAddress("Laputop.zenith", &LapZenith, &b_LapZenith);
    chain->SetBranchAddress("Laputop.fit_status", &LapFitStatus, &b_LapFitStatus);
    chain->SetBranchAddress("LaputopSmall.azimuth", &LSAzimuth, &b_LSAzimuth);
    chain->SetBranchAddress("LaputopSmall.zenith", &LSZenith, &b_LSZenith);
    chain->SetBranchAddress("LaputopSmall.fit_status", &LSFitStatus, &b_LSFitStatus);
    chain->SetBranchAddress("LaputopParams.s125", &s125, &b_s125);
    chain->SetBranchAddress("LaputopSmallParams.s125", &ss125, &b_ss125);

  }

}

