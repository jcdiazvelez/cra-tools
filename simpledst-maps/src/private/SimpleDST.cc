
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

    // Configuration-specific reconstructions
    reco = "ShowerPlane";
    if (config == "IT81-2013")
      reco = "LaputopStandard";
    else if (config == "ITpass2")
      reco = "LaputopSmall";

    std::string recoazstr(reco + ".azimuth");
    std::string recozstr(reco + ".zenith");
    std::string recoexstr(reco + ".exists");
    const char* recoaz = recoazstr.c_str();
    const char* recoz = recozstr.c_str();
    const char* recoex = recoexstr.c_str();
    chain->SetBranchAddress(recoaz, &Azimuth, &b_Azimuth);
    chain->SetBranchAddress(recoz,  &Zenith, &b_Zenith);
    chain->SetBranchAddress(recoex, &isReco, &b_isReco);

    // IceTop Filters
    if (config == "ITpass2"){
      chain->SetBranchAddress("QFilterMask.IceTopSTA3_13", &isSTA3, &b_isSTA3);
      chain->SetBranchAddress("QFilterMask.IceTopSTA5_13", &isSTA5, &b_isSTA5);
    }
    if (config == "IT81-2013"){
      chain->SetBranchAddress("FilterMask.IceTopSTA3_13", &isSTA3, &b_isSTA3);
      chain->SetBranchAddress("FilterMask.IceTopSTA5_13", &isSTA5, &b_isSTA5);
    }
    if (config == "IT81-2012"){
      chain->SetBranchAddress("QFilterMask.IceTopSTA3_12", &isSTA3, &b_isSTA3);
      chain->SetBranchAddress("QFilterMask.IceTopSTA5_12", &isSTA5, &b_isSTA5);
    }
    if (config == "IT81-2011"){
      chain->SetBranchAddress("FilterMask.IceTopSTA3_11", &isSTA3, &b_isSTA3);
      chain->SetBranchAddress("FilterMask.IceTopSTA8_11", &isSTA8, &b_isSTA8);
      chain->SetBranchAddress("FilterMask.IceTopSTA3_InIceSMT_11", 
              &isSTA3ii, &b_isSTA3ii);
    }
    if (config == "IT73"){
      chain->SetBranchAddress("FilterMask.IceTopSTA3_10", &isSTA3, &b_isSTA3);
      chain->SetBranchAddress("FilterMask.IceTopSTA8_10", &isSTA8, &b_isSTA8);
      chain->SetBranchAddress("FilterMask.IceTopSTA3_InIceSMT_10", 
              &isSTA3ii, &b_isSTA3ii);
    }
    if (config == "IT59"){
      chain->SetBranchAddress("FilterMask.IceTopSTA3_09", &isSTA3, &b_isSTA3);
      chain->SetBranchAddress("FilterMask.IceTopSTA8_09", &isSTA8, &b_isSTA8);
      chain->SetBranchAddress("FilterMask.IceTopSTA3_InIceSMT_09", 
              &isSTA3ii, &b_isSTA3ii);
    }
    

    //Laputop values
    if (config == "ITpass2"){
        chain->SetBranchAddress("LaputopParams.s125", &s125, &b_s125);
        chain->SetBranchAddress("LaputopSmallParams.s125", &ss125, &b_ss125);
    } else {
        chain->SetBranchAddress("LaputopStandardParams.s125", &s125, &b_s125);
        chain->SetBranchAddress("LaputopSmallShowerParams.s125", &ss125, &b_ss125);

        // ShowerLLH values
        chain->SetBranchAddress("ShowerLLH_proton.energy", &pEnergy, &b_pEnergy);
        chain->SetBranchAddress("ShowerLLH_iron.energy", &fEnergy, &b_fEnergy);
        chain->SetBranchAddress("maxLLH_proton.value", &pLLH, &b_pLLH);
        chain->SetBranchAddress("maxLLH_iron.value", &fLLH, &b_fLLH);

    }

  }

}

