#include <NStations.h>
#include <string>
#include <TBranch.h>
#include <TChain.h>

void
NStations::SetupChain(TChain* chain, std::string config)
{
  std::string detector = config.substr(0,2);

  if (detector == "ITpass2") {

    chain->SetBranchAddress("value", &nStations, &b_nStations);
  }

}
