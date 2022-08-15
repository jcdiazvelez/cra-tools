
#ifndef __NStations_h__
#define __NStations_h__

#include <Rtypes.h>

class TChain;
class TBranch;

class NStaions {

 public:

  NStations() { }
  NStations(TChain* chain, std::string config) { SetupChain(chain, config); }
  ~NStations() { }

  Int_t value;

  TBranch* b_value;


  void SetupChain(TChain* chain, std::string config);

};

#endif // __NStations_h__
