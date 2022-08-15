#ifndef __LaputopSmall_h__
#define __LaputopSmall_h__

#include <Rtypes.h>

class TChain;
class TBranch;

class LaputopSmall {

 public:

  LaputopSmall() { }
  LaputopSmall(TChain* chain, std::string config) { SetupChain(chain, config); }
  ~LaputopSmall() { }

  Double_t azimuth;
  Double_t zenith;
  Double_t time;

  TBranch* b_azimuth;
  TBranch* b_zenith;
  TBranch* b_time;

  void SetupChain(TChain* chain, std::string config);

};

#endif // __LaputopSmall_h__

