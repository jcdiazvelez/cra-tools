
#ifndef __SimpleDST_h__
#define __SimpleDST_h__

#include <Rtypes.h>

class TChain;
class TBranch;

class SimpleDST {

 public:

  SimpleDST() { }
  SimpleDST(TChain* chain, std::string config) { SetupChain(chain, config); }
  ~SimpleDST() { }

  // Joint parameters
  Double_t ModJulDay;
  // IceCube parameters
  UShort_t NChannels;
  Bool_t IsGoodLineFit;
  Bool_t isReco;
  Float_t LLHAzimuthDeg;
  Float_t LLHZenithDeg;
  Float_t RLogL;
  UInt_t NDirHits; 
  UInt_t LDir; 

  // IceTop parameters
  Double_t SPAzimuth;
  Double_t SPZenith;
  Int_t SPFitStatus;
  Double_t LapAzimuth;
  Double_t LapZenith;
  Int_t LapFitStatus;
  Double_t LSAzimuth;
  Double_t LSZenith;
  Int_t LSFitStatus;
  Int_t nStations;
  Double_t s125;
  Double_t ss125;

  // Joint parameters
  TBranch* b_ModJulDay;
  //IceCube parameters
  TBranch* b_nchan;
  TBranch* b_isGoodLineFit;
  TBranch* b_isReco;
  TBranch* b_LLHAzimuth;
  TBranch* b_LLHZenith;
  TBranch* b_RLogL;
  TBranch* b_NDirHits;
  TBranch* b_LDir;

  //IceTop parameters
  TBranch* b_SPAzimuth;
  TBranch* b_SPZenith;
  TBranch* b_SPFitStatus;
  TBranch* b_LapAzimuth;
  TBranch* b_LapZenith;
  TBranch* b_LapFitStatus;
  TBranch* b_LSAzimuth;
  TBranch* b_LSZenith;
  TBranch* b_LSFitStatus;
  TBranch* b_nStations;
  TBranch* b_s125;
  TBranch* b_ss125;

  void SetupChain(TChain* chain, std::string config);

};

#endif // __SimpleDST_h__

