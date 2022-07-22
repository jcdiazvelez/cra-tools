{
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

//gROOT->LoadMacro("SpaceAngleRad.C");
// Unit definitions and forward declarations

const Double_t r2d = TMath::RadToDeg();

//Output file
TFile* outfile = new TFile("converted.root", "recreate");


//Opening files in a chain

char files[300];
  
//sprintf(files, "/data/sim/IceCube/2010/filtered/level2/CORSIKA-in-ice/6443/00000-00999/Level2_IC79_corsika.006443.000001.root");
sprintf(files, "Level2_IC79_corsika.006443.000000.root");


cout << "Opening  files in " << files <<  endl;

TChain *chain = new TChain("tree");
chain->Add(files);  


//Create Tree

TTree *MCEvents = new TTree("CutDST", "CutDST");
MCEvents->SetAutoSave(10000);

double llhazimuth;
MCEvents->Branch("LLHAzimuthDeg", &llhazimuth, "llhazimuth/D");
MCEvents->Branch("NChannels", &nChannel, "nChannel/i");
double rlogl;
MCEvents->Branch("RLogL", &rlogl, "rlog/D");
double pAzimuth;

double tFlux;
double pslope;
double rigidity;
double Nevents;
double tEnergy;

Int_t entries = chain->GetEntries();
Int_t Ntrees = chain->GetNtrees();


cout << "Number of opened files: " << Ntrees <<  endl;
cout << "Loop over " << entries << " entries" << endl;

for(int i = 0; i < entries; i++)
{
  //  cout << "i = "  << i << endl;
  chain->GetEvent(i);

  mjd           = chain->GetLeaf("time")->GetValue()/3600/24.;
  llhazimuth    = chain->GetLeaf("reco2.direction.azimuth_")->GetValue();
  llhzenith     = chain->GetLeaf("reco2.direction.zenith_")->GetValue();
  cogz          = chain->GetLeaf("cog.cog.z_")->GetValue();
  nChannel      = chain->GetLeaf("nchan")->GetValue();
  ndirhits      = chain->GetLeaf("ndir")->GetValue();
  rlogl         = chain->GetLeaf("rlogl")->GetValue();

  //Angles in degree
  llhazimuth *= r2d;
  llhzenith *= r2d;

  MCEvents->Fill();

}
cout << MCEvents->GetEntries() << endl;
outfile->cd();
outfile->Add(MCEvents);
outfile->Write();
outfile->Close();

}
