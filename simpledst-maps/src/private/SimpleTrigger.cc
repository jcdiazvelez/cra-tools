
#include <SimpleTrigger.h>
#include <string>
#include <TBranch.h>
#include <TChain.h>

void
SimpleTrigger::SetupChain(TChain* chain)
{
  chain->SetBranchAddress("TriggID_1006", &TriggID_1006, &b_TriggID_1006);
}


