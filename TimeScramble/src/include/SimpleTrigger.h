
#ifndef __SimpleTrigger_h__
#define __SimpleTrigger_h__

#include <Rtypes.h>

class TChain;
class TBranch;

class SimpleTrigger {

 public:

  SimpleTrigger() { }
  SimpleTrigger(TChain* chain) { SetupChain(chain); }
  ~SimpleTrigger() { }

  Bool_t TriggID_1006;

  TBranch* b_TriggID_1006;

  void SetupChain(TChain* chain);

};

#endif // __SimpleTrigger_h__

