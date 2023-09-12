/**
    copyright  (C) 2004
    the icecube collaboration
    $Id: dst.cxx 16031 2006-02-20 12:38:45Z troy $

    @version $Revision: 1.2 $
    @date $Date: 2006-02-20 06:38:45 -0600 (lun, 20 feb 2006) $
    @author juancarlos@icecube.wisc.edu
*/

#include "TDST.h"

#include <iostream>
#include <cmath>
#include "TBuffer.h"
#include "TMemberInspector.h"

using std::cout;

#ifdef I3_USE_CINT
ClassImp(TDST)
#endif

TDST::~TDST() { }


TDST::TDST():
	reco1(),
	reco2(),
	nchan(0),
	nhit(0),
	time(0),
	dt(0),
	duration(0),
	triggertag(0),
	runId(0),
	eventId(0),
	subEventId(0),
	rlogl(0.),
	logE(0.),
	cog(0,0,0),
	sdcog(0),
	linefitspeed(0.),
	weight(1.0),
	diplopiaweight(1.0),
	TimeScale(1.0)
{
}


TDSTReco::~TDSTReco() { }


TDSTReco::TDSTReco():
	direction(NAN,NAN),
	distance(NAN),
	reco_id(0)
{
}

