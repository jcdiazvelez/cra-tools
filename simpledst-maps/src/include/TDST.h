/**
    copyright  (C) 2004
    the icecube collaboration
    $Id: dst.h 16031 2006-02-20 12:38:45Z troy $

    @version $Revision: 1.2 $
    @date $Date: 2006-02-20 06:38:45 -0600 (lun, 20 feb 2006) $
    @author juancarlos@icecube.wisc.edu
    @brief TOject dataclass for storing DST data in a ROOT tree
*/

#ifndef TDST_H_INCLUDED
#define TDST_H_INCLUDED

#include "units.h"
#include "Direction.h"
#include "Position.h"
#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#endif
#include <map>

using namespace std;

class TDSTReco {

	public:
		TDSTReco();
		virtual ~TDSTReco();

		Direction direction;   
		Float_t distance;  
		UShort_t reco_id;   

	    inline void SetDirection(double theta, double phi) { direction = Direction(theta,phi); }
        inline void SetDistance(double d) { distance = d; }
        inline void SetID(UShort_t id) { reco_id = id; }

};



class TDST {
	public:

		TDST();
		virtual ~TDST();
					  

		// 1st reconstruction - LF                 
		TDSTReco reco1;   

		// 2nd reconstruction - DP
		TDSTReco reco2;   

		// Number of DOMs (after cleaning) 
		UShort_t nchan;   

		// Number of DOMs (after cleaning) 
		UShort_t nstring;   

		// Number of Hits (from FE or Portia)  (3+3 B)
		UInt_t nhit;   

		// Event time  
		Double_t time;   

		// Delta_t 
		Double_t dt;   

		// Event duration 
		Double_t duration;   

		// Event mjd
		Double_t mjd;   

		// Trigger Tag   (2 B)
		UInt_t triggertag;   
		// Run ID (2 B)

		UInt_t runId;   
		// Event ID (2 B)

		UInt_t eventId;   

		UInt_t subEventId;   

		// Year (1 B)
		UInt_t year;   

		// NDir number of direct hits for selected reco
		UInt_t ndir;   

		// NDir direct length for selected reco
		UInt_t ldir;   

		// NDir direct length for selected reco
		Float_t rlogl;   

		// log of Energy
		Float_t logE;   

		Position cog;   
		Float_t sdcog;   
		Float_t linefitspeed;   

		// Event Weight  (for Simulated events)
		Double_t weight;   
		Double_t diplopiaweight;   
		Double_t TimeScale;   

};

//I3_POINTER_TYPEDEFS(TDSTReco);
//I3_POINTER_TYPEDEFS(TDST);

#endif
