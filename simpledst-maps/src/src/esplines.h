#ifndef ESPLINES
#define ESPLINES
#if __cplusplus > 199711L
#include <photospline/splinetable.h>
#include <photospline/bspline.h>

bool ICenergyCut(
	unsigned NChannels, 
	const photospline::splinetable<> &t, 
	double zenith, 
	double emin, 
	double emax);
#endif // __cplusplus
#endif // ESPLINES
