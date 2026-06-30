#ifndef ESPLINES
#define ESPLINES
#if __cplusplus > 199711L
#include <photospline/splinetable.h>
#include <photospline/bspline.h>
#include <SimpleDST.h>

/**
 *	Spline-based energy cuts for In-Ice SimpleDST
 *
 */
bool ICenergyCut(SimpleDST dst, photospline::splinetable<> &spline, double zenith, double emin, double emax);
#endif // __cplusplus


/**
 *	IceTop NStations cut
 *
 */
bool ITNstatCut(SimpleDST dst, double smin, double smax);

/**
 *	IceTop s125 cut
 *
 */
bool ITs125Cut(SimpleDST dst, double smin, double smax);
#endif // ESPLINES
