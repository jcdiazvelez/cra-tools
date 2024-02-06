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

/**
 *	Spline-based energy reconstruction for In-Ice SimpleDST
 *
 */
double ICenergy(SimpleDST dst, photospline::splinetable<> &spline, double zenith);
#endif // __cplusplus


/**
 *	IceTop energy cut
 *
 */
int ITenergyCut(SimpleDST dst, double emin, double emax);

/**
 *	IceTop s125 cut
 *
 */
int ITs125Cut(SimpleDST dst, double smin, double smax);
#endif // ESPLINES
