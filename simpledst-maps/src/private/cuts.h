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

/**
 *  GB IceTop stations cut --> low and high
 *
 */
int ITstats(Simble DST dst, double loen, double hien);
#endif //
