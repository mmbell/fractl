#include <cstddef>

#include <stdlib.h>
#include <cmath>
#include "Interps.hh"

// This is a C++ version of Jim Leise code from the Cedric source tree
// https://github.com/NCAR/RadarAnalysisTools/blob/master/cedric/trunk/EXTEND.f
// -- translated by f2c (version 20100827) with a few hand fixes.

/* Table of constant values */

// static int c__2 = 2;

// This one comes from the f2c library
// https://github.com/juanjosegarciaripoll/f2c/blob/master/lib/pow_ii.c

int pow_ii(int ap, int bp)
{
  int pow, x, n;
  unsigned long u;

  x = ap;
  n = bp;

  if (n <= 0) {
    if (n == 0 || x == 1)
      return 1;
    if (x != -1)
      return x == 0 ? 1/x : 0;
    n = -n;
  }
  u = n;
  for(pow = 1; ; )
    {
      if(u & 01)
	pow *= x;
      if(u >>= 1)
	x *= x;
      else
	break;
    }
  return(pow);
}

LeiseInterp::LeiseInterp()
{
}

LeiseInterp::~LeiseInterp()
{
}

bool LeiseInterp::interpolate(double *y, int n1, int n2, int n3, double /* flag */)
{
  /* System generated locals */
  int i__1, i__2, i__3, i__4, i__5, i__6;
  double r__1, r__2, r__3;

  /* Builtin functions */
  //  double log(doubledouble);
  //  int pow_ii(int *, int *);

  /* Local variables */
  static int i__, j, k, n;
  static double x;
  static int k1, m1, kk, kn, ln, net[5];
  static double yms;
  static int ndim, kord[5], kmax;
  static double yscl, ymax;
  static int ntot, kmark;
  static double ymean, ymark;
  static int knext;
  static double ylast;
  static int istop, jstop, kount, kstop, mstop, kstrt;
  static double ynext;
  static int mstrt;
  static double factor;

  /*                                             JIM LEISE 10/80 */
  /*    ************************************************************ */
  /*    HELLO, */
  /*    I AM A MULTIDIMENSIONAL INTERPOLATION/EXTRAPOLATION SCHEME */
  /*    WHICH NEEDS NO EXTRA ARRAY SPACE.  BAD OR MISSING POINTS Y(K) */
  /*    ARE SET BY THE USER:  Y(K)=FLAG.  THIS ROUTINE THEN REPLACES */
  /*    SUCH FLAGGED VALUES WITH AN ITERATIVE REPLACMENT METHOD. */
  /*          Y(N1,N2,N3)=MULTIDIMENSIONAL ARRAY (INPUT OR OUTPUT). */
  /*          FLAG       =BAD DATA VALUE OR INDICATOR. */
  /*    ************************************************************ */


  /*    START THE 3-D ARITHMETIC. */
  /* Parameter adjustments */
  --y;

  /* Function Body */
  ndim = 1;
  if (n2 > 1) {
    ndim = 2;
  }
  if (n3 > 1) {
    ndim = 3;
  }
  kord[0] = max(1,n1);
  kord[1] = max(1,n2);
  kord[2] = max(1,n3);
  kord[3] = kord[0];
  kord[4] = kord[1];
  ntot = kord[0] * kord[1] * kord[2];
  kmax = max(kord[0],kord[1]);
  kmax = max(kmax,kord[2]);
  if (kmax <= 1) {
    return 0;
  }
  net[0] = 1;
  net[1] = kord[0];
  net[2] = kord[0] * kord[1];
  net[3] = net[0];
  net[4] = net[1];

  /*    COMPUTE THE MEAN OF Y AND THE MAXIMUM OF ABS(Y). */
  ymax = (float)0.;
  ymean = (float)0.;
  kount = 0;
  i__1 = ntot;
  for (k = 1; k <= i__1; ++k) {
    //    if (y[k] == flag) {
    if (std::isnan(y[k])) {
      goto L10;
    }
    ymean += y[k];
    ++kount;
    /* Computing MAX */
    r__2 = ymax, r__3 = (r__1 = y[k], abs(r__1));
    ymax = max(r__2,r__3);
  L10:
    ;
  }
  if (kount == 0 || kount == ntot) {
    return 0;
  }
  ymean /= kount;
  ymax *= (float)1.5;
  if (ymax == (float)0.) {
    ymax = (float)1.;
  }

  /*    COMPUTE A SCALE FOR SHIFTING THE DATA TO POSITIVE VALUES. */
  /*    NEGATIVE NUMBERS ARE RESERVED FOR THE INTERPOLATION/ */
  /*    EXTRAPOLATION PROCEDURE.  INITIALLY, THE FLAGGED VALUES ARE */
  /*    PRESET WITH THE MEAN (SHIFTED NEGATIVE). */
  yscl = ymax * (float)3.;
  yms = -ymean - ymax;
  i__1 = ntot;
  for (k = 1; k <= i__1; ++k) {
    //    if (y[k] == flag) {
    if (std::isnan(y[k])) {
      goto L20;
    }
    y[k] += ymax;
    goto L30;
  L20:
    y[k] = yms;
  L30:
    ;
  }

  /*    ***** START THE MAIN LOOP ***** */
  x = (double) kmax;
  /*    FIRST, COMPUTE THE LARGEST POWER OF 2 .LE.KMAX. */
  k1 = log(x) / (float).693;
  //  k1 = pow_ii(c__2, k1);
  k1 = pow_ii(2, k1);
 L40:
  if (k1 <= 0) {
    goto L300;
  }
  i__1 = ndim;
  for (n = 1; n <= i__1; ++n) {

    /*    THE 3-D ARITHMETIC. */
    m1 = k1 * net[n - 1];
    istop = kord[n];
    jstop = kord[n + 1];
    i__2 = istop;
    for (i__ = 1; i__ <= i__2; ++i__) {
      i__3 = jstop;
      for (j = 1; j <= i__3; ++j) {
	kstrt = (i__ - 1) * net[n] + 1 + (j - 1) * net[n + 1];
	kstop = kstrt + (kord[n - 1] - 1) * net[n - 1];
	if (kstrt + m1 > kstop) {
	  goto L90;
	}
	kn = kstrt - net[n - 1];
	i__4 = k1;
	for (kk = 1; kk <= i__4; ++kk) {
	  kn += net[n - 1];
	  ln = kn + (kstop - kn) / m1 * m1;
	  if (ln <= kn) {
	    goto L80;
	  }

	  /*    INITIALIZE THE ITERATION. */
	  kmark = ((r__1 = y[kn], abs(r__1)) + ymax) / yscl;
	  knext = ((r__1 = y[kn + m1], abs(r__1)) + ymax) / yscl;
	  ymark = (r__1 = y[kn], abs(r__1)) - kmark * yscl;
	  ynext = (r__1 = y[kn + m1], abs(r__1)) - knext * yscl;
	  if (y[kn] > (float)0.) {
	    goto L50;
	  }
	  factor = kmark + (float)1.;
	  y[kn] = -factor * yscl - (kmark * ymark + ynext) / factor;

	  /*   FILL THE CENTER MOD(M1). */
	L50:
	  mstrt = kn + m1;
	  mstop = ln - m1;
	  if (mstrt > mstop) {
	    goto L70;
	  }
	  i__5 = mstop;
	  i__6 = m1;
	  for (k = mstrt; i__6 < 0 ? k >= i__5 : k <= i__5; k += 
		 i__6) {
	    ylast = ymark;
	    ymark = ynext;
	    kmark = knext;
	    knext = ((r__1 = y[k + m1], abs(r__1)) + ymax) / 
	      yscl;
	    ynext = (r__1 = y[k + m1], abs(r__1)) - knext * yscl;
	    if (y[k] > (float)0.) {
	      goto L60;
	    }
	    factor = kmark + (float)2.;
	    y[k] = -factor * yscl - (kmark * ymark + ylast + 
				     ynext) / factor;
	  L60:
	    ;
	  }
	L70:
	  if (y[ln] > (float)0.) {
	    goto L80;
	  }
	  factor = knext + (float)1.;
	  y[ln] = -factor * yscl - (knext * ynext + ymark) / factor;
	L80:
	  ;
	}
      L90:
	;
      }
    }
    /* L100: */
  }

  /*   NORMALIZE BEFORE CHANGING SCALES. */
  i__1 = ntot;
  for (k = 1; k <= i__1; ++k) {
    if (y[k] > (float)0.) {
      goto L200;
    }
    kmark = (-y[k] + ymax) / yscl;
    y[k] += kmark * yscl;
  L200:
    ;
  }
  /*    UPDATE THE SAMPLING INCREMENT. */
  k1 /= 2;
  goto L40;

  /*    RETURN THE DATA TO ITS ORIGINAL DYNAMIC RANGE. */
 L300:
  i__1 = ntot;
  for (k = 1; k <= i__1; ++k) {
    y[k] = (r__1 = y[k], abs(r__1)) - ymax;
    /* L310: */
  }

  return 0;
} /* extend_ */

