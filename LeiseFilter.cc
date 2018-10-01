#include <cstddef>
#include "Filters.hh"

LeiseFilter::LeiseFilter()
{
}

LeiseFilter::~LeiseFilter()
{
}

// This is a C++ version of Jim Leise code from the Cedric source tree
// https://github.com/NCAR/RadarAnalysisTools/blob/master/cedric/trunk/T5FLTR.f
// -- translated by f2c (version 20100827) with a few hand fixes.

// The ns parameter is what was passed in the COMMON section in the
// original Fortran code. If not passed in it will default to all 0's

bool LeiseFilter::filter(double *y, int n1, int n2, int n3,
			 int nstep, int *ns)
{
  /* System generated locals */
  int i__1, i__2, i__3, i__4, i__5, i__6, i__7;

  /* Local variables */
  static int i__, j, k, m, n, k1, m1, m2, m3, kn, ln;
  static double ym1, ym2;
  static int net[5], nns[3];
  static double ykn, yln, ykn1, yln1;
  static int main, ndim, kord[5];
  static double ysave;
  static int mstep, istop, jstop, kstop, mstop, kstrt, mstrt, mpyrmd;

  /*                                               JIM LEISE 8/80 */
  /*    ***************************************************************** */
  /*    HELLO, */
  /*    I AM A MULTIDIMENSIONAL LOW-PASS FILTER WHICH NEEDS NO EXTRA */
  /*    ARRAY SPACE.  THUS, THE FILTERED ANSWER IS RETURNED IN THE SAME */
  /*    ARRAY Y(N1,N2,N3) THAT THE DATA IS INPUT.  THE CENTRAL FILTER */
  /*    IS A LINEAR 5-PT FILTER AND THE BOUNDARY FILTER IS COMPUTED */
  /*    USING A MIRROR EXTENSION.  THUS, THE TOTAL FILTER IS LINEAR. */

  /*          ********** NSTEP CONTROL FOR 1-DIM ********** */
  /*        STEP RESTRICTION:  5*2**(NSTEP-1) .LE. MAX(N1,N2,N3) */
  /*         PASSBAND .LE. 2**(NSTEP+2)  POINTS/CYCLE */
  /*         STOPBAND .GE. 2**(NSTEP)  POINTS/CYCLE. */

  /*          ********** MULTIDIMENSIONAL USE ********** */
  /*    PARAMETER CONTROL FOR THE THREE DIMENSIONS CAN BE DOUBLEIZED */
  /*    VIA COMMON/FLTRPL/ WHERE NS CORRESPONDS TO NSTEP.  IF THIS */
  /*    COMMON IS NOT USED, THE VALUES OF NS ARE DEFAULTED TO NSTEP */
  /*    -I.E. NSTEP IS USED IN PLACE OF ANY ZEROS. */
  /*    ****************************************************************** */

  /*    INITIALIZATION OF NS FOR CSD APPLICATIONS (4/15/82) */
  /* Parameter adjustments */
  --y;

  /* Function Body */

  /*    INITIALIZE THE 3-D ARITHMETIC. */
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
  net[0] = 1;
  net[1] = kord[0];
  net[2] = kord[0] * kord[1];
  net[3] = net[0];
  net[4] = net[1];

  /*    DEFAULT PARAMETER TRANSFER. */
  mpyrmd = 0;
  i__1 = ndim;

  // Default ns values if ns was not passed in
  int default_ns[ndim] = {0};
  if (ns == NULL) {
    ns = default_ns;
  }

  for (n = 1; n <= i__1; ++n) {
    nns[n - 1] = ns[n - 1];
    if (ns[n - 1] == 0) {
      nns[n - 1] = nstep;
    }
    if (kord[n - 1] < 5) {
      nns[n - 1] = 0;
    }
    /* L10: */
    /* Computing MAX */
    i__2 = mpyrmd, i__3 = nns[n - 1] + nns[n - 1] - 1;
    mpyrmd = max(i__2,i__3);
  }

  // -------------------

  if (mpyrmd <= 0) {
    return false;
  }
  mstep = (mpyrmd + 1) / 2;

  /*    ***** START THE MAIN LOOP ***** */
  k1 = 1;
  i__2 = mpyrmd;
  for (main = 1; main <= i__2; ++main) {
    i__3 = ndim;
    for (n = 1; n <= i__3; ++n) {
      /*    SAMPLING CHECKS. */
      if (k1 * 10 > kord[n - 1]) {
	/* Computing MIN */
	i__1 = nns[n - 1];
	nns[n - 1] = min(i__1,main);
      }
      if (main >= nns[n - 1] && mpyrmd - main >= nns[n - 1]) {
	goto L40;
      }

      /*    THE 3-D ARITHMETIC. */
      m1 = k1 * net[n - 1];
      m2 = m1 + m1;
      m3 = m2 + m1;
      istop = kord[n];
      jstop = kord[n + 1];
      i__1 = istop;
      for (i__ = 1; i__ <= i__1; ++i__) {
	i__4 = jstop;
	for (j = 1; j <= i__4; ++j) {
	  kstrt = (i__ - 1) * net[n] + 1 + (j - 1) * net[n + 1];
	  kstop = kstrt + (kord[n - 1] - 1) * net[n - 1];
	  kn = kstrt - net[n - 1];
	  i__5 = k1;
	  for (k = 1; k <= i__5; ++k) {
	    kn += net[n - 1];
	    ln = kn + (kstop - kn) / m1 * m1;

	    /*    FILTER THE ENDS USING A MIRROR EXTENSION. */
	    ykn = y[kn] * (double).875 + y[kn + m1] * (double).1875
	      - y[kn + m2] * (double).0625;
	    yln = y[ln] * (double).875 + y[ln - m1] * (double).1875
	      - y[ln - m2] * (double).0625;
	    ykn1 = y[kn] * (double).1875 + y[kn + m1] * (double)
	      .625 + y[kn + m2] * (double).25 - y[kn + m3] *
	      (double).0625;
	    yln1 = y[ln] * (double).1875 + y[ln - m1] * (double)
	      .625 + y[ln - m2] * (double).25 - y[ln - m3] *
	      (double).0625;

	    /*    DO THE CENTRAL 5-PT FILTER. */
	    ym2 = y[kn];
	    ym1 = y[kn + m1];
	    mstrt = kn + m2;
	    mstop = ln - m2;

	    i__6 = mstop;
	    i__7 = m1;
	    for (m = mstrt; i__7 < 0 ? m >= i__6 : m <= i__6; m +=
		   i__7) {
	      ysave = y[m];
	      y[m] = y[m] * (double).625
		+ (ym1 + y[m + m1])
		* (double).25 - (ym2 + y[m + m2]) * (double) .0625;
	      ym2 = ym1;
	      /* L20: */
	      ym1 = ysave;
	    }

	    y[kn + m1] = ykn1;
	    y[ln - m1] = yln1;
	    y[kn] = ykn;
	    y[ln] = yln;

	    /* L30: */
	  }
	}
      }
    L40:
      ;
    }
    /*    UPDATE THE SAMPLING INCREMENT. */
    k1 += k1;
    if (main >= mstep) {
      k1 /= 4;
    }
    /* L50: */
  }

  return true;
} /* t5fltr_ */
