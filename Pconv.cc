/* PCONVG.F -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
   on Microsoft Windows system, link with libf2c.lib;
   on Linux or Unix systems, link with .../path/to/libf2c.a -lm
   or, if you install libf2c.a in a standard place, with -lf2c -lm
   -- in that order, at the end of the command line, as in
   cc *.o -lf2c -lm
   Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

   http://www.netlib.org/f2c/libf2c.zip
*/
#include <limits>

#include "Fractl.hh"

int Fractl::pconvg(double *u, double *v, double *rbuf, int nx, 
		   int ny, int nder, double *xydeli,
		   double *csp, int iactc)
{
  /* Initialized data */

  static int kpos[5] = { 0,1,2,0,3 };
  static double wx[125]	/* was [5][5][5] */ = { (float)0.,(float)0.,(float)0.,
						(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(
													     float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)
						0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(
														float)0.,(float)0.,(float)0.,(float)-1.,(float)1.,(float)0.,(
																					     float)0.,(float)0.,(float)-1.,(float)1.,(float)0.,(float)0.,(
																													  float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)
						0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(
														float)0.,(float)0.,(float)0.,(float)-1.5,(float)2.,(float)-.5,(
																					       float)0.,(float)0.,(float)-.5,(float)0.,(float).5,(float)0.,(
																													    float)0.,(float).5,(float)-2.,(float)1.5,(float)0.,(float)0.,(
																																					  float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)
						0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(
														float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)
						0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(
														float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)0.,(float)
						0.,(float)0.,(float)0.,(float)-2.08333333,(float)4.,(float)-3.,(
														float)1.333333,(float)-.25,(float)-.25,(float)-.8333333,(float)
						1.5,(float)-.5,(float).08333333,(float).08333333,(float)-.6666667,
						(float)0.,(float).6666667,(float)-.08333333,(float)-.08333333,(
													       float).5,(float)-1.5,(float).8333333,(float).25,(float).25,(float)
						-1.333333,(float)3.,(float)-4.,(float)2.08333333 };

  /* System generated locals */
  int u_dim1, u_offset, v_dim1, v_offset, rbuf_dim1, rbuf_offset, i__1, 
    i__2, i__3, i__4;

  /* Local variables */
  static int i__, j, k, m, ii, jj, kl, ip, kr;
  static double xx;
  static int khi, igo, idx, klo;
  static double sum;
  static int idel;
  static double corr;
  static int mpos, npts;


  /*        CALCULATES THE CONVERGENCE USING U,V COMPONENTS */
  /*        NDER CONSECUTIVE VALUES MUST BE PRESENT TO COMPUTE */
  /*             THE PARTIAL ALONG A LINEAR AXIS. PARTIAL WILL */
  /*             BE PRESENT AT NDER LOCATIONS SINCE EDGE WEIGHTING */
  /*             IS PERFORMED. */

  /*        MODIFIED TO PERFORM EDGE WEIGHTING-- C. MOHR 8/87 */

  /* Parameter adjustments */
  rbuf_dim1 = nx;
  rbuf_offset = 1 + rbuf_dim1;
  rbuf -= rbuf_offset;
  v_dim1 = nx;
  v_offset = 1 + v_dim1;
  v -= v_offset;
  u_dim1 = nx;
  u_offset = 1 + u_dim1;
  u -= u_offset;
  --xydeli;
  csp -= 4;

  /* Function Body */

  /*        PARTIAL ACROSS THE X-DIMENSION */

  i__1 = ny;
  for (jj = 1; jj <= i__1; ++jj) {
    j = jj;
    klo = 0;
    i__2 = nx;
    for (ii = 1; ii <= i__2; ++ii) {
      i__ = ii;
      if ( ! std::isnan(u[i__ + j * u_dim1]) ) {
	if (klo == 0) {
	  klo = i__;
	}
	if (i__ != nx) {
	  goto L20;
	}
	i__ = nx + 1;
      }
      if (klo == 0) {
	goto L16;
      }
      khi = i__ - 1;
      k = i__ - klo;
      if (k <= 2) {
	/*                                        MUST BE AT LEAST 3 CONSEC. PTS. */
	rbuf[klo + j * rbuf_dim1] = numeric_limits<double>::quiet_NaN();
	rbuf[khi + j * rbuf_dim1] = numeric_limits<double>::quiet_NaN();
      } else {
	/*                                        GOT AT LEAST 3, 5 FALLS BACK TO */
	/*                                        3 IF NDER=5. */
	npts = min(k,nder);
	if (npts == 4) {
	  npts = 3;
	}
	idel = npts - kpos[npts - 1];
	i__3 = khi;
	for (m = klo; m <= i__3; ++m) {
	  kl = m - klo;
	  kr = khi - m;
	  if (kl < idel) {
	    /*                                        LEFT  OF CENTER WEIGHTING */
	    mpos = kl + 1;
	  } else if (kr < idel) {
	    /*                                        RIGHT OF CENTER WEIGHTING */
	    mpos = npts - kr;
	  } else {
	    /*                                              CENTER WEIGHTING */
	    mpos = kpos[npts - 1];
	  }
	  igo = m - mpos;
	  sum = (float)0.;
	  i__4 = npts;
	  for (ip = 1; ip <= i__4; ++ip) {
	    idx = ip + igo;
	    sum += u[idx + j * u_dim1] * wx[ip + (mpos + npts * 5)
					    * 5 - 31];
	    if (ip == mpos && iactc == 1) {

	      /*     CORRECTION FOR COPLANE SPACE */

	      xx = csp[4] + (idx - 1) * csp[6];
	      if (xx != (float)0.) {
		corr = u[idx + j * u_dim1] / xx;
	      } else {
		rbuf[m + j * rbuf_dim1] = numeric_limits<double>::quiet_NaN();
		goto L23;
	      }
	    }
	    /* L9: */
	  }
	  sum *= xydeli[1];
	  if (iactc == 1) {
	    rbuf[m + j * rbuf_dim1] = sum + corr;
	  } else {
	    rbuf[m + j * rbuf_dim1] = sum;
	  }
	L23:
	  /* L10: */
	  ;
	}
      }
      if (i__ > nx) {
	goto L20;
      }
      klo = 0;

    L16:
      rbuf[i__ + j * rbuf_dim1] = numeric_limits<double>::quiet_NaN();

    L20:
      ;
    }
  }

  /*        PARTIAL ACROSS THE Y-DIMENSION */

  i__2 = nx;
  for (ii = 1; ii <= i__2; ++ii) {
    i__ = ii;
    klo = 0;
    i__1 = ny;
    for (jj = 1; jj <= i__1; ++jj) {
      j = jj;
      if ( ! std::isnan(v[i__ + j * v_dim1]) ) {
	if (klo == 0) {
	  klo = j;
	}
	if (j != ny) {
	  goto L40;
	}
	j = ny + 1;
      }
      if (klo == 0) {
	goto L36;
      }
      khi = j - 1;
      k = j - klo;
      if (k <= 2) {
	/*                                        MUST BE AT LEAST 3 CONSEC. PTS. */
	rbuf[i__ + klo * rbuf_dim1] = numeric_limits<double>::quiet_NaN();
	rbuf[i__ + khi * rbuf_dim1] = numeric_limits<double>::quiet_NaN();
      } else {
	/*                                        GOT AT LEAST 3, 5 FALLS BACK TO */
	/*                                        3 IF NDER=5. */
	npts = min(k,nder);
	if (npts == 4) {
	  npts = 3;
	}
	idel = npts - kpos[npts - 1];
	i__3 = khi;
	for (m = klo; m <= i__3; ++m) {
	  kl = m - klo;
	  kr = khi - m;
	  if (kl < idel) {
	    /*                                        LEFT  OF CENTER WEIGHTING */
	    mpos = kl + 1;
	  } else if (kr < idel) {
	    /*                                        RIGHT OF CENTER WEIGHTING */
	    mpos = npts - kr;
	  } else {
	    /*                                              CENTER WEIGHTING */
	    mpos = kpos[npts - 1];
	  }
	  igo = m - mpos;
	  sum = (float)0.;
	  i__4 = npts;
	  for (ip = 1; ip <= i__4; ++ip) {
	    idx = ip + igo;
	    sum += v[i__ + idx * v_dim1] * wx[ip + (mpos + npts * 
						    5) * 5 - 31];
	    /* L29: */
	  }
	  if ( ! std::isnan(rbuf[i__ + m * rbuf_dim1]) ) {
	    rbuf[i__ + m * rbuf_dim1] = -(rbuf[i__ + m * 
					       rbuf_dim1] + sum * xydeli[2]);
	  }
	  /* L30: */
	}
      }
      if (j > ny) {
	goto L40;
      }
      klo = 0;

    L36:
      rbuf[i__ + j * rbuf_dim1] = numeric_limits<double>::quiet_NaN();

    L40:
      ;
    }
  }

  return 0;
} /* pconvg */
