/* Dwitter.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <cmath>
#include <limits>

// #include "f2c.h"
#include "Fractl.hh"
#include "Utils.hh"

/* Table of constant values */

static double c_b2 = (float) 0.0;
static double c_b11 = 2.0;
static int c__1 = 1;

int Fractl::dwiter(double *euc, double *evc, double *cvgc,
		     double *wc, double *cvgp, double *wp,
		     double *eup, double *evp,
		     int nx, int ny, double *xydeli,
		     double dzh, double rhoc, double rhop, double *c2rd, 
		     int *iter, double *dmn,
		     int *kst, int l, double zlev)
{
  double *ws = new double[nx * ny]();	// should be initalized to 0.0
  
  /* Initialized data */

  static int itmax = 5;
  static int nquad = 0;
  static int minpts = 15;

  /* System generated locals */
  int ws_dim1, ws_offset, euc_dim1, euc_offset, evc_dim1, evc_offset, 
    cvgc_dim1, cvgc_offset, wc_dim1, wc_offset, cvgp_dim1, 
    cvgp_offset, wp_dim1, wp_offset, eup_dim1, eup_offset, evp_dim1, 
    evp_offset, i__1, i__2, i__3;
  double r__1;
  double d__1;

  /* Builtin functions */
  double pow_dd(double *, double *), sqrt(double);

  /* Local variables */
  static int i__, j, ic, jc, im1, jm1, ip1, jp1;
  static double cnt;
  static int nxm1, nym1;
  static double dfdx, dmin__, dgdy, dmax__, derr, wcvg, dsqr, cxfac, cyfac, 
    witer;
  extern /* Subroutine */
    int coprx_(double *,
	       double *,
	       int *);
  static double dmprev, errmax;
  static int itrmax;
  static double rhonrm;

  /* --ITERATIVEREFINEMENT OF WC AT A LEVEL TO MINIMIZE THE */
  /*     DIFFERENCE BETWEEN ADJUSTED VALUES OF THE HORIZONTAL CONVERGENCE */
  /*     AND VERTICAL CONVERGENCE COMPUTED USING THE VALUE OF WC (CALLED WP) */
  /*     FROM THE PREVIOUS LEVEL. */

  /*     WS- SCRATCH BUFFER TO HOLD PREVIOUS ITERATE OF WC */
  /*     +++ MODIFIED BY THIS ROUTINE +++ */
  /*     EUC- U ERROR TERM FROM THE SYNTHESIS */
  /*     EVC- V   "    "    "    "     " */
  /*     CVGC- CONVERGENCE AT CURRENT LEVEL */
  /*     WC- W ESTIMATE AT THIS LEVEL   +++ MODIFIED BY THIS ROUTINE +++ */
  /*     CVGP- CONVERGENCE AT PREVIOUS LEVEL */
  /*     WP- W FROM PREVIOUS LEVEL */
  /*     +++ MODIFIED BY THIS ROUTINE +++ */
  /*     EUP- U ERROR TERM FROM PREVIOUS LEVEL */
  /*     EVP- V   "    "    "    "     " */
  /*     NX- NUMBER OF POINTS ALONG X-AXIS */
  /*     NY- NUMBER OF POINTS ALONG Y-AXIS */
  /*     XYDELI- INVERSE OF HORIZONTAL AXES SPACINGS (X AND Y IN KM) */
  /*     DZH- (1./DELZ) * SIGN (DIRECTION OF INTEGRATION) */
  /*     RHOC- DENSITY AT CURRENT LEVEL */
  /*     RHOP-    "     " PREVIOUS  " */
  /*     C2RD- USER SPECIFIED ITERATION CONTROLS: */
  /*     1- TARGET DIFFERENCR BETWEEN CONSECUTIVE ITERATIONS; */
  /*     2- MAXIMUM NUMBER OF ITERATIONS */
  /*     BAD- MISSING DATA FLAG */
  /*     ITER- ACTUAL NUMBER OF ITERATIONS              (OUTPUT) */
  /*     DMN- ACTUAL DIFFERENCE BETWEEN ITERATIONS        " */
  /*     KST- STATUS FLAG:                                " */
  /*     0- CONVERGED */
  /*     1- DIVERGED OR MAXED OUT ON ITERATIONS */
  /*     2- NO DATA THIS LEVEL */


  /*     PARAMETERIZATION FOR BOUNDED DATA FILL */

  /* Parameter adjustments */
  evp_dim1 = nx;
  evp_offset = 1 + evp_dim1;
  evp -= evp_offset;
  eup_dim1 = nx;
  eup_offset = 1 + eup_dim1;
  eup -= eup_offset;
  wp_dim1 = nx;
  wp_offset = 1 + wp_dim1;
  wp -= wp_offset;
  cvgp_dim1 = nx;
  cvgp_offset = 1 + cvgp_dim1;
  cvgp -= cvgp_offset;
  wc_dim1 = nx;
  wc_offset = 1 + wc_dim1;
  wc -= wc_offset;
  cvgc_dim1 = nx;
  cvgc_offset = 1 + cvgc_dim1;
  cvgc -= cvgc_offset;
  evc_dim1 = nx;
  evc_offset = 1 + evc_dim1;
  evc -= evc_offset;
  euc_dim1 = nx;
  euc_offset = 1 + euc_dim1;
  euc -= euc_offset;
  ws_dim1 = nx;
  ws_offset = 1 + ws_dim1;
  ws -= ws_offset;
  --xydeli;
  --c2rd;

  /* Function Body */

  *iter = 1;
  errmax = c2rd[1];
  itrmax = c2rd[2];
  nxm1 = nx - 1;
  nym1 = ny - 1;
  dmprev = (float)1e8;
  cxfac = xydeli[1] * (float)-.5;
  cyfac = xydeli[2] * (float)-.5;
  rhonrm = (float)1. / rhoc;
  i__1 = nx * ny;
  // confld_(&ws[ws_offset], &i__1, &c_b2);
  fill_array(&ws[ws_offset],  i__1, c_b2);

  /*     DETERMINE WHICH LOCATIONS ARE */
  /*     GOOD AND SET WC=0.0 THERE */
  i__1 = nym1;
  for (j = 2; j <= i__1; ++j) {
    jm1 = j - 1;
    jp1 = j + 1;
    i__2 = nxm1;
    for (i__ = 2; i__ <= i__2; ++i__) {
      wc[i__ + j * wc_dim1] = (float)0.;
      if (std::isnan(wp[i__ + j * wp_dim1])
	  || std::isnan(euc[i__ + j * euc_dim1])
	  || std::isnan(eup[i__ + j * eup_dim1])
	  || std::isnan(evc[i__ + j * evc_dim1])
	  || std::isnan(evp[i__ + j * evp_dim1])
	  || std::isnan(cvgc[i__ + j * cvgc_dim1])
	  || std::isnan(cvgp[i__ + j * cvgp_dim1])) {
	wc[i__ + j * wc_dim1] = numeric_limits<double>::quiet_NaN();
      }
      if (! std::isnan(wc[i__ + j * wc_dim1])) {
	i__3 = jp1;
	for (jc = jm1; jc <= i__3; jc += 2) {
	  if (std::isnan(wp[i__ + jc * wp_dim1])
	      || std::isnan(euc[i__ + jc * euc_dim1])
	      || std::isnan(eup[i__ + jc * eup_dim1])
	      || std::isnan(evc[i__ + jc * evc_dim1])
	      || std::isnan(evp[i__ + jc * evp_dim1])) {
	    goto L3;
	  }
	  /* L1: */
	}
	im1 = i__ - 1;
	ip1 = i__ + 1;
	i__3 = ip1;
	for (ic = im1; ic <= i__3; ic += 2) {
	  if (std::isnan(wp[ic + j * wp_dim1])
	      || std::isnan(euc[ic + j * euc_dim1])
	      || std::isnan(eup[ic + j * eup_dim1])
	      || std::isnan(evc[ic + j * evc_dim1])
	      || std::isnan(evp[ic + j * evp_dim1])) {
	    goto L3;
	  }
	  /* L2: */
	}
	goto L4;
      L3:
	wc[i__ + j * wc_dim1] = numeric_limits<double>::quiet_NaN();
      L4:
	;
      }
      /* L5: */
    }
  }

 L10:
  /*     ITERATION LOOP */
  cnt = (float)0.;
  *dmn = (float)0.;
  dsqr = (float)0.;
  dmax__ = (float)1e-35;
  dmin__ = (float)1e35;
  i__2 = nym1;
  for (j = 2; j <= i__2; ++j) {
    i__1 = nxm1;
    for (i__ = 2; i__ <= i__1; ++i__) {
      if ( std::isnan(wc[i__ + j * wc_dim1]) ) {
	goto L50;
      }
      dfdx = (eup[i__ + 1 + j * eup_dim1] * rhop * wp[i__ + 1 + j * 
						      wp_dim1] + euc[i__ + 1 + j * euc_dim1] * rhoc * ws[i__ + 
													 1 + j * ws_dim1] - (eup[i__ - 1 + j * eup_dim1] * rhop * 
															     wp[i__ - 1 + j * wp_dim1] + euc[i__ - 1 + j * euc_dim1] * 
															     rhoc * ws[i__ - 1 + j * ws_dim1])) * cxfac;
      dgdy = (evp[i__ + (j + 1) * evp_dim1] * rhop * wp[i__ + (j + 1) *
							wp_dim1] + evc[i__ + (j + 1) * evc_dim1] * rhoc * ws[
													     i__ + (j + 1) * ws_dim1] - (evp[i__ + (j - 1) * evp_dim1] 
																	 * rhop * wp[i__ + (j - 1) * wp_dim1] + evc[i__ + (j - 1) 
																						    * evc_dim1] * rhoc * ws[i__ + (j - 1) * ws_dim1])) * 
	cyfac;
      wcvg = dfdx + dgdy;
      witer = (wp[i__ + j * wp_dim1] * rhop + dzh * (wcvg + (rhoc * 
							     cvgc[i__ + j * cvgc_dim1] + rhop * cvgp[i__ + j * 
												     cvgp_dim1]))) * rhonrm;
      /*     IF(WITER.GT.100.) THEN */
      /* +++++DEBUGS */
      /*     TYPE *, 'WP(I,J),RHOP,DZH,WCVG,RHOC,CVGC(I,J),CVGP(I,J),RHONRM', */
      /*     X         WP(I,J),RHOP,DZH,WCVG,RHOC,CVGC(I,J),CVGP(I,J),RHONRM */
      /*     TYPE *, 'WC(I,J),EUP(I+1,J),WP(I+1,J),EUC(I+1,J),WS(I+1,J)', */
      /*     X         WC(I,J),EUP(I+1,J),WP(I+1,J),EUC(I+1,J),WS(I+1,J) */
      /*     TYPE *, 'CXFAC,CYFAC,WITER', CXFAC,CYFAC,WITER */
      /*     END IF */
      /* +++++ */
      cnt += (float)1.;
      derr = (r__1 = wc[i__ + j * wc_dim1] - witer, std::abs(r__1));
      *dmn += derr;
      if (std::abs(derr) < dmin__) {
	dmin__ = std::abs(derr);
      }
      if (std::abs(derr) > dmax__) {
	dmax__ = std::abs(derr);
      }
      /* Computing 2nd power */
      r__1 = derr;
      dsqr += r__1 * r__1;
      wc[i__ + j * wc_dim1] = witer;
    L50:
      ;
    }
  }

  /*     TEST CONVERGENCE */

  if (cnt == (float)0.) {
    goto L92;
  }
  *dmn /= cnt;
  dsqr /= cnt;
  d__1 = (double) (*dmn);
  //  dsqr -= pow_dd(&d__1, &c_b11);
  dsqr -= pow(d__1, c_b11);
  dsqr = sqrt(dsqr);
  if ( (*dmn <= errmax && itrmax > 0)
       || ((*dmn <= errmax) && (itrmax < 0) && (*iter >= abs(itrmax)))) {
    goto L90;
  }
  if (*iter >= abs(itrmax)) {
    goto L91;
  }
  if (*dmn - dmprev >= errmax && itrmax > 0) {
    goto L91;
  }
  /*     RESET VALUES AND DO IT AGAIN */
  ++(*iter);
  dmprev = *dmn;
  i__1 = nym1;
  for (j = 2; j <= i__1; ++j) {
    i__2 = nxm1;
    for (i__ = 2; i__ <= i__2; ++i__) {
      if ( ! std::isnan(wc[i__ + j * wc_dim1]) ) {
	ws[i__ + j * ws_dim1] = wc[i__ + j * wc_dim1];
      }
      /* L55: */
    }
  }
  goto L10;
  /*     FINISHED WITH ITERATIONS */
 L90:
  *kst = 0;
  goto L95;
 L91:
  *kst = 1;
  goto L95;
 L92:
  *kst = 2;
  return 0;
 L95:

  /*     FILL IN W-VALUES AT DATA BOUNDARIES */
  /*     WS IS DECISION FILED FOR BNDFIL */

  i__2 = nx * ny;
  // confld_(&ws[ws_offset], &i__2, &c_b2);
  fill_array(&ws[ws_offset], i__2, c_b2);
  i__2 = ny;
  for (j = 1; j <= i__2; ++j) {
    i__1 = nx;
    for (i__ = 1; i__ <= i__1; ++i__) {
      if ( ! std::isnan(wc[i__ + j * wc_dim1]) ) {
	goto L100;
      }
      if (std::isnan(wp[i__ + j * wp_dim1])
	  || std::isnan(euc[i__ + j * euc_dim1])
	  || std::isnan(eup[i__ + j * eup_dim1])
	  || std::isnan(evc[i__ + j * evc_dim1])
	  || std::isnan(evp[i__ + j * evp_dim1])
	  || std::isnan(cvgc[i__ + j * cvgc_dim1])
	  || std::isnan(cvgp[i__ + j * cvgp_dim1])) {
	goto L100;
      }

      /*     WE WANT ESTIMATE OF W HERE */

      ws[i__ + j * ws_dim1] = numeric_limits<double>::quiet_NaN();

    L100:
      /* L105: */
      ;
    }
    /* L110: */
  }

  /*     AT THIS POINT DATA FROM PREVIOUS */
  /*     LEVEL ARE NO LONGER NEEDED --- */
  /*     USE FOR DATA FILL */

  i__2 = nx * ny;
  // coprx_(&wp[wp_offset], &wc[wc_offset], &i__2);
  copy_array(&wp[wp_offset], &wc[wc_offset], i__2);  
  bndfil(&wc[wc_offset],
	 &wp[wp_offset],
	 &ws[ws_offset],
	 nx,
	 ny,
	 c__1,
	 nx,
	 c__1,
	 ny,
	 itmax,
	 nquad,
	 minpts,
	 numeric_limits<double>::quiet_NaN());

  return 0;
} /* dwiter_ */

int Fractl::bndfil(double *c__, double *a, double *b,
		   int ni, int nj,
		   int ibeg, int iend,
		   int jbeg, int jend,
		   int itmax, int nquad,
		   int minpts, double bad)
{
    /* Initialized data */

    static double eps = (float)1e-5;

    /* System generated locals */
    int a_dim1, a_offset, c_dim1, c_offset, b_dim1, b_offset, i__1, i__2, 
	    i__3, i__4, i__5;

    /* Local variables */
    static int i__, j, k, l, i1, j1, j2, i2;
    static double t1, t2, t3, am[3], bm[3], cm[3], dm[3];
    static int jo, io, kq, ix, iy;
    static double deno, anum;
    static int iquad[4];
    static double ptsmin;


/*        PERFORMS LEAST-SQUARES DATA FILLING OF A BOUNDED REGION */

    /* Parameter adjustments */
    b_dim1 = ni;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = ni;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = ni;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    ptsmin = (double) (minpts);
    i__1 = jend;
    for (jo = jbeg; jo <= i__1; ++jo) {
	i__2 = iend;
	for (io = ibeg; io <= i__2; ++io) {
	    if (a[io + jo * a_dim1] != bad || b[io + jo * b_dim1] != bad) {
	      // goto L50;
	      continue;
	    }
	    for (l = 1; l <= 3; ++l) {
		am[l - 1] = (float)0.;
		bm[l - 1] = (float)0.;
		cm[l - 1] = (float)0.;
		dm[l - 1] = (float)0.;
/* L15: */
	    }
	    for (l = 1; l <= 4; ++l) {
		iquad[l - 1] = 0;
/* L16: */
	    }
	    i__3 = itmax;
	    for (l = 1; l <= i__3; ++l) {
/* Computing MAX */
		i__4 = 1, i__5 = jo - l;
		j1 = max(i__4,i__5);
/* Computing MIN */
		i__4 = nj, i__5 = jo + l;
		j2 = min(i__4,i__5);
/* Computing MAX */
		i__4 = 1, i__5 = io - l;
		i1 = max(i__4,i__5);
/* Computing MIN */
		i__4 = ni, i__5 = io + l;
		i2 = min(i__4,i__5);
		i__4 = j2;
		for (j = j1; j <= i__4; ++j) {
		    iy = j - jo;
		    i__5 = i2;
		    for (i__ = i1; i__ <= i__5; ++i__) {
			ix = i__ - io;
			if (abs(ix) != l && abs(iy) != l) {
			  // goto L20;
			  continue;
			}
			if (a[i__ + j * a_dim1] == bad) {
			  // goto L20;
			  continue;
			}
			if (ix >= 0 && iy > 0) {
			    iquad[0] = 1;
			}
			if (ix > 0 && iy <= 0) {
			    iquad[1] = 1;
			}
			if (ix <= 0 && iy < 0) {
			    iquad[2] = 1;
			}
			if (ix < 0 && iy >= 0) {
			    iquad[3] = 1;
			}
			am[0] += (float)1.;
			am[1] += ix;
			am[2] += iy;
			bm[1] += ix * ix;
			bm[2] += ix * iy;
			cm[2] += iy * iy;
			dm[0] += a[i__ + j * a_dim1];
			dm[1] += ix * a[i__ + j * a_dim1];
			dm[2] += iy * a[i__ + j * a_dim1];
// L20:
			;
		    }
		}
		kq = 0;
		for (k = 1; k <= 4; ++k) {
/* L25: */
		    kq += iquad[k - 1];
		}
		if (kq < nquad) {
		  // goto L30;
		  continue;
		}
		if (am[0] < ptsmin) {
		  // goto L30;
		  continue;
		}
		bm[0] = am[1];
		cm[0] = am[2];
		cm[1] = bm[2];
		t1 = bm[1] * cm[2] - bm[2] * cm[1];
		t2 = bm[0] * cm[2] - bm[2] * cm[0];
		t3 = bm[0] * cm[1] - bm[1] * cm[0];
		deno = am[0] * t1 - am[1] * t2 + am[2] * t3;
		if (deno <= eps) {
		  // goto L30;
		  continue;
		}
		anum = dm[0] * t1 - dm[1] * t2 + dm[2] * t3;
		c__[io + jo * c_dim1] = anum / deno;
		goto L50;
// L30:
		;
	    }
L50:
	    ;
	}
    }
    return 0;
} /* bndfil_ */
