// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR, CSU  (c) 1990 - 2018
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Colorado State University
// ** BSD licence applies
// ** DISCLAIMER: THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS
// ** OR IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
// ** WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

#include <iostream>
#include <string>
#include "Fractl.hh"
#include "Params.hh"
#include "Filters.hh"
#include "Interps.hh"

// Eigen linear algebra library
// http://eigen.tuxfamily.org/index.php?title=Main_Page
#include <Eigen/Dense>

double sumTimea = 0;
int cntTimeb = 0;
double sumTimeb = 0;
double sumTimec = 0;
double sumTimee = 0;
double sumTimef = 0;

bool Fractl::calcWinds()
{
  // One approach is to iterate:
  //   calculate V and U winds
  //   calculate W winds
  //   recalculate V and U winds using the new W winds
  //   recalculate W with the new V, U.
  //   etc.
  //
  // However, tests show that there is negligible gain
  // in accuracy after the first iteration.
  // So the loop limit is hardcoded at 1.

  for (int imain = 0; imain < 1; imain++) {
    bool showStat = true;
    // Based on current w, calc u, v.
    if (gridType == Params::GRID_MISH)
      calcAllVUOnMish(
		numNbrMax,       // max num nearest nbrs
		pointVec,        // all observations
		radarKdTree,     // nearest nbr tree for pointVec
		cellMat);
    else
      calcAllVU(
		numNbrMax,       // max num nearest nbrs
		pointVec,        // all observations
		radarKdTree,     // nearest nbr tree for pointVec
		cellMat);

      
    printRunTime("calcAllVU", &timea);

    // Run low-pass filtering if requested

    Filter *filter = FilterFactory::createFilter(uvFilter);		// TOD indices
    if (filter != NULL) {
      filter->filter_U(cellMat, nradx, nrady, nradz, 1, NULL);
      filter->filter_V(cellMat, nradx, nrady, nradz, 1, NULL);
      delete filter;
    }

    // Run interpolation for missing data if requested
    // TODO: add old RadarWind interpolation and code Interp::interpolate

    Interp *interp = InterpFactory::createInterp(uvInterp);
    if (interp != NULL) {
      interp->interpolate_U(cellMat, nradx, nrady, nradz, std::numeric_limits<double>::quiet_NaN());
      interp->interpolate_V(cellMat, nradx, nrady, nradz, std::numeric_limits<double>::quiet_NaN());
      delete interp;
    }

    // NO LONGER USED:
    // Interpolate missing values from neighbors.
    //
    // The problem is that some cells don't have enough valid
    // neighbors to calculate their winds, so their values are
    // left as NaN.
    //
    // The idea was to fill the NaN cells
    // by interpolating the wind values from nearby cells.
    //
    // But what happens is that generally if a cell doesn't have
    // enough good neighbors for a valid wind calculation,
    // any neighbor cells we might use for interpolation are
    // pretty flakey.
    //xxx del interpolation?
    // Interpolate NaN values in cellMat.uu, cellMat.vv.
    //interpMissing( bugs, nradz, nrady, nradx, cellMat, tmpmat);
    //printRunTime("interpMissing for V,U", &timea);

    // Based on u, v, calc w

    if (gridType == Params::GRID_MISH)    
      calcAllWOnMish(cellMat);
    else
      calcAllW(cellMat);
    
    printRunTime("calcAllW", &timea);

    // Run low-pass filtering on W if requested

    filter = FilterFactory::createFilter(wFilter);
    if (filter != NULL) {
      filter->filter_W(cellMat, nradx, nrady, nradz, 1, NULL);  // TODO indices
      delete filter;
    }

    // Calc deltas using verification data, if any.
    // Write outTxt file.
    checkVerif(imain, cellMat);
    printRunTime("checkVerif", &timea);
  } // for imain
  return true;
}

//======================================================================

// Calculate V and U winds at all cells in the z,y,x grid.

void Fractl::calcAllVU(
		       long numNbrMax,              // max num nearest nbrs
		       vector<Point *> *pointVec,   // all observations
		       KD_tree * radarKdTree,       // nearest nbr tree for pointVec
		       Cell ***& cellMat)           // we set Cell.uu, vv
{

  KD_real * centerLoc = new KD_real[ndim];
  long okCount = 0;

  for (long iz = 0; iz < nradz; iz++) {
    for (long iy = 0; iy < nrady; iy++) {
      for (long ix = 0; ix < nradx; ix++) {

        centerLoc[0] = zgridmin + iz * zgridinc;
        centerLoc[1] = ygridmin + iy * ygridinc;
        centerLoc[2] = xgridmin + ix * xgridinc;

        if (bugs >= Params::DEBUG_NORM) {
          cout << setprecision(5);
          cout << endl << "calcAllVU: iz: " << iz
	       << "  iy: " << iy
	       << "  ix: " << ix
	       << "  z: " << centerLoc[0]
	       << "  y: " << centerLoc[1]
	       << "  x: " << centerLoc[2]
	       << endl;
        }

        Cell * pcell = & cellMat[iz][iy][ix];
        calcCellVU(
		   centerLoc,             // query point
		   numNbrMax,             // max num nearest nbrs
		   pointVec,              // all observations
		   radarKdTree,           // nearest nbr tree for pointVec
		   pcell);                 // Cell.vv, uu are set.

        bool ok = false;
        if (isOkDouble( pcell->vv)
	    && isOkDouble( pcell->uu))
	  {
	    ok = true;
	    okCount++;
	  }

        if (bugs >= Params::DEBUG_VERBOSE) {
          cout << setprecision(7);
          cout << "calcAllVU: ok:" << ok
	       << "  iz: " << iz
	       << "  iy: " << iy
	       << "  ix: " << ix
	       << "  loc:"
	       << "  " << centerLoc[0]
	       << "  " << centerLoc[1]
	       << "  " << centerLoc[2]
	       << "  W: " << cellMat[iz][iy][ix].ww
	       << "  V: " << cellMat[iz][iy][ix].vv
	       << "  U: " << cellMat[iz][iy][ix].uu << endl;
        }

      } // for ix
    } // for iy
  } // for iz

  cout << "---------------- okCount: " << okCount << endl;

  cout << "sumTimea: " << sumTimea << endl;
  cout << "cntTimeb: " << cntTimeb << "  sumTimeb: " << sumTimeb << endl;
  cout << "sumTimec: " << sumTimec << endl;
  cout << "sumTimee: " << sumTimee << endl;
  cout << "sumTimef: " << sumTimef << endl;

  delete[] centerLoc;

} // end calcAllVU

//======================================================================

// Try to calc V and U on the Samurai Mish grid

void Fractl::calcAllVUOnMish(
			     long numNbrMax,              // max num nearest nbrs
			     vector<Point *> *pointVec,   // all observations
			     KD_tree * radarKdTree,       // nearest nbr tree for pointVec
			     Cell ***& cellMat)           // we set Cell.uu, vv
{

  KD_real * centerLoc = new KD_real[ndim];
  long okCount = 0;
  
  // for (long iz = 0; iz < nradz; iz++)
  //   for (long iy = 0; iy < nrady; iy++)
  //     for (long ix = 0; ix < nradx; ix++)

  //       centerLoc[0] = zgridmin + iz * zgridinc;
  //       centerLoc[1] = ygridmin + iy * ygridinc;
  //       centerLoc[2] = xgridmin + ix * xgridinc;

  // This is the way Samurai iterates on the Mish.
  // *Pos is the position in meter of the observation
  // i*   is the index of the position in node table
  
  for (int ii = -1; ii < (nradx); ii++) {
    for (int imu = -1; imu <= 1; imu += 2) {
      double iPos = xgridmin + xgridinc * (ii + (0.5 * sqrt(1. / 3.) * imu + 0.5));
	    
      for (int ji = -1; ji < (nrady); ji++) {
	for (int jmu = -1; jmu <= 1; jmu += 2) {
	  double jPos = ygridmin + ygridinc * (ji + (0.5 * sqrt(1. / 3.) * jmu + 0.5));
		
	  for (int ki = -1; ki < (nradz); ki++) {
	    for (int kmu = -1; kmu <= 1; kmu += 2) {
	      double kPos = zgridmin + zgridinc * (ki + (0.5 * sqrt(1. / 3.) * kmu + 0.5));
	
	      int ix = (ii + 1) * 2 + (imu + 1) / 2;
	      int iy = (ji + 1) * 2 + (jmu + 1) / 2;
	      int iz = (ki + 1) * 2 + (kmu + 1) / 2;

	      // centerLoc[0] = zgridmin + iz * zgridinc;
	      // centerLoc[1] = ygridmin + iy * ygridinc;
	      // centerLoc[2] = xgridmin + ix * xgridinc;

	      centerLoc[0] = kPos;
	      centerLoc[1] = jPos;
	      centerLoc[2] = iPos;
  
	      if (bugs >= Params::DEBUG_NORM) {
		std::cout << setprecision(5);
		std::cout << endl << "calcAllVU: iz: " << iz
			  << "  iy: " << iy
			  << "  ix: " << ix
			  << "  z: " << centerLoc[0]
			  << "  y: " << centerLoc[1]
			  << "  x: " << centerLoc[2]
			  << endl;
	      }             

	      Cell * pcell = & cellMat[iz][iy][ix];
	      calcCellVU(centerLoc,             // query point
			 numNbrMax,             // max num nearest nbrs
			 pointVec,              // all observations
			 radarKdTree,           // nearest nbr tree for pointVec
			 pcell);                 // Cell.vv, uu are set.

	      bool ok = false;
	      if (isOkDouble( pcell->vv) && isOkDouble( pcell->uu)) {
		ok = true;
		okCount++;
	      }

	      if (bugs >= Params::DEBUG_VERBOSE) {
		cout << setprecision(7);
		cout << "calcAllVU: ok:" << ok
		     << "  iz: " << iz
		     << "  iy: " << iy
		     << "  ix: " << ix
		     << "  loc:"
		     << "  " << centerLoc[0]
		     << "  " << centerLoc[1]
		     << "  " << centerLoc[2]
		     << "  W: " << cellMat[iz][iy][ix].ww
		     << "  V: " << cellMat[iz][iy][ix].vv
		     << "  U: " << cellMat[iz][iy][ix].uu << endl;
	      }
	    } 
	  } // ki
	}
      } // ji
    }
  }  // ii
	
  cout << "---------------- okCount: " << okCount << endl;

  cout << "sumTimea: " << sumTimea << endl;
  cout << "cntTimeb: " << cntTimeb << "  sumTimeb: " << sumTimeb << endl;
  cout << "sumTimec: " << sumTimec << endl;
  cout << "sumTimee: " << sumTimee << endl;
  cout << "sumTimef: " << sumTimef << endl;

  delete[] centerLoc;

} // end calcAllVUOnMish


//======================================================================

// Calculate V and U winds for a single location.
//
// At a given location centerLoc, find the nearest nbrs
// in pointVec, and use the radial velocities to calculate
// the winds V, U, (not W).
//
// Let vx, vy, vz be estimates of the wind velocity,
// and vr be the radial velocity.
//
// For each radar m, we use the linear model:
//    vr_est[m] = cos(theta) cos(elev) vx
//      + sin(theta) cos(elev) vy
//      + sin(elev) vz
//    vr_est[m] = conx vx + cony * vy + conz * vz
//
// Let E[m] = error = vr_true[m] - vr_est[m]
//
// We want to minimize Q = sum_m E[m]^2
//
//
// Let
//   sumxx = sum(conx * conx)
//   sumxr = sum(conx * vr)
//   etc.
//
// Normal equations, derived from
// 0 = d(Q)/d(vx), 0 = d(Q)/d(vy), 0 = d(Q)/d(vz),
// give:
// vx sumxx + vy sumxy + vz sumxz = sumxr
// vx sumxy + vy sumyy + vz sumyz = sumyr
// vx sumxz + vy sumyz + vz sumzz = sumzr
//
// Or for 2 dim,
// vx sumxx + vy sumxy = sumxr - vz sumxz
// vx sumxy + vy sumyy = sumyr - vz sumyz
//
// Using Cramer's rule on the 2-dim case:
// Let demom = sumxx sumyy - sumxy^2
// vx_hat = ( (sumxr - vz sumxz) sumyy - (sumyr - vz sumyz) sumxy ) / denom
// vy_hat = ( (sumyr - vz sumyz) sumxx - (sumxr - vz sumxz) sumxy ) / denom

void Fractl::calcCellVU(
			KD_real * centerLoc,           // query point: z, y, x
			long numNbrMax,                // max num nearest nbrs
			vector<Point *> *pointVec,     // all observations
			KD_tree * radarKdTree,         // nearest nbr tree for pointVec
			Cell * pcell)                  // we fill vv, uu.
{
  struct timeval timea;
  addDeltaTime( &timea, NULL);

  pcell->uu = numeric_limits<double>::quiet_NaN();
  pcell->vv = numeric_limits<double>::quiet_NaN();

  bool showDetail = testDetail(
			       centerLoc[0],         // z
			       centerLoc[1],         // y
			       centerLoc[2]);        // x

  if (bugs >= Params::DEBUG_EXTRA && showDetail) {
    cout << setprecision(7);
    cout << "calcCellVU.entry: showDetail:" << endl;
    cout << "    centerLoc: z: " << centerLoc[0] << endl;
    cout << "    centerLoc: y: " << centerLoc[1] << endl;
    cout << "    centerLoc: x: " << centerLoc[2] << endl;
    cout << "    ww: " << pcell->ww << endl;
    cout << "    vv: " << pcell->vv << endl;
    cout << "    uu: " << pcell->uu << endl;
    cout << "    meanNbrDbz: " << pcell->meanNbrDbz << endl;
    cout << "    meanNbrNcp: " << pcell->meanNbrNcp << endl;
    cout << "    meanNbrElevDeg: " << pcell->meanNbrElevDeg << endl;
    cout << "    meanNbrKeepDist: " << pcell->meanNbrKeepDist << endl;
    cout << "    meanNbrOmitDist: " << pcell->meanNbrOmitDist << endl;
    cout << "    condNum: " << pcell->conditionNumber << endl;
  }

  if (numNbrMax == 0) throwerr("cell numNbrMax == 0");

  int nbrIxs[numNbrMax];
  KD_real nbrDistSqs[numNbrMax];              // nbr dist^2
  for (long inbr = 0; inbr < numNbrMax; inbr++) {
    nbrIxs[inbr] = -1;
  }

  Statistic nbrDbzStat;
  Statistic nbrNcpStat;
  Statistic nbrElevDegStat;
  Statistic nbrKeepDistStat;
  Statistic nbrOmitDistStat;

  Point** nearPts = new Point *[numNbrMax];

  addDeltaTime( &timea, &sumTimea);
  struct timeval timeb;
  addDeltaTime( &timeb, NULL);

  // Find nearest nbrs
  radarKdTree->nnquery(
		       centerLoc,        // query point
		       numNbrMax,        // desired num nearest nbrs
		       KD_EUCLIDEAN,     // Metric
		       1,                // MinkP
		       nbrIxs,           // out: parallel array, indices of nearest nbrs
		       nbrDistSqs);      // out: parallel array, squares of distances of nbrs

  cntTimeb++;
  addDeltaTime( &timeb, &sumTimeb);
  struct timeval timec;
  addDeltaTime( &timec, NULL);

  if (bugs >= Params::DEBUG_VERBOSE && showDetail) {
    cout << setprecision(7);
    cout << "  calcCellVU: showDetail:  nearPts for centerLoc: z: "
	 << centerLoc[0] << "  y: " << centerLoc[1]
	 << "  x: " << centerLoc[2] << endl;
  }

  int numNbrActual = 0;

  // Keep points within base + factor * aircraft dist

  for (long inbr = 0; inbr < numNbrMax; inbr++) {
    if (nbrIxs[inbr] < 0) throwerr("nbrIxs < 0");

    Point * nearPt = pointVec->at( nbrIxs[inbr]);
    double aircraftDist = calcDistPtAircraft( nearPt);
    double localDist = calcDistLocPt( centerLoc, nearPt);
    double maxDist = maxDistBase + maxDistFactor * aircraftDist;

    // Use local distance constraint from grid point center
    double roi = sqrt(xgridinc * xgridinc + ygridinc * ygridinc + zgridinc * zgridinc);
    if (roi < maxDist) maxDist = roi;

    const char * msg;

    if (localDist < maxDist) {
      nbrDbzStat.addOb( nearPt->dbz);
      nbrNcpStat.addOb( nearPt->ncp);
      nbrElevDegStat.addOb( nearPt->elevRad * 180 / M_PI);
      nbrKeepDistStat.addOb( localDist);
      nearPts[numNbrActual++] = nearPt;
      msg = "KEEP";
    } // if localDist < maxDist

    else {
      nbrOmitDistStat.addOb( localDist);
      msg = "OMIT";
    }

    if (bugs >= Params::DEBUG_VERBOSE && showDetail) {
      cout << setprecision(5);
      cout << "    " << msg
	   << "  inbr: " << inbr
	   << "  dist: " << localDist
	   << "  pt:"
	   << "  coordz: " << nearPt->coordz
	   << "  coordy: " << nearPt->coordy
	   << "  coordx: " << nearPt->coordx
	   << "  vg: " << nearPt->vg
	   << endl;
    }

  } // for inbr

  if (bugs >= Params::DEBUG_EXTRA && showDetail) {
    cout << "  calcCellVU: showDetail: cell z: "
	 << centerLoc[0] << "  y: " << centerLoc[1]
	 << "  x: " << centerLoc[2] << "  numNbrActual: " << numNbrActual << endl;
  }

  if (numNbrActual >= 2) {    // if numNbrActual is ok

    pcell->meanNbrDbz = nbrDbzStat.dsum / nbrDbzStat.numGood;
    pcell->meanNbrNcp = nbrNcpStat.dsum / nbrNcpStat.numGood;
    pcell->meanNbrElevDeg = nbrElevDegStat.dsum / nbrElevDegStat.numGood;
    pcell->meanNbrKeepDist = nbrKeepDistStat.dsum / nbrKeepDistStat.numGood;
    pcell->meanNbrOmitDist = nbrOmitDistStat.dsum / nbrOmitDistStat.numGood;

    if (bugs >= Params::DEBUG_EXTRA && showDetail) {
      for (long inbr = 0; inbr < numNbrActual; inbr++) {
	Point * pt = nearPts[inbr];
	cout << setprecision(7);
	cout << "  calcCellVU: showDetail: nbr: inbr: " << inbr
	     << "  vg: " << pt->vg
	     << "  dbz: " << pt->dbz
	     << "  ncp: " << pt->ncp
	     << "  theta deg: " << (pt->thetaRad * 180 / M_PI)
	     << "  elev deg: " << (pt->elevRad * 180 / M_PI)
	     << setprecision(15)
	     << "  deltaTime: "
	     << (pt->rayTime - nearPts[0]->rayTime)
	     << endl;
      }
      for (long inbr = 0; inbr < numNbrActual; inbr++) {
	Point * pt = nearPts[inbr];
	cout << "  showDetail.from.aircraft.to.nbr: set arrow "
	     << (inbr + 1)
	     << " from " << pt->aircraftx << "," << pt->aircrafty
	     << " to " << pt->coordx << "," << pt->coordy
	     << endl;
      }
    } // if showDetail

    addDeltaTime( &timec, &sumTimec);
    struct timeval timee;
    addDeltaTime( &timee, NULL);


    // 3d solve
    //    vr = cos(theta) cos(elev) vx
    //      + sin(theta) cos(elev) vy
    //      + sin(elev) vz

    // 2d solve
    //    vr - sin(elev) vz = cos(theta) cos(elev) vx
    //      + sin(theta) cos(elev) vy


    // Find wwind = W wind estimate near the pt.
    double wwind = pcell->ww;
    bool isCellOk = false;

    if (useEigen) {

      Eigen::MatrixXd amat( numNbrActual, 2);
      Eigen::VectorXd bvec( numNbrActual);

      for (long inbr = 0; inbr < numNbrActual; inbr++) {
	Point * pt = nearPts[inbr];
	amat( inbr, 0) = cos( pt->thetaRad) * cos( pt->elevRad);
	amat( inbr, 1) = sin( pt->thetaRad) * cos( pt->elevRad);

	// Find wwind = W wind estimate near the pt.
	double wwind = pcell->ww;

	bvec( inbr) = pt->vg - wwind * sin( pt->elevRad);
      }
      if (bugs >= Params::DEBUG_EXTRA && showDetail) {
	cout << "\n  calcCellVU: eigen: amat:\n" << amat << endl;
	cout << "\n  calcCellVU: eigen: bvec:\n" << bvec << endl;
      }

      Eigen::JacobiSVD<Eigen::MatrixXd> svd(
					    amat, Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::VectorXd singVals = svd.singularValues();
      Eigen::MatrixXd thinV = svd.matrixV();
      long slen = singVals.size();
      long vlen = thinV.rows();

      pcell->conditionNumber = numeric_limits<double>::infinity();
      if (slen > 0 && singVals[slen-1] != 0)
	pcell->conditionNumber = singVals[0] / singVals[slen-1];

      if (bugs >= Params::DEBUG_EXTRA && showDetail) {
	cout << "\n  calcCellVU: eigen: singVals:\n" << singVals << endl;
	cout << "\n  calcCellVU: num nonzeroSingularValues:\n"
	     << svd.nonzeroSingularValues() << endl;
	cout << "  calcCellVU: conditionNumber: "
	     << pcell->conditionNumber << endl;
      }

      // Using a conditionNumberCutoff > 10 causes some cells
      // to have extreme values for U and V.
      double conditionNumberCutoff = 100;     // xxxx
      if (pcell->conditionNumber < conditionNumberCutoff)
	isCellOk = true;
      if (isCellOk) {
	Eigen::VectorXd xvec = svd.solve( bvec);
	Eigen::VectorXd errvec = amat * xvec - bvec;
	double maxAbsErr = errvec.array().abs().maxCoeff();
	double meanSqErr = errvec.dot( errvec) / numNbrActual;   // dot product

	pcell->vv = xvec(1);       // v
	pcell->uu = xvec(0);       // u

	double ustd = 0;
	double vstd = 0;
	
	for (long s = 0; s < slen; s++) {
	  double thin0 = thinV(0, s);
	  double thin1 = thinV(1, s);
	  double sing = singVals[s];
	  ustd += pow( (thinV(s, 0) / singVals[s]) , 2.0);
	  vstd += pow( (thinV(s, 1) / singVals[s]) , 2.0);
	}
	pcell->ustd = sqrt(ustd);
	pcell->vstd = sqrt(vstd);
      } // if isCellOk
    } // if useEigen


    // Else use cramer's method
    else {
      double sxx = 0;
      double sxy = 0;
      double syy = 0;
      double sxz = 0;
      double syz = 0;

      for (long inbr = 0; inbr < numNbrActual; inbr++) {
	Point * pt = nearPts[inbr];
	double xval = cos( pt->thetaRad) * cos( pt->elevRad);
	double yval = sin( pt->thetaRad) * cos( pt->elevRad);
	double zval = pt->vg - wwind * sin( pt->elevRad);

	sxx += xval * xval;
	sxy += xval * yval;
	syy += yval * yval;
	sxz += xval * zval;
	syz += yval * zval;
      }
      double detbase = sxx*syy - sxy*sxy;
      double det1    = sxz*syy - syz*sxy;
      double det2    = sxx*syz - sxy*sxz;
      if (bugs >= Params::DEBUG_EXTRA && showDetail) {
	cout << "  calcCellVU: cramer: showDetail: "
	     << "  sxx: " << sxx
	     << "  sxy: " << sxy
	     << "  syy: " << syy
	     << "  sxz: " << sxz
	     << "  syz: " << syz << endl;
	cout << "  calcCellVU: cramer: showDetail: "
	     << "  detbase: " << detbase
	     << "  det1: " << det1
	     << "  det2: " << det2 << endl;
      }

      // Using a detCutoff <= 0.01 causes some cells
      // to have extreme values for U and V.
      // A detCutoff of 0.1 gives roughly similar results to using
      // a conditionNumber cutoff of 10.
      double detCutoff = 0.1;    // xxxxxxxxx

      if (fabs(detbase) > detCutoff) {
	pcell->uu = det1 / detbase;
	pcell->vv = det2 / detbase;
      }
    }

    addDeltaTime( &timee, &sumTimee);
    delete[] nearPts;
  } // else numNbrActual is ok

  struct timeval timef;
  addDeltaTime( &timef, NULL);

  if (bugs >= Params::DEBUG_EXTRA && showDetail) {
    cout << setprecision(7);
    cout << "calcCellVU.exit: showDetail:" << endl;
    cout << "    centerLoc: z: " << centerLoc[0] << endl;
    cout << "    centerLoc: y: " << centerLoc[1] << endl;
    cout << "    centerLoc: x: " << centerLoc[2] << endl;
    cout << "    ww: " << pcell->ww << endl;
    cout << "    vv: " << pcell->vv << endl;
    cout << "    uu: " << pcell->uu << endl;
    cout << "    meanNbrDbz: " << pcell->meanNbrDbz << endl;
    cout << "    meanNbrNcp: " << pcell->meanNbrNcp << endl;
    cout << "    meanNbrElevDeg: " << pcell->meanNbrElevDeg << endl;
    cout << "    meanNbrKeepDist: " << pcell->meanNbrKeepDist << endl;
    cout << "    meanNbrOmitDist: " << pcell->meanNbrOmitDist << endl;
    cout << "    condNum: " << pcell->conditionNumber << endl;
  }

} // end calcCellVU

// xxx all alloc: delete

//======================================================================


// Calculate the W (vertical) winds for all cells in the z,y,x grid.
//
// Consider the vertical column of cells for some
// given horizontal coordinates iy, ix.
// By the mass continuity equation, the total volume of air coming
// into the column must equal the total volume of air leaving the column.
//
// 0 = totalFlow = sum_iz (density_iz * sideFlow_iz) + endFlow
// where
//   density_iz = the density at level iz
//   sideFlow_iz = flow through the sides of the cell at level iz
//   endFlow = flow in the bottom and out the top end of the column
//           = densityBot * bottomFlow - densityTop * topFlow
//
// We assume the end flows are 0, so
//   bottendFlow = 0
//   topFlow = 0
//   endFlow = 0
//
// sideFlow_iz = uFlow_iz + vFlow_iz
// uFlow_iz = (flow in from cell ix-1) - (flow out to cell ix+1)
// uFlow_iz = 0.5 * (u[ix-1] + u[ix]) - 0.5 * (u[ix] + u[ix+1])
//          = 0.5 * (u[ix-1] - u[ix+1])
//
// sideFlow_iz = 0.5 * (u[ix-1] - u[ix+1] + v[ix-1] - v[ix+1])
//
// 0 = totalFlow = sum_iz (density_iz * 0.5 * (
//       v[iy-1] - v[iy+1]
//     + u[ix-1] - u[ix+1]))
//
// But in real life the totalFlow sum is not 0.
//
// Now we ask what modifications to the U and V values
// would make the totalFlow 0.
// We want to modify those U,V values with larger elevation angles more.
// xxx future: handle geometric uncertainty values too.
//
// At each layer, we have 4 terms to adjust: 2 U terms and 2 V terms.
// Let us subtract h[iz] from each term on level iz.
//
// Find h[iz] such that:
//   0 = sum_iz (density_iz * 0.5 * (
//         (v[iy-1]-h[iz]) - (v[iy+1]+h[iz])
//       + (u[ix-1]-h[iz]) - (u[ix+1])+h[iz]))
//
//     = totalFlow - sum_iz (density_iz * 0.5 * 4 * h[iz])
//
//   totalFlow = 2 * sum_iz (density_iz * h[iz])
//
// Let h[iz] be weighted by the elevation angle, so
//   h[iz] = H * wgt[iz]
//   wgt[iz] = cos(elevationAngle)
//
//   totalFlow = 2 * sum_iz (density[iz] * H * wgt[iz])
//   H = totalFlow / (2 * sum_iz (density[iz] * wgt[iz]))

void Fractl::calcAllWOnMish(
		      Cell ***& cellMat)           // We set Cell.ww
{
  long maxx = (xgridmax - xgridmin) * 2 / xgridinc + 1;
  long maxy = (ygridmax - ygridmin) * 2 / ygridinc + 1;  
  long maxz = (zgridmax - zgridmin) * 2 / zgridinc + 1;
  
  double * density = new double[maxz + 3];	 // range is 0..maxz
  double * wgts = new double[maxz + 3];
  
  // double * density = new double[nradz];
  // double * wgts = new double[nradz];
  
  //  for (long iz = 0; iz < nradz; iz++)
  //    density[iz] = calcDensity( zgridmin + iz * zgridinc);

  // Use samurai index computation since we need both position and index
  
  for (int ki = -1; ki < (nradz); ki++) {
    for (int kmu = -1; kmu <= 1; kmu += 2) {
      double kPos = zgridmin + zgridinc * (ki + (0.5 * sqrt(1. / 3.) * kmu + 0.5));
      int ik = (ki + 1) * 2 + (kmu + 1) / 2;
      density[ik] = calcDensity(kPos);
    }
  }

  // Omit the edges of the region as we use ix-1, ix+1, iy-1, iy+1.
  
  // for (long iy = 1; iy < nrady - 1; iy++)
  //   for (long ix = 1; ix < nradx - 1; ix++)

  for (long iy = 1; iy < maxy - 1; iy++) {
    for (long ix = 1; ix < maxx - 1; ix++) {

      // for (long iz = 0; iz < nradz; iz++)
      for (long iz = 0; iz < maxz; iz++) {
        // xxx Future:
        // Find c = the nearest cell to this one
        // through which the aircraft flew.
        // Let cosElev = horizDistToC / slantDistToC
        // wgt = cosElev

        wgts[iz] = 1.0;
      }

      // Calc totalFlow = sum of everything flowing into the column,
      // not counting the top or bottom faces.

      double totalFlow = 0;
      double sumWgt = 0;

      // for (long iz = 0; iz < nradz; iz++)
      for (long iz = 0; iz < maxz; iz++) {
	
        if ( isOkDouble( cellMat[iz][iy][ix-1].uu)
	     && isOkDouble( cellMat[iz][iy][ix+1].uu)
	     && isOkDouble( cellMat[iz][iy-1][ix].vv)
	     && isOkDouble( cellMat[iz][iy+1][ix].vv))
	  {
	    totalFlow += density[iz] * 0.5
	      * ( cellMat[iz][iy][ix-1].uu - cellMat[iz][iy][ix+1].uu
		  +  cellMat[iz][iy-1][ix].vv - cellMat[iz][iy+1][ix].vv);
	    sumWgt += density[iz] * wgts[iz];
	  }
      } // for iz

      double hcon;
      if (fabs(sumWgt) < epsilon) hcon = 0;
      else hcon = totalFlow / (2 * sumWgt);

      // Calc W wind = totalFlow, starting at the bottom,
      // using the modified U, V winds.
      double wwind = baseW;
      // for (long iz = 0; iz < nradz; iz++)
      for (long iz = 0; iz < maxz; iz++) {
        if ( isOkDouble( cellMat[iz][iy][ix-1].uu)
	     && isOkDouble( cellMat[iz][iy][ix+1].uu)
	     && isOkDouble( cellMat[iz][iy-1][ix].vv)
	     && isOkDouble( cellMat[iz][iy+1][ix].vv))
	  {
	    wwind += density[iz] * 0.5
	      * (  cellMat[iz][iy][ix-1].uu - cellMat[iz][iy][ix+1].uu
		   + cellMat[iz][iy-1][ix].vv - cellMat[iz][iy+1][ix].vv
		   - 4 * hcon * wgts[iz]);
	    if (! isOkDouble( wwind))
	      throwerr("calcAllW: invalid w wind");
	    cellMat[iz][iy][ix].ww = wwind;
	  } else cellMat[iz][iy][ix].ww = numeric_limits<double>::quiet_NaN();
      } // for iz
      if (fabs(wwind - baseW) > epsilon) {
        cout << setprecision(15);
        cout << "calcAllW: iy: " << iy << "  ix: " << ix
	     << "  baseW: " << baseW << "  wwind: " << wwind << endl;
        cout.flush();
        throwerr("wwind error");
      }
    } // for ix
  } // for iy

  //xxx maybe omit:
  // We only calculated the inner points for iy, ix.
  // But we also will need the edge values when
  // we calc the next iteration of U, V values.
  // So set the edges from the nearest interior points.
  // May be NaN.

  // for (long iz = 0; iz < nradz; iz++)
  for (long iz = 0; iz < maxz; iz++) {
    // Corners
    cellMat[iz][0][0].ww             = cellMat[iz][1][1].ww;
    cellMat[iz][0][maxx-1].ww       = cellMat[iz][1][maxx-2].ww;
    cellMat[iz][maxy-1][0].ww       = cellMat[iz][maxy-2][1].ww;
    cellMat[iz][maxy-1][maxx-1].ww = cellMat[iz][maxy-2][maxx-2].ww;

    // Side edges
    // for (long iy = 1; iy < nrady - 1; iy++)
    for (long iy = 1; iy < maxy - 1; iy++) {    
      cellMat[iz][iy][0].ww       = cellMat[iz][iy][1].ww;
      cellMat[iz][iy][maxx-1].ww = cellMat[iz][iy][maxx-2].ww;
    }
    // Top and bottom edges
    // for (long ix = 1; ix < nradx - 1; ix++)
    for (long ix = 1; ix < maxx - 1; ix++) {
      cellMat[iz][0][ix].ww       = cellMat[iz][1][ix].ww;
      cellMat[iz][maxy-1][ix].ww = cellMat[iz][maxy-2][ix].ww;
    }
  } // for iz

  delete[] density;
  delete[] wgts;

} // end calcAllWOnMish

//==================================================================

void Fractl::calcAllW(
		      Cell ***& cellMat)           // We set Cell.ww
{
  double * density = new double[nradz];
  double * wgts = new double[nradz];
  for (long iz = 0; iz < nradz; iz++)
    density[iz] = calcDensity( zgridmin + iz * zgridinc);

  // Omit the edges of the region as we use ix-1, ix+1, iy-1, iy+1.
  for (long iy = 1; iy < nrady - 1; iy++) {
    for (long ix = 1; ix < nradx - 1; ix++) {

      for (long iz = 0; iz < nradz; iz++) {
        // xxx Future:
        // Find c = the nearest cell to this one
        // through which the aircraft flew.
        // Let cosElev = horizDistToC / slantDistToC
        // wgt = cosElev

        wgts[iz] = 1.0;
      }

      // Calc totalFlow = sum of everything flowing into the column,
      // not counting the top or bottom faces.

      double totalFlow = 0;
      double sumWgt = 0;

      for (long iz = 0; iz < nradz; iz++) {
        if ( isOkDouble( cellMat[iz][iy][ix-1].uu)
	     && isOkDouble( cellMat[iz][iy][ix+1].uu)
	     && isOkDouble( cellMat[iz][iy-1][ix].vv)
	     && isOkDouble( cellMat[iz][iy+1][ix].vv))
	  {
	    totalFlow += density[iz] * 0.5
	      * ( cellMat[iz][iy][ix-1].uu - cellMat[iz][iy][ix+1].uu
		  +  cellMat[iz][iy-1][ix].vv - cellMat[iz][iy+1][ix].vv);
	    sumWgt += density[iz] * wgts[iz];

	    bool showDetail = testDetail(
					 zgridmin + iz * zgridinc,      // z
					 ygridmin + iy * ygridinc,      // y
					 xgridmin + ix * xgridinc);      // x
	    // detailSpec);                   // z, y, x, delta

	    if (showDetail) {
	      cout << setprecision(7);
	      cout << "calcAllW: showDetail:" << endl
		   << "    iz: " << iz << endl
		   << "    iy: " << iy << endl
		   << "    ix: " << ix << endl
		   << "    den: " << density[iz] << endl
		   << "    wgt: " << wgts[iz] << endl
		   << "    U-: " << cellMat[iz][iy][ix-1].uu << endl
		   << "    U+: " << cellMat[iz][iy][ix+1].uu << endl
		   << "    V-: " << cellMat[iz][iy-1][ix].vv << endl
		   << "    V+: " << cellMat[iz][iy+1][ix].vv << endl;
	    }
	  }
      } // for iz

      double hcon;
      if (fabs(sumWgt) < epsilon) hcon = 0;
      else hcon = totalFlow / (2 * sumWgt);

      // Calc W wind = totalFlow, starting at the bottom,
      // using the modified U, V winds.
      double wwind = baseW;
      for (long iz = 0; iz < nradz; iz++) {
        if ( isOkDouble( cellMat[iz][iy][ix-1].uu)
	     && isOkDouble( cellMat[iz][iy][ix+1].uu)
	     && isOkDouble( cellMat[iz][iy-1][ix].vv)
	     && isOkDouble( cellMat[iz][iy+1][ix].vv))
	  {
	    wwind += density[iz] * 0.5
	      * (  cellMat[iz][iy][ix-1].uu - cellMat[iz][iy][ix+1].uu
		   + cellMat[iz][iy-1][ix].vv - cellMat[iz][iy+1][ix].vv
		   - 4 * hcon * wgts[iz]);
	    if (! isOkDouble( wwind))
	      throwerr("calcAllW: invalid w wind");
	    cellMat[iz][iy][ix].ww = wwind;

	    bool showDetail = testDetail(
					 zgridmin + iz * zgridinc,      // z
					 ygridmin + iy * ygridinc,      // y
					 xgridmin + ix * xgridinc);      // x

	    if (showDetail) {
	      cout << setprecision(7);
	      cout << "calcAllW: showDetail:" << endl
		   << "  iz: " << iz << endl
		   << "  iy: " << iy << endl
		   << "  ix: " << ix << endl
		   << "  den: " << density[iz] << endl
		   << "  U-: " << cellMat[iz][iy][ix-1].uu << endl
		   << "  U+: " << cellMat[iz][iy][ix+1].uu << endl
		   << "  V-: " << cellMat[iz][iy-1][ix].vv << endl
		   << "  V+: " << cellMat[iz][iy+1][ix].vv << endl
		   << "  hcon: " << hcon << endl
		   << "  wgt: " << wgts[iz] << endl
		   << "  con: " << (4 * hcon * wgts[iz]) << endl
		   << "  wwind: " << wwind << endl;
	    }
	  }
        else cellMat[iz][iy][ix].ww = numeric_limits<double>::quiet_NaN();
      } // for iz
      if (fabs(wwind - baseW) > epsilon) {
        cout << setprecision(15);
        cout << "calcAllW: iy: " << iy << "  ix: " << ix
	     << "  baseW: " << baseW << "  wwind: " << wwind << endl;
        cout.flush();
        throwerr("wwind error");
      }
    } // for ix
  } // for iy


  //xxx maybe omit:
  // We only calculated the inner points for iy, ix.
  // But we also will need the edge values when
  // we calc the next iteration of U, V values.
  // So set the edges from the nearest interior points.
  // May be NaN.

  for (long iz = 0; iz < nradz; iz++) {
    // Corners
    cellMat[iz][0][0].ww             = cellMat[iz][1][1].ww;
    cellMat[iz][0][nradx-1].ww       = cellMat[iz][1][nradx-2].ww;
    cellMat[iz][nrady-1][0].ww       = cellMat[iz][nrady-2][1].ww;
    cellMat[iz][nrady-1][nradx-1].ww = cellMat[iz][nrady-2][nradx-2].ww;

    // Side edges
    for (long iy = 1; iy < nrady - 1; iy++) {
      cellMat[iz][iy][0].ww       = cellMat[iz][iy][1].ww;
      cellMat[iz][iy][nradx-1].ww = cellMat[iz][iy][nradx-2].ww;
    }
    // Top and bottom edges
    for (long ix = 1; ix < nradx - 1; ix++) {
      cellMat[iz][0][ix].ww       = cellMat[iz][1][ix].ww;
      cellMat[iz][nrady-1][ix].ww = cellMat[iz][nrady-2][ix].ww;
    }
  } // for iz

  delete[] density;
  delete[] wgts;

} // end calcAllWOnMish

//==================================================================

// Calc distance from between a Point and the aircraft

double Fractl::calcDistPtAircraft( Point * pta) {
  double sumsq = 0;
  double delta;
  delta = pta->coordz - pta->aircraftz;
  sumsq += delta * delta;
  delta = pta->coordy - pta->aircrafty;
  sumsq += delta * delta;
  delta = pta->coordx - pta->aircraftx;
  sumsq += delta * delta;
  return sqrt( sumsq);
}

//==================================================================

// Calc distance from a location to Point, in km.

double Fractl::calcDistLocPt( double * loc, Point * pta) {
  double sumsq = 0;
  double delta;
  delta = loc[0] - pta->coordz;
  sumsq += delta * delta;
  delta = loc[1] - pta->coordy;
  sumsq += delta * delta;
  delta = loc[2] - pta->coordx;
  sumsq += delta * delta;
  return sqrt( sumsq);
}

//==================================================================

// Calc distance between two Points, in km.

double Fractl::calcDistPtPt( Point * pta, Point * ptb) {
  double sumsq = 0;
  double delta;
  delta = pta->coordz - ptb->coordz;
  sumsq += delta * delta;
  delta = pta->coordy - ptb->coordy;
  sumsq += delta * delta;
  delta = pta->coordx - ptb->coordx;
  sumsq += delta * delta;
  return sqrt( sumsq);
}
