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

// Given multiple doppler observations, calculate the U, V, W winds.
//
// For sample invocations see:
//   mkwind.sh     Invoke using Eldora aircraft data
//   mkchill.sh    Invoke using CSU-Chill data
//   testa.gdb     Debug using Eldora Aircraft data
//   testb.gdb     Debug using CSU-Chill data
//
//==================================================================
//
// Future to do:
//
// Examine: how many neighbors to use?
//
// Do 3D or 4D linear interpolation to calculate scalar means
// for dbz, ncp, etc.
//
// Maybe weight neighbors by their distance from cell centers.
// Output geometric uncertainty factors: func of angle, vg uncertainty
//
// Make and study artificial quantization errors on vg: 2^(-8) or 2^(-16)
//
// Have the radarKdTree use a different metric for horizontal distance
// than for vertical distance.
//   cvs/libs/kd/src/include/kd/kd.hh
//   cvs/libs/kd/src/kd/kd.cc
//   Add parm: double * dimWeights
//     rnnEuclidean:
//       d = dimWeights[j] * (querpoint[j] - _points[_perm[i]][j]);
//
// ==========
//
// Speed this up.  Example run times are:
//
// ./mkwind.sh  10        0.5         0.5         0.5         zetaBeltrami   0,0,0     0,0       0                3,1,-7,0.1  /d1/steves/tda/tdwind
//
//   runTime:                start: 1e-06
//   runTime:                 init: 0.000253
//   runTime:              readDir: 0.003311
//   runTime:       readRadarFiles: 43.08
//   runTime:         build KD mat: 0.27992
//   runTime:        build KD tree: 3.0377
//   runTime:        alloc cellMat: 0.13854
//   runTime:            calcAllVU: 542.88
//   runTime:             calcAllW: 0.70067
//   runTime:           checkVerif: 34.11723
//
// Speed this up by:
// 1. Replace the eigen linear algebra library with a closed
// form solution.  Could use either:
//   Miller's implementation of Cramer's rule
//   Bullock's rotation of axes in 2D
//
// 2. Use a parallel approach.  Use pthreads.
// Use a thread pool of, say, the number of available processors - 1,
// assuming numProc > 1.
// In calcAllVU have a separate thread for each z layer.
//
//
//==================================================================
//
//
// Some rough file sizes:
//
//   Bell's synthetic:   240 files, 371,000 valid obs per file
//   garden city vortex: 269 files, 437,000 valid obs per file
//   brodzik:            365 files, 585,000 valid obs per file
//
// If we have 1000 files, at 500,000 obs per file, that's 5e8 obs.
// Each obs has:
//   lat, lon, alt, vg, dbz, ncp
//   about 6 * 4 = 24 bytes
//   Plus some memory allocation overhead.
//
// So the total memory in this case is about 1.2e10 bytes.
//
//
//==================================================================
//
//
// Eldora fields:
//   name: VT   longName: Radial Velocity, Combined
//   name: ZZ   longName: Radar Reflectivity Factor, Combined
//   name: VV   longName: Radial Velocity, Combined
//   name: NCP  longName: Normalized Coherent Power, Combined
//   name: SW   longName: Spectral Width, Combined
//   name: DBZ  longName: Radar Reflectivity Factor, Combined
//   name: VR   longName: Radial Velocity, Combined
//   name: VG   longName: Ground relative velocity
//   name: GG   longName: Ground Gates
//   name: SWZ  longName: Ratio
//
//
// CSU-Chill fields:
//   name: DZ  longName: reflectivity
//   name: VE  longName: radial velocity
//   name: DR  longName: differential reflectivity
//   name: DP  longName: differential phase
//   name: RH  longName: correlation H-to-V
//   name: LH  longName: Linear depolarization ratio, v-rx
//   name: LV  longName: Linear depolarization ratio, v-tx, h-rx
//   name: NC  longName: normalized coherent power
//   name: CH  longName: mag correlation HH to VH
//   name: CV  longName: mag correlation VV to HV
//   name: XH
//   name: XV
//
//==================================================================
//





#include <algorithm>
#include <cerrno>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <map>
#include <ctype.h>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <cstring>
#include <vector>
#include <random>

#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>

#include "Bbox.hh"
#include "Cell.hh"
#include "FileSpec.hh"
#include "Point.hh"
#include "Statistic.hh"
#include "Fractl.hh"
#include "Params.hh"

#include "FractlConfig.h" // version number

using namespace std;


//==================================================================

// TODO: Pull option parsing from Fractl class and do it in main.

// Driver

int main( int argc, char * argv[])
{
  int retVal = 0;
  // Fractl * rwind = new Fractl( argc, argv);
  Fractl rwind;
  if ( ! rwind.run(argc, argv) )
    retVal = 1;

  return retVal;
}


//==================================================================

const GeographicLib::TransverseMercatorExact Fractl::tranMerc =
  GeographicLib::TransverseMercatorExact::UTM();


//==================================================================

// Constructor: does everything.    --- Why?
//
// Acquire command line parms.
// Call readRadarFile to read swp radar files.
// Or if TESTMODE_ALPHA, generate synthetic data.
//
// Generate a grid of Cells.
// Call calcAllVU to calculate V and U winds.
// Call calcAllW to calculate W winds.
// Call checkVerif to calc verification deltas and statistics (if testing),
// and to write a text format of the cellMat to the outTxt file.
// Call writeNetcdf to write the variables in the cellMat
// to the Netcdf file outNc.

// TODO: These belong in a counter/timer object

int cntTimeh = 0;
double sumTimeh = 0;
int cntTimei = 0;
double sumTimei = 0;
int cntTimej = 0;
double sumTimej = 0;
int cntTimek = 0;
double sumTimek = 0;
int cntTimel = 0;
double sumTimel = 0;
int cntTimem = 0;
double sumTimem = 0;
int cntTimen = 0;
double sumTimen = 0;

Fractl::Fractl()
{
  verbose = false;
  preGridded = false;

  bugs = Params::DEBUG_OFF;
  testMode = Params::MODE_NONE;
  synWinds = NULL;
  radFiles = NULL;

  // Init timing

  // struct timeval timea;
  if (gettimeofday( &timea, NULL) != 0) throwerr("gettimeofday err");
  printRunTime("start", &timea);
  minDbz = -20;           // Default. Overwrite with
  minNcp = 0.3;           // -minDbz and -minNcp

  zgridmin = MISS_PARM;
  zgridmax = MISS_PARM;
  zgridinc = MISS_PARM;

  ygridmin = MISS_PARM;
  ygridmax = MISS_PARM;
  ygridinc = MISS_PARM;

  xgridmin = MISS_PARM;
  xgridmax = MISS_PARM;
  xgridinc = MISS_PARM;

  projName = "";
  projLat0 = MISS_PARM;
  projLon0 = MISS_PARM;
  radarAlt = -1;

  baseW = MISS_PARM;
  epsilon = MISS_PARM;
  maxDeltaAltKm = MISS_PARM;
  maxAbsElevDeg = MISS_PARM;
  minRadialDistKm = MISS_PARM;
  numNbrMax = MISS_PARM;
  maxDistBase = MISS_PARM;
  maxDistFactor = MISS_PARM;
  forceOk = false;
  useEigen = false;
  //xxx also spec max dist of observation from cell center
  inDir = "";
  fileRegex = "";
  fileList = "";
  radialName = "";
  dbzName = "";
  ncpName = "";
  outTxt = "";
  outNc = "";
  detailSpec = NULL;

  cellMat = NULL;

  ndim = 3;

  long randSeed = 1;
  srandom( randSeed);

  // Set up geodesic model and transverseMercator model
  earthRadiusMeter = GeographicLib::Constants::WGS84_a<double>();

  // flattening of WGS84 ellipsoid (1/298.257223563).
  flattening = GeographicLib::Constants::WGS84_f<double>();

  geodesic = new GeographicLib::Geodesic(earthRadiusMeter, flattening);

  pointVec = NULL;
  radarKdTree = NULL;

  aircraftBbox = new Bbox();   // overall bounding box of aircraft locs
  pointBbox = new Bbox();      // overall bounding box of point locs

  timeMin = numeric_limits<double>::quiet_NaN();
  timeMax = numeric_limits<double>::quiet_NaN();

} // end constructor


bool Fractl::run(int argc, char *argv[])
{
  parseArgs(argc, argv);

  // Given projLon0 = 148.0, projLat0 = 16.5,
  // calculate basex = 0, basey = 1.82424e+06 / 1000 = 1842

  latLonToYX(
    // tranMerc,          // TransverseMercatorExact
    projLon0,          // central meridian of projection
    0,                 // basey: coord y base, km
    0,                 // basex: coord x base, km
    projLat0,          // input lat
    projLon0,          // input lon
    basey,             // output value
    basex);            // output value

  if (verbose) {
    cout << setprecision(15);
    cout << "readRadarFile: basex: " << basex << endl;
    cout << "readRadarFile: basey: " << basey << endl;
  }

  if ( ! loadObservations() )
    return false;
  buildKdTree();
  findLimits();
  allocateCellMat();
  calcWinds();

  writeNetcdf();
  freeCellMat();

  struct rusage ruse;
  getrusage( RUSAGE_SELF, &ruse);
  cout << setprecision(5);
  cout << "final   ru_maxrss: " << ruse.ru_maxrss << endl;
  cout << "final   ru_ixrss: " << ruse.ru_ixrss << endl;
  cout << "final   ru_idrss: " << ruse.ru_idrss << endl;
  cout << "final   ru_isrss: " << ruse.ru_isrss << endl;
  return true;
}

//==================================================================

Fractl::~Fractl() {}

//==================================================================

//======================================================================



//======================================================================

// Returns the air density at a given height,
// using the US Standard Atmosphere.
//
// Hypsometric equation for US Standard Atmosphere (1976), from
// Seymour L. Hess, Introduction to Theoretical Meteorology, 1979, p 82-83.
//
// Caution: the current wikipedia article assumes constant temperature;
// Seymour Hess does not.
//
// http://scipp.ucsc.edu/outreach/balloon/atmos/1976%20Standard%20Atmosphere.htm


double Fractl::calcDensity(
  double height)                // height in km MSL
{
  double t1 = 288.15;           // kelvin, at mean sea level
  double t2;                    // kelvin

  if (height <= 11) t2 = t1 - 6.5 * height;
  else t2 = 216.65;

  double Rd = 287.053072047065;        // J kg-1 K-1
  double gval = 9.80665;               // m s-2
  double gamma = 0.0065;               // K m-1
  double expon = gval / (Rd * gamma);
  double p1 = 101325;                  // Pa = N m-2 = kg m-1 s-2
                                       // at mean sea level
  double p2 = p1 * pow( t2 / t1, expon);

  // Ideal gas law
  // p v = n r t
  // density = n r / v = p / t
  double density1 = 1.2250;            // kg m-3, at mean sea level
  double density2 = density1 * (p2 / p1) * (t1 / t2);
  return density2;
} // end calcDensity


//======================================================================

// Calc verification deltas and statistics (if testing),
// and to write a text format of the cellMat to the outTxt file.

void Fractl::checkVerif(
  long imain,
  Cell ***& cellMat)           // we set Cell.uu, vv
{

  Statistic statVerifDist;
  Statistic statTruew;
  Statistic statTruev;
  Statistic statTrueu;
  Statistic statCalcw;
  Statistic statCalcv;
  Statistic statCalcu;
  Statistic statDiffw;
  Statistic statDiffv;
  Statistic statDiffu;
  Statistic statAbsDiffw;
  Statistic statAbsDiffv;
  Statistic statAbsDiffu;

  ofstream * ostm = new ofstream();
  ostm->open( outTxt.c_str());

  for (long iz = 0; iz < nradz; iz++) {
    for (long iy = 0; iy < nrady; iy++) {
      for (long ix = 0; ix < nradx; ix++) {

        double centerLocZ = zgridmin + iz * zgridinc;
        double centerLocY = ygridmin + iy * ygridinc;
        double centerLocX = xgridmin + ix * xgridinc;

        bool showDetail = testDetail(
          centerLocZ,           // z
          centerLocY,           // y
          centerLocX);           // x

        if (showDetail) {
          cout << endl << "===== checkVerif: showDetail"
            << "  centerLoc: z: " << centerLocZ
            << "  y: " << centerLocY
            << "  x: " << centerLocX << endl << endl;
        }

        // Get the true wind values
        double verVels[3];    // W, V, U

        if (testMode == Params::MODE_ALPHA || testMode == Params::MODE_BETA) {
          // Get verVelW, verVelV, verVelU
          double synThetaRad = 0;
          double synElevRad = 0;
          calcSyntheticWinds(
            showDetail,
            centerLocZ,         // locz
            centerLocY,         // locy
            centerLocX,         // locx
            verVels);           // returned W, V, U
        }

        else if (testMode == Params::MODE_ZETA_BELTRAMI) {
          calcBeltramiFlow(
            showDetail,
            centerLocZ,         // locz
            centerLocY,         // locy
            centerLocX,         // locx
            verVels);           // returned W, V, U

        }

        else if (testMode == Params::MODE_ZETA) {
          for (int ii = 0; ii < 3; ii++) {
            verVels[ii] = numeric_limits<double>::quiet_NaN();
          }
        }

        else if (testMode == Params::MODE_GAMMA) {
          verVels[0] = synWinds[0];
          verVels[1] = synWinds[1];
          verVels[2] = synWinds[2];
        }

        else throwerr("invalid testMode");

        double wdiff = cellMat[iz][iy][ix].ww - verVels[0];   // W
        double vdiff = cellMat[iz][iy][ix].vv - verVels[1];   // V
        double udiff = cellMat[iz][iy][ix].uu - verVels[2];   // U

        statTruew.addOb( verVels[0]);    // W
        statTruev.addOb( verVels[1]);    // V
        statTrueu.addOb( verVels[2]);    // U

        statCalcw.addOb( cellMat[iz][iy][ix].ww);
        statCalcv.addOb( cellMat[iz][iy][ix].vv);
        statCalcu.addOb( cellMat[iz][iy][ix].uu);

        statDiffw.addOb( wdiff);
        statDiffv.addOb( vdiff);
        statDiffu.addOb( udiff);

        statAbsDiffw.addOb( fabs( wdiff));
        statAbsDiffv.addOb( fabs( vdiff));
        statAbsDiffu.addOb( fabs( udiff));

        bool okFlag = false;
        if ( isOkDouble( cellMat[iz][iy][ix].ww)
          && isOkDouble( cellMat[iz][iy][ix].vv)
          && isOkDouble( cellMat[iz][iy][ix].uu))
        {
          okFlag = true;
        }

        Cell * pcell = & cellMat[iz][iy][ix];
        ostringstream msgstm;
        msgstm << setprecision(7);
        msgstm << "ckv:ok: " << okFlag
          << "  izyx: " << iz << "  " << iy << "  " << ix
          << "  locZYX:"
          << "  " << centerLocZ
          << "  " << centerLocY
          << "  " << centerLocX
          << "  verifWVU:"
          << "  " << verVels[0]
          << "  " << verVels[1]
          << "  " << verVels[2]
          << "  calcWVU:"
          << "  " << pcell->ww
          << "  " << pcell->vv
          << "  " << pcell->uu
          << "  diffWVU:"
          << "  " << (pcell->ww - verVels[0])
          << "  " << (pcell->vv - verVels[1])
          << "  " << (pcell->uu - verVels[2])
          << "  std:"
          << "  " << pcell->ustd
          << "  " << pcell->vstd
          << "  meanNbrElevDeg: " << pcell->meanNbrElevDeg
          << "  meanNbrKeepDist: " << pcell->meanNbrKeepDist
          << "  meanNbrOmitDist: " << pcell->meanNbrOmitDist
          << "  condNum: " << pcell->conditionNumber
          << endl;
        (*ostm) << msgstm.str();
        if (showDetail) cout << msgstm.str();

      } // for ix
    } // for iy
  } // for iz


  cout << setprecision(7);
  cout << endl;
  cout << "checkVerif: imain: " << imain << endl;
  cout << "statVerifDist: ";  statVerifDist.print( 7);
  cout << endl;
  cout << "statTruew: ";  statTruew.print( 7);
  cout << "statTruev: ";  statTruev.print( 7);
  cout << "statTrueu: ";  statTrueu.print( 7);
  cout << endl;
  cout << "statCalcw: ";  statCalcw.print( 7);
  cout << "statCalcv: ";  statCalcv.print( 7);
  cout << "statCalcu: ";  statCalcu.print( 7);
  cout << endl;
  cout << "statDiffw: ";  statDiffw.print( 7);
  cout << "statDiffv: ";  statDiffv.print( 7);
  cout << "statDiffu: ";  statDiffu.print( 7);
  cout << endl;
  cout << "statAbsDiffw: ";  statAbsDiffw.print( 7);
  cout << "statAbsDiffv: ";  statAbsDiffv.print( 7);
  cout << "statAbsDiffu: ";  statAbsDiffu.print( 7);
  cout << endl;

  ostm->close();
  delete ostm;


} // end checkVerif

//==================================================================

// Return a list of filenames in the given dir that
// match the specified regex.

vector<FileSpec *>* Fractl::readDir(
  string dirName,
  string fileRegex)
{
  vector<FileSpec *>* resvec = new vector<FileSpec *>();
  DIR * fdir = opendir( dirName.c_str());
  if (fdir == NULL) badparms("cannot open inDir");
  while (true) {
    struct dirent * dent = readdir( fdir);
    if (dent == NULL) break;
    string fname = dent->d_name;

    regex_t regexWk;
    if (0 != regcomp( &regexWk, fileRegex.c_str(), REG_EXTENDED))
      badparms("invalid fileRegex");
    if (0 == regexec( &regexWk, fname.c_str(), 0, NULL, 0)) {
      string fpath = dirName + "/" + fname;
      struct stat filestat;
      if (0 == stat( fpath.c_str(), &filestat)) {
        if (S_ISREG( filestat.st_mode)) {
          resvec->push_back( new FileSpec( fpath));
        }
      }
    }
    regfree(&regexWk);
  }
  return resvec;
}


//==================================================================


// Read a text file with one file spec per line.
// Return a vector of FileSpec.
// A # in the first column starts a comment.
// Each line has the format
//   fileName altKmMsl latDeg lonDeg

vector<FileSpec *>* Fractl::readFileList(
  string fileListName)
{
  vector<FileSpec *>* resvec = new vector<FileSpec *>();

  ifstream istm( fileListName.c_str());
  if ( ! istm.is_open() ) throwerr("fileList file not found");
  string inLine;

  while (true) {
    getline( istm, inLine, '\n');
    if (istm.eof()) break;
    if (inLine.size() > 0 && inLine.at(0) != '#') {
      vector<string> tokens;
      splitString( inLine, " \t", tokens);
      if (tokens.size() != 4) throwerr("wrong num toks");

      string fpath = tokens[0];
      double altKmMsl = parseDouble( "altKmMsl", tokens.at(1));
      double latDeg = parseDouble( "latDeg", tokens.at(2));
      double lonDeg = parseDouble( "lonDeg", tokens.at(3));

      resvec->push_back( new FileSpec( fpath, altKmMsl, latDeg, lonDeg));
    }
  }
  return resvec;
} // end readFileList

//==================================================================


//==================================================================


// Check that the grid specification is valid.
// This is called 3 times: for z, y, and x grids.

void Fractl::checkGrid(
  const char * msg,
  long num,
  double gridmin,
  double gridmax,
  double gridinc)
{
  string errmsg = "";
  if (num <= 0) errmsg = "grid num <= 0";
  if (fabs( gridmin - gridinc * (gridmin / gridinc)) > epsilon)
    errmsg = "grid min is not a multiple of gridinc";
  if (fabs( gridmax - gridinc * (gridmax / gridinc)) > epsilon)
    errmsg = "grid max is not a multiple of gridinc";
  long ntest = lround( (gridmax - gridmin) / gridinc) + 1;
  if (ntest != num) errmsg = "grid num mismatch";
  if (errmsg != "") {
    cout << setprecision(15);
    cout << endl;
    cout << "grid mismatch.  msg: " << msg << endl;
    cout << "  epsilon: " << epsilon << endl;
    cout << "  num: " << num << endl;
    cout << "  gridmin: " << gridmin << endl;
    cout << "  gridmax: " << gridmax << endl;
    cout << "  gridinc: " << gridinc << endl;
    cout << "  ntest: " << ntest << endl;
    throwerr("grid mismatch");
  }
} // end checkGrid


//==================================================================


// Given a user specified specification triple synWinds,
// calculate the synthetic winds at the specified point.
// The synWinds triple contains: (zspec, yspec, xspec).
// If the spec > SYNFUNC_START (which is way negative),
// just set the wind to the spec value.
// Otherwise the spec is some function id like SYNFUNC_FUNZ,
// which means windComponent = sin(locz).

void Fractl::calcSyntheticWinds(
  bool showDetail,
   double locz,          // wind location
  double locy,
  double locx,
  double * vels)        // returned W, V, U
{
  for (int ii = 0; ii < 3; ii++) {
    double prm = synWinds[ii];
    if (prm > SYNFUNC_START) vels[ii] = prm;
    else if (prm == SYNFUNC_SINZ) vels[ii] = sin( locz);
    else if (prm == SYNFUNC_SINY) vels[ii] = sin( locy);
    else if (prm == SYNFUNC_SINX) vels[ii] = sin( locx);
    else throwerr("unknown synWinds");
  }
  if (showDetail) {
    cout << setprecision(5);
    cout << "      calcsyn:" << endl
      << "        locz: " << locz << "  locy: " << locy
      << "  locx: " << locx << endl
      << "        velW: " << vels[0] << "  velV: " << vels[1]
      << "  velU: " << vels[2] << endl;
  }
} // end calcSyntheticWinds

//==================================================================

// Beltrami flow equations
// See publications by Michael Bell

void Fractl::calcBeltramiFlow(
  bool showDetail,
  double zz,            // location
  double yy,
  double xx,
  double * vels)        // returned W, V, U
{

  // Beltrami flow equations
                                    // Michael Bell's notation:
  double ubase = 10;                // U: mean U wind
  double vbase = 10;                // V: mean V wind
  double wbase = 10;                // A: peak vertical velocity
  double mc = 2 * M_PI / 32000;     // m: vertical wavelength
  double kc = 2 * M_PI / 16000;     // k: horizontal wavelength
  double lc = kc;                   // l

  double amp = wbase / (kc*kc + lc*lc);
  double wavenum = sqrt(kc*kc + lc*lc + mc*mc);
  double tm = 0.;     // time?
  double nu = 15.11e-6;

  zz *= 1000.;         // z in meters
  yy *= 1000.;         // y in meters
  xx *= 1000.;         // x in meters

  double uu = ubase
    - amp
      * exp(-nu * wavenum * wavenum * tm)
      * ( wavenum * lc * cos( kc * (xx - ubase*tm))
            * sin( lc * (yy-vbase*tm))
            * sin(mc*zz)
          + mc * kc * sin( kc*(xx-ubase*tm))
            * cos(lc * (yy-vbase*tm))*cos(mc*zz));

  double dudx = - amp
    * exp( -nu*wavenum*wavenum*tm)
    * ( - wavenum * kc * lc * sin( kc * (xx - ubase*tm))
      * sin( lc*(yy-vbase*tm))
      * sin( mc*zz)
    + mc * kc * kc *cos( kc*(xx-ubase*tm))
      * cos( lc * (yy-vbase*tm)) * cos( mc*zz));

  double dudy = - amp
    * exp( -nu*wavenum*wavenum*tm)
    * ( wavenum * lc * cos(kc * (xx - ubase*tm))
        * lc * cos( lc*(yy-vbase*tm)) * sin( mc*zz)
      - mc * lc * kc * sin(kc*(xx-ubase*tm))
        * sin( lc * (yy-vbase*tm)) * cos(mc*zz));

  double vv = vbase
    + amp
      * exp(-nu * wavenum * wavenum * tm)
      * ( wavenum * kc * sin( kc * (xx - ubase*tm))
            * cos( lc * (yy-vbase*tm))
            * sin( mc * zz)
          - mc * lc * cos( kc*(xx-ubase*tm))
            * sin( lc*(yy-vbase*tm))*cos(mc*zz));

  double dvdx = amp
    * exp( -nu*wavenum*wavenum*tm)
    * ( wavenum * kc * kc * cos( kc*(xx - ubase*tm))
          * cos( lc*(yy-vbase*tm)) * sin( mc*zz)
        + mc * lc * kc * sin( kc*(xx-ubase*tm))
          * sin( lc * (yy-vbase*tm)) * cos( mc*zz));

  double dvdy = amp
    * exp( -nu*wavenum*wavenum*tm)
    * ( -wavenum * kc * lc * sin( kc * (xx - ubase*tm))
      * sin( lc * (yy-vbase*tm)) * sin( mc*zz)
    - mc * lc * lc * cos( kc*(xx-ubase*tm))
      * cos( lc*(yy-vbase*tm)) * cos(mc*zz));

  double ww = wbase
    * cos( kc*(xx-ubase*tm))
    * cos( lc*(yy-vbase*tm))
    * sin(mc*zz)
    * exp( -nu*wavenum*wavenum*tm);

  double vort = 1e5 * (dvdx - dudy);
  double div = 1e5 * (dudx + dvdy);

  vels[0] = ww;
  vels[1] = vv;
  vels[2] = uu;

} // end calcBeltramiFlow


//==================================================================


// Given synthetic W, V, U, calculate the radial velocity
// the radar would have observed.

double Fractl::calcRadialVelocity(
  bool showDetail,
  double thetaRad,      // polar coord angle from observer, radians
  double elevRad,       // elevation angle from observer, radians
  double * vels)        // W, V, U
{
  double velW = vels[0];
  double velV = vels[1];
  double velU = vels[2];

  // Calc radial velocity
  // The view angle can be expressed as a vector having length 1, as
  //   view = (x=cos(theta)*cos(elev), y=sin(theta)*cos(elev), z=sin(elev))
  //
  // The radial velocity is the dot product
  //   view dot (U, V, W)

  double velRadial
    =   cos(thetaRad) * cos(elevRad) * velU
      + sin(thetaRad) * cos(elevRad) * velV
      + sin(elevRad) * velW;

  if (showDetail) {
    cout << setprecision(5);
    cout << "      calcRadial:" << endl
      << "        thetaRad: " << thetaRad
      << "  thetaDeg: " << (thetaRad * 180 / M_PI) << endl
      << "        elevRad: " << elevRad
      << "  elevDeg: " << (elevRad * 180 / M_PI) << endl
      << "        velW: " << vels[0] << "  velV: " << vels[1]
      << "  velU: " << vels[2] << endl
      << "        velRadial: " << velRadial
      << endl;
  }
  return velRadial;
} // end calcRadialVelocity

//==================================================================


// Convert lat,lon coords to Y,X using our transverse mercator projection.

void Fractl::latLonToYX(
  double projLon0Deg,     // central meridian of projection
  double base_y,           // coord y base, km
  double base_x,           // coord x base, km
  double latDeg,          // latitude
  double lonDeg,          // longitude
  double & coordy,        // output coord y, km
  double & coordx)        // output coord x, km
{
  tranMerc.Forward(
    projLon0Deg,
    latDeg,
    lonDeg,
    coordx,             // output value
    coordy);            // output value
  coordy = 0.001 * coordy - base_y;      // convert meters to km
  coordx = 0.001 * coordx - base_x;      // convert meters to km
}

//==================================================================

// Convert Y, X to lat,lon coords using our transverse mercator projection.

void Fractl::yxToLatLon(
  double projLon0,        // central meridian of projection
  double coordy,          // coord y, km
  double coordx,          // coord x, km
  double & latDeg,        // output latitude
  double & lonDeg)        // output longitude
{
  tranMerc.Reverse(
    projLon0,
    1000 * (basex + coordx),
    1000 * (basey + coordy),
    latDeg,             // output value
    lonDeg);            // output value
}


//==================================================================

// Return true if (z,y,x) is within the region specified
// by detailSpec.  The detailSpec contains:
//   [0]: zcenter
//   [1]: ycenter
//   [2]: xcenter
//   [3]: radius about the center

bool Fractl::testDetail(
  double zloc,
  double yloc,
  double xloc)
{
  bool showDetail = false;
  if (detailSpec != NULL
    && zloc >= detailSpec[0] - detailSpec[3]
    && zloc <= detailSpec[0] + detailSpec[3]
    && yloc >= detailSpec[1] - detailSpec[3]
    && yloc <= detailSpec[1] + detailSpec[3]
    && xloc >= detailSpec[2] - detailSpec[3]
    && xloc <= detailSpec[2] + detailSpec[3])
  {
    showDetail = true;
  }
  return showDetail;
}

//==================================================================

// Not used.
// Someday we will use this for testing.

long Fractl::getLongRandom(
  long vmin,
  long vlim)
{
//  int32_t rr = 0;
//  if (0 != random_r( randInfo, &rr)) throwerr("random_r error");
  std::random_device r;
  std::default_random_engine e1(r());
  std::uniform_int_distribution<int> uniform_dist(vmin, vlim);
  long res = uniform_dist(e1);
//  double frac = rr / (((double) RAND_MAX) + 1);0
//  long res = vmin + (long) (frac * (vlim - vmin));
  if (res < 0 || res >= vlim) throwerr("invalid getRandLong res");
  return res;
}

//==================================================================

// Return true if val is a normal value, not NaN or +/-inf.

bool Fractl::isOkDouble( double val) {
  bool res = false;
  if (val == 0 || isnormal( val)) res = true;
  return res;
}

//==================================================================

// Return true if val is a normal value, not NaN or +/-inf.

bool Fractl::isOkFloat( float val) {
  bool res = false;
  if (val == 0 || isnormal( val)) res = true;
  return res;
}

//==================================================================

// Print the elapsed run time since the previous call, in seconds.

void Fractl::printRunTime(
  const string& str,
  struct timeval * ptva)
{
  struct timeval tvb;
  if (gettimeofday( &tvb, NULL) != 0) throwerr("gettimeofday err");
  double deltaSec = tvb.tv_sec - ptva->tv_sec
    + 1.e-6 * (tvb.tv_usec - ptva->tv_usec);
  cout << "runTime: " << setw(20) << str << ": " << deltaSec << endl;
  ptva->tv_sec = tvb.tv_sec;
  ptva->tv_usec = tvb.tv_usec;
}

//==================================================================

// Add the elapsed run time since the previous call, in seconds.

void Fractl::addDeltaTime(
  struct timeval * ptva,
  double * psum)
{
  struct timeval tvb;
  if (gettimeofday( &tvb, NULL) != 0) throwerr("gettimeofday err");
  double deltaSec = tvb.tv_sec - ptva->tv_sec
    + 1.e-6 * (tvb.tv_usec - ptva->tv_usec);
  ptva->tv_sec = tvb.tv_sec;
  ptva->tv_usec = tvb.tv_usec;
  if (psum != NULL) (*psum) += deltaSec;
}

//==================================================================

// Split a string into a vector of tokens.
// Why does C++ make this so weird?

void Fractl::splitString(
  const string& str,
  const string& delimiters,
  vector<string>& tokens)       // appended
{
  // Scan to find beginning of first token.
  long begPos = str.find_first_not_of(delimiters, 0);

  while (begPos != string::npos) {
    // Find end of this token
    long endPos = str.find_first_of(delimiters, begPos);
    tokens.push_back( str.substr( begPos, endPos - begPos));
    // Find start of next token
    begPos = str.find_first_not_of( delimiters, endPos);
  }
}

//==================================================================

// format a time, stored as double seconds since 1970.

string Fractl::formatTime( double dtm)
{
  time_t itime = (time_t) dtm;
  struct tm stm;
  gmtime_r( &itime, &stm);

  char bufa[1000];
  strftime( bufa, 1000, "%Y-%m-%dT%H:%M:%S", &stm);

  long ifrac = 1000 * (dtm - itime);
  char bufb[100];
  snprintf( bufb, 100, "%03d", ifrac);
  string stg = bufa;
  stg += ".";
  stg += bufb;
  return stg;
}

//==================================================================

// Converts string to long.

long Fractl::parseLong( string msg, string stg)
{
  char * endptr;
  const char * stgc = stg.c_str();
  long ires = strtol( stgc, &endptr, 10);
  if (endptr != stgc + std::strlen( stgc)) throwerr("invalid integer");
  return ires;
}


//==================================================================

// Converts string to double.

double Fractl::parseDouble( string msg, string stg)
{
  char * endptr;
  const char * stgc = stg.c_str();
  double dres = strtod( stgc, &endptr);
  if (endptr != stgc + std::strlen( stgc)) throwerr("invalid number");
  return dres;
}

//==================================================================

// Converts string (y or n) to bool.

bool Fractl::parseBool( string msg, string stg)
{
  bool bres = false;
  if (stg == "n") bres = false;
  else if (stg == "y") bres = true;
  else throwerr("invalid bool");
  return bres;
}

//==================================================================

// Throws an exception

void Fractl::throwerr( const char * msg, ...)
{
  int nbufa = 10000;
  char bufa[10000];

  va_list arglist;
  va_start( arglist, msg);
  vsnprintf( bufa, nbufa, msg, arglist);
  va_end( arglist);

  cerr << "Fractl: throwerr: " << bufa << endl;
  cerr.flush();
  throw bufa;
}

//==================================================================


//==================================================================

void Fractl::dumpArgs()
{
  cout << setprecision(15);
  cout << "bugs: " << bugs << endl;
  cout << "testMode: " << testMode << endl;
  cout << "zgrid: " << zgridmin << "," << zgridmax << "," << zgridinc << endl;
  cout << "ygrid: " << ygridmin << "," << ygridmax << "," << ygridinc << endl;
  cout << "xgrid: " << xgridmin << "," << xgridmax << "," << xgridinc << endl;
  cout << "projName: " << projName << endl;
  cout << "projLat0: " << projLat0 << endl;
  cout << "projLon0: " << projLon0 << endl;
  cout << "baseW: " << baseW << endl;
  cout << "epsilon: " << epsilon << endl;
  cout << "maxDeltaAltKm: " << maxDeltaAltKm << endl;
  cout << "maxAbsElevDeg: " << maxAbsElevDeg << endl;
  cout << "minRadialDistKm: " << minRadialDistKm << endl;
  cout << "numNbrMax: " << numNbrMax << endl;
  cout << "maxDistBase: " << maxDistBase << endl;
  cout << "maxDistFactor: " << maxDistFactor << endl;
  cout << "forceOk: " << forceOk << endl;
  cout << "useEigen: " << useEigen << endl;
  cout << "inDir: " << inDir << endl;
  cout << "fileRegex: " << fileRegex << endl;
  cout << "fileList: " << fileList << endl;
  cout << "radialName: " << radialName << endl;
  cout << "dbzName: " << dbzName << endl;
  cout << "ncpName: " << ncpName << endl;
  cout << "outTxt: " << outTxt << endl;
  cout << "outNc: " << outNc << endl;
  if (detailSpec == NULL)
    cout << "detailSpec: (none)" << endl;
  else {
    cout << "detailSpec:"
      << "  z: " << detailSpec[0] << "  y: " << detailSpec[1]
      << "  x: " << detailSpec[2] << "  delta: " << detailSpec[3] << endl;
  }
}

//==================================================================

bool Fractl::checkArgs()
{
  if (testMode == Params::MODE_NONE) badparms("parameter not specified: -testMode");
  if (synWinds == NULL) badparms("parameter not specified: -synWinds");
  if (radFiles == NULL) badparms("parameter not specified: -radFiles");
  if (zgridinc == MISS_PARM) badparms("parm not specified: -zgrid");
  if (ygridinc == MISS_PARM) badparms("parm not specified: -ygrid");
  if (xgridinc == MISS_PARM) badparms("parm not specified: -xgrid");
  if (projName == "") badparms("parm not specified: -projName");
  if (projLat0 == MISS_PARM) badparms("parm not specified: -projLat0");
  if (projLon0 == MISS_PARM) badparms("parm not specified: -projLon0");
  if (baseW == MISS_PARM) badparms("parm not specified: -baseW");
  if (epsilon == MISS_PARM) badparms("parm not specified: -epsilon");
  if (maxDeltaAltKm == MISS_PARM)
    badparms("parm not specified: -maxDeltaAltKm");
  if (maxAbsElevDeg == MISS_PARM)
    badparms("parm not specified: -maxAbsElevDeg");
  if (minRadialDistKm == MISS_PARM)
    badparms("parm not specified: -minRadialDistKm");
  if (numNbrMax == MISS_PARM) badparms("parm not specified: -numNbrMax");
  if (maxDistBase == MISS_PARM) badparms("parm not specified: -maxDistBase");
  if (maxDistFactor == MISS_PARM)
    badparms("parm not specified: -maxDistFactor");
  if (testMode == Params::MODE_ZETA && forceOk)
    throwerr("forceOk not allowed with testMode zeta");

  if (inDir == "") {
    if (fileList == "") badparms("must spec either -inDir or -fileList");
    if (fileRegex != "") badparms("parm -inDir required for -fileRegex");
  }
  else {
    if (fileList != "") badparms("cannot spec both -inDir and -fileList");
    if (fileRegex == "") badparms("parm -fileRegex required with -inDir");
  }

  if (radialName == "") badparms("parameter not specified: -radialName");
  if (dbzName == "") badparms("parameter not specified: -dbzName");
  if (ncpName == "") badparms("parameter not specified: -ncpName");
  if (outTxt == "") badparms("parameter not specified: -outTxt");
  if (outNc == "") badparms("parameter not specified: -outNc");

  if (projName != "transverseMercator") badparms("unknown projName");
  printRunTime("init", &timea);

  return true;
}

//==================================================================

bool Fractl::loadObservations()
{
  bool retVal;

  if (testMode == Params::MODE_ALPHA)
    retVal = fillWithSyntheticWinds();
  else
    retVal = fillWithObservations();

  if(verbose) {
    cout << "cntTimeh: " << cntTimeh << "  sumTimeh: " << sumTimeh << endl;
    cout << "cntTimei: " << cntTimei << "  sumTimei: " << sumTimei << endl;
    cout << "cntTimej: " << cntTimej << "  sumTimej: " << sumTimej << endl;
    cout << "cntTimek: " << cntTimek << "  sumTimek: " << sumTimek << endl;
    cout << "cntTimel: " << cntTimel << "  sumTimel: " << sumTimel << endl;
    cout << "cntTimem: " << cntTimem << "  sumTimem: " << sumTimem << endl;
    cout << "cntTimen: " << cntTimen << "  sumTimen: " << sumTimen << endl;
  }
  return retVal;
}


//==================================================================

bool Fractl::fillWithSyntheticWinds()
{
  // Fill pointVec with all the observations.
  // Find the overall bounding box.

  if (zgridmin == MISS_PARM) badparms("parm not specified: -zgrid");
  if (ygridmin == MISS_PARM) badparms("parm not specified: -ygrid");
  if (xgridmin == MISS_PARM) badparms("parm not specified: -xgrid");
  if (zgridmax == MISS_PARM) badparms("parm not specified: -zgrid");
  if (ygridmax == MISS_PARM) badparms("parm not specified: -ygrid");
  if (xgridmax == MISS_PARM) badparms("parm not specified: -xgrid");
  if (zgridinc == MISS_PARM) badparms("parm not specified: -zgrid");
  if (ygridinc == MISS_PARM) badparms("parm not specified: -ygrid");
  if (xgridinc == MISS_PARM) badparms("parm not specified: -xgrid");
  nradz = lround( (zgridmax - zgridmin) / zgridinc) + 1;
  nrady = lround( (ygridmax - ygridmin) / ygridinc) + 1;
  nradx = lround( (xgridmax - xgridmin) / xgridinc) + 1;
  // checkGrid( "radz", epsilon, nradz, zgridmin, zgridmax, zgridinc);
  // checkGrid( "rady", epsilon, nrady, ygridmin, ygridmax, ygridinc);
  // checkGrid( "radx", epsilon, nradx, xgridmin, xgridmax, xgridinc);
  checkGrid( "radz", nradz, zgridmin, zgridmax, zgridinc);
  checkGrid( "rady", nrady, ygridmin, ygridmax, ygridinc);
  checkGrid( "radx", nradx, xgridmin, xgridmax, xgridinc);

  pointVec = new vector<Point *>();

  for (long iz = 0; iz < nradz; iz++) {
    for (long iy = 0; iy < nrady; iy++) {
      for (long ix = 0; ix < nradx; ix++) {

	double centerLoc[3];
	centerLoc[0] = zgridmin + iz * zgridinc;
	centerLoc[1] = ygridmin + iy * ygridinc;
	centerLoc[2] = xgridmin + ix * xgridinc;
	bool showDetail = testDetail(
				     centerLoc[0],         // z
				     centerLoc[1],         // y
				     centerLoc[2]);         // x
	//detailSpec);          // z, y, x, delta
	if (showDetail) {
	  cout << "const.init: showDetail:"
	       << "  iz: " << iz << "  iy: " << iy << "  ix: " << ix << endl;
	  cout << "const: showDetail: centerLoc z: " << centerLoc[0] << endl;
	  cout << "const: showDetail: centerLoc y: " << centerLoc[1] << endl;
	  cout << "const: showDetail: centerLoc x: " << centerLoc[2] << endl;
	}

	aircraftBbox->addOb( centerLoc[0], centerLoc[1], centerLoc[2]);
	pointBbox->addOb( centerLoc[0], centerLoc[1], centerLoc[2]);

	// Make numNbrMax points, each one exactly at the cell center.
	for (long inbr = 0; inbr < numNbrMax; inbr++) {
	  // We must use different thetas so the problem has a solution.
	  double synThetaRad = inbr * M_PI / 4;
	  double synElevRad = 0;
	  double synVels[3];    // W, V, U
	  calcSyntheticWinds(
			     showDetail,
			     // synWinds,           // user specified winds
			     centerLoc[0],       // locz
			     centerLoc[1],       // locy
			     centerLoc[2],       // locx
			     synVels);           // returned W, V, U

	  double velRadial = calcRadialVelocity(
						showDetail,
						synThetaRad,        // polar coord angle from observer, radians
						synElevRad,         // elevation angle from observer, radians
						synVels);           // W, V, U

	  pointVec->push_back( new Point(
					 0,             // ifile
					 0,             // iray
					 0,             // ipt
					 centerLoc[0],  // aircraftz
					 centerLoc[1],  // aircrafty
					 centerLoc[2],  // aircraftx
					 centerLoc[0],  // altKmMsl
					 centerLoc[1],  // coordy
					 centerLoc[2],  // coordx
					 1329336556,    // UTC seconds since 1970 jan 1.
					 synThetaRad,   // angle from aircraft in horiz plane from north
					 synElevRad,    // angle from the aircraft between horiz and pt
					 velRadial,     // radial velocity
					 30.,           // dbz: reflectivity
					 1.0));         // valNcp: radar net coherent power
	} // for inbr
      } // for ix
    } // for iy
  } // for iz
  return true;
}

//==================================================================

bool Fractl::fillWithObservations()
{
  // Fill pointVec with all the observations.
  // Find the overall bounding box.

  // Read dir to get list of radar file names.
  vector<FileSpec *>* fspecList;
  if (inDir != "")
    fspecList = readDir( inDir, fileRegex);
  else
    fspecList = readFileList( fileList);
  printRunTime("readDir", &timea);
  cout << endl << "fsubsetList:" << endl;

  // Get a list of files to load

  vector<FileSpec *>* fsubsetList = new vector<FileSpec *>();
  for (long ifile = 0; ifile < fspecList->size(); ifile++) {
    FileSpec * fspec = fspecList->at( ifile);
    if (radFiles[0] == 0 && radFiles[1] == 0
        || ifile >= radFiles[0] && ifile < radFiles[1])
      {
        fsubsetList->push_back( fspec);
	if (bugs >= Params::DEBUG_VERBOSE)
	  cout << "  fsubsetList: ifile: " << ifile
	       << "  fspec: " << fspec->fpath << endl;
      }
  }

  if (fsubsetList->size() <= 0)
    throwerr("No input file found in %s", inDir.c_str());

  // Fill pointVec with all the observations from all files.
  // Find the overall bounding box.
  pointVec = new vector<Point *>();

  double maxAbsErrHoriz = 0;
  double maxAbsErrVert = 0;

  Statistic statVg;
  Statistic statDbz;
  Statistic statNcp;

  for (long ifile = 0; ifile < fsubsetList->size(); ifile++) {
    FileSpec * fspec = fsubsetList->at( ifile);
    cout << "main: start fpath: " << fspec->fpath << endl;

    bool fOk = true;

    if (preGridded)
      fOk = readPreGriddedFile(
		  ifile,
		  fspec,
		  aircraftBbox,
		  pointBbox,
		  &statVg,
		  &statDbz,
		  &statNcp,
		  pointVec,                 // appended
		  &maxAbsErrHoriz,
		  &maxAbsErrVert,
		  timeMin,                  // returned
		  timeMax);                 // returned
    else
      readRadarFile(
		  forceOk,
		  ifile,
		  fspec,
		  maxAbsElevDeg,
		  aircraftBbox,
		  pointBbox,
		  &statVg,
		  &statDbz,
		  &statNcp,
		  pointVec,                 // appended
		  &maxAbsErrHoriz,
		  &maxAbsErrVert,
		  timeMin,                  // returned
		  timeMax);                 // returned

    if ( ! fOk ) {
      std::cerr << "Issue reading " << fspec->fpath << " aborting." << std::endl;
      return false;
    }
  } // for ifile
  delete fsubsetList;
  printRunTime("readRadarFiles", &timea);

  // xxx del maxAbsErrHoriz, vert
  cout << setprecision(5);
  cout << "final maxAbsErrHoriz: " << maxAbsErrHoriz << endl;
  cout << "final maxAbsErrVert: " << maxAbsErrVert << endl;

  cout << endl;
  cout << "aircraftBbox:" << endl;
  aircraftBbox->print();

  cout << endl;
  cout << "pointBbox:" << endl;
  pointBbox->print();
  cout << endl;

  cout << endl;
  cout << "statVg:" << endl;
  statVg.print( 7);
  cout << endl;

  cout << endl;
  cout << "statNcp:" << endl;
  statNcp.print( 7);
  cout << endl;

  cout << endl;
  cout << setprecision(20);
  cout << "timeMin: " << timeMin << endl;
  cout << "timeMax: " << timeMax << endl;
  cout << endl;
  return true;
}

//==================================================================

bool Fractl::buildKdTree()
{
  // Build the KD tree for 3D radar observations
  long npt = pointVec->size();
  cout << "npt: " << npt << endl;

  KD_real **radarKdMat = new KD_real*[npt];     // npt x ndim

  for (long ii = 0; ii < npt; ii++) {
    radarKdMat[ii] = new KD_real[ndim];
    Point * pt = pointVec->at(ii);
    radarKdMat[ii][0] = pt->coordz;        // z
    radarKdMat[ii][1] = pt->coordy;        // y
    radarKdMat[ii][2] = pt->coordx;        // x
    if (ii % 1000 == 0) {
      cout << setprecision(5);
    }
  }
  printRunTime("build KD mat", &timea);

  radarKdTree = new KD_tree(
    (const KD_real **) radarKdMat,
    npt,
    ndim);
  printRunTime("build KD tree", &timea);
  return true;
}

//==================================================================

bool Fractl::findLimits()
{
  // If the limits for zgrid,ygrid,xgrid were not specified,
  // then try to get them from the bounding box.
  if (pointBbox != NULL) {
    if (zgridmin == MISS_PARM) {
      zgridmin = zgridinc * floor( pointBbox->coordzMin / zgridinc);
      zgridmax = zgridinc * ceil( pointBbox->coordzMax / zgridinc);
    }
    if (ygridmin == MISS_PARM) {
      ygridmin = ygridinc * floor( pointBbox->coordyMin / ygridinc);
      ygridmax = ygridinc * ceil( pointBbox->coordyMax / ygridinc);
    }
    if (xgridmin == MISS_PARM) {
      xgridmin = xgridinc * floor( pointBbox->coordxMin / xgridinc);
      xgridmax = xgridinc * ceil( pointBbox->coordxMax / xgridinc);
    }
  } // if pointBbox != NULL

  if (zgridmin == MISS_PARM) badparms("parm not specified: -zgrid");
  if (ygridmin == MISS_PARM) badparms("parm not specified: -ygrid");
  if (xgridmin == MISS_PARM) badparms("parm not specified: -xgrid");

  // Check the grids and calc nradz, nrady, nradx.
  nradz = lround( (zgridmax - zgridmin) / zgridinc) + 1;
  nrady = lround( (ygridmax - ygridmin) / ygridinc) + 1;
  nradx = lround( (xgridmax - xgridmin) / xgridinc) + 1;

  checkGrid( "radz", nradz, zgridmin, zgridmax, zgridinc);
  checkGrid( "rady", nrady, ygridmin, ygridmax, ygridinc);
  checkGrid( "radx", nradx, xgridmin, xgridmax, xgridinc);

  cout << setprecision(5);
  cout << "main:" << endl;
  cout << "radz: nradz: " << nradz
       << "  zgridmin: " << zgridmin
       << "  zgridmax: " << zgridmax
       << "  zgridinc: " << zgridinc
       << endl;
  cout << "rady: nrady: " << nrady
       << "  ygridmin: " << ygridmin
       << "  ygridmax: " << ygridmax
       << "  ygridinc: " << ygridinc
       << endl;
  cout << "radx: nradx: " << nradx
       << "  xgridmin: " << xgridmin
       << "  xgridmax: " << xgridmax
       << "  xgridinc: " << xgridinc
       << endl;
}

//==================================================================

bool Fractl::allocateCellMat()
{
  //xxx all new: delete
  // Allocate cellMat.
  // Init Cell.ww = 0.

  if (gridType == Params::GRID_MESH) {
    cellMat = new Cell**[nradz];
    for (long iz = 0; iz < nradz; iz++) {
      cellMat[iz] = new Cell*[nrady];
      for (long iy = 0; iy < nrady; iy++) {
	cellMat[iz][iy] = new Cell[nradx];
	for (long ix = 0; ix < nradx; ix++) {
	  cellMat[iz][iy][ix].ww = 0;
	}
      }
    }
  } else {
    long maxx = (xgridmax - xgridmin) * 2 / xgridinc + 4;
    long maxy = (ygridmax - ygridmin) * 2 / ygridinc + 4;  
    long maxz = (zgridmax - zgridmin) * 2 / zgridinc + 4;

    cellMat = new Cell**[maxx];
    for (long iz = 0; iz < maxz; iz++) {
      cellMat[iz] = new Cell*[maxy];
      for (long iy = 0; iy < maxy; iy++) {
	cellMat[iz][iy] = new Cell[maxx];
	for (long ix = 0; ix < maxx; ix++) {
	  cellMat[iz][iy][ix].ww = 0;
	}
      }
    }    
  }
  
  printRunTime("alloc cellMat", &timea);
  return true;
}

//==================================================================

void Fractl::freeCellMat()
{
  long maxy, maxz;
  
  if (gridType == Params::GRID_MESH) {
    maxy = nrady;
    maxz = nradz;
  } else {
    maxy = (ygridmax - ygridmin) * 2 / ygridinc + 4;
    maxz = (zgridmax - zgridmin) * 2 / zgridinc + 4;
  }
  for (long iz = 0; iz < maxz; iz++) {
    for (long iy = 0; iy < maxy; iy++) {
      delete[] cellMat[iz][iy];
    }
    delete[] cellMat[iz];
  }
  
  delete[] cellMat;
  cellMat = NULL;
}
