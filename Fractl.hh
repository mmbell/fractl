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

#ifndef FRACTL_H
#define FRACTL_H

#include <vector>
#include <iostream>

#include "Bbox.hh"
#include "Statistic.hh"
#include "Cell.hh"
#include "Point.hh"
#include "FileSpec.hh"

// See kdtree source at:
//   cvs/libs/kd/src/include/kd/kd.hh
//   cvs/libs/kd/src/kd/kd.cc

#include <kd/kd.hh>

// GeographicLib.  See doc at:
//   http://geographiclib.sourceforge.net
//   http://geographiclib.sourceforge.net/html/organization.html
//   http://geographiclib.sourceforge.net/html/geodesic.html
//   http://geographiclib.sourceforge.net/html/classGeographicLib_1_1TransverseMercatorExact.html

#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/TransverseMercatorExact.hpp>

// TODO revisit these globals.

extern int cntTimeh;
extern double sumTimeh;
extern int cntTimei;
extern double sumTimei;
extern int cntTimej;
extern double sumTimej;
extern int cntTimek;
extern double sumTimek;
extern int cntTimel;
extern double sumTimel;
extern int cntTimem;
extern double sumTimem;
extern int cntTimen;
extern double sumTimen;

using namespace std;

class Fractl {

public:

  static const int HANDLE_BBOX      = 1;
  static const int HANDLE_GETSTATS  = 2;

  static const int MISS_PARM = -99999;

#if 0
  static const int TESTMODE_ALPHA          = 1;
  static const int TESTMODE_BETA           = 2;
  static const int TESTMODE_ZETA           = 3;
  static const int TESTMODE_ZETA_BELTRAMI  = 4;
#endif

  static const long SYNFUNC_START = -999900;
  static const long SYNFUNC_SINX  = -999901;
  static const long SYNFUNC_SINY  = -999902;
  static const long SYNFUNC_SINZ  = -999903;

  void badparms( const string msg, ...);

  Fractl();
  ~Fractl();

  bool run(int argc, char *argv[]);

  bool parseArgs(int argc, char *argv[]);
  bool parseFile(const char *fileName, int &argc, char ***argv);

  // TODO Make these private
  bool parseSynWinds(char *param, double *synWinds);
  bool parseRadFiles(char *param, long *radFiles);
  bool parseGridSpec(char *param, double &min, double &max, double &step);
  bool parseDetailSpec(char *param, double *spec);

  bool checkArgs();
  void dumpArgs();

  bool loadObservations();
  bool fillWithSyntheticWinds();
  bool fillWithObservations();

  bool readPreGriddedFile(long ifile,
			  FileSpec *fspec,

                     Bbox * aircraftBbox,              // overall bounding box of aircraft locs
		     Bbox * pointBbox,                 // overall bounding box of point locs
		     Statistic * statVg,               // statistics for radial velocity
		     Statistic * statDbz,              // statistics for DBZ
		     Statistic * statNcp,              // statistics for net coherent power
		     std::vector<Point *> *pointVec,
		     double *maxAbsErrHoriz,
		     double *maxAbsErrVert,
		     double &timeMin,                  // returned
                     double &timeMax);


  void readRadarFile(
		     long forceOk,        // if != 0, force all point ncp and dbz to ok
		     long ifile,
		     FileSpec * fspec,
		     double maxElevDeg,
		     Bbox * aircraftBbox,              // overall bounding box of aircraft locs
		     Bbox * pointBbox,                 // overall bounding box of point locs
		     Statistic * statVg,               // statistics for radial velocity
		     Statistic * statDbz,              // statistics for DBZ
		     Statistic * statNcp,              // statistics for net coherent power
		     std::vector<Point *> *pointVec,
		     double *maxAbsErrHoriz,
		     double *maxAbsErrVert,
		     double &timeMin,                  // returned
		     double &timeMax);                  // returned

  void calcAllVU(
		 long maxNumNbr,              // max num nearest nbrs
		 std::vector<Point *> *pointVec,   // all observations
		 KD_tree * radarKdTree,       // nearest nbr tree for pointVec
		 Cell *** & cellMat);          // we set Cell.uu, vv

  void calcCellVU(
		  KD_real * centerLoc,           // query point
		  long maxNumNbr,                // max num nearest nbrs
		  std::vector<Point *> *pointVec,     // all observations
		  KD_tree * radarKdTree,         // nearest nbr tree for pointVec
		  Cell * pcell);                  // we fill vv, uu.

  void calcAllW(
		Cell *** & cellMat);          // we set Cell.ww

  double calcDensity( double height);


  void checkVerif(
		  long imain,
		  Cell *** & cellMat);

  std::vector<FileSpec *>* readDir(
			      string dirName,
			      string filePrefix);


  std::vector<FileSpec *>* readFileList(
				   string fileListName);


  void writeNetcdf();

  void checkGrid(
		 const char * msg,
		 long num,
		 double minVal,
		 double maxVal,
		 double incVal);


  void calcSyntheticWinds(
			  bool showDetail,
			  double locz,          // wind location
			  double locy,
			  double locx,
			  double * vels);       // returned W, V, U


  void calcBeltramiFlow(
			bool showDetail,
			double zz,            // location
			double yy,
			double xx,
			double * vels);       // returned W, V, U


  double calcRadialVelocity(
			    bool showDetail,
			    double thetaRad,      // polar coord angle from observer, radians
			    double elevRad,       // elevation angle from observer, radians
			    double * vels);       // W, V, U


  void latLonToYX(
		  double projLon0Deg,     // central meridian of projection
		  double base_y,           // coord y base, km
		  double base_x,           // coord x base, km
		  double latDeg,          // latitude
		  double lonDeg,          // longitude
		  double & coordy,        // output coord y, km
		  double & coordx);       // output coord x, km


  void yxToLatLon(
		  double projLon0Deg,     // central meridian of projection
		  double coordy,          // coord y, km
		  double coordx,          // coord x, km
		  double & latDeg,        // output latitude
		  double & lonDeg);       // output longitude


  double calcDistLocPt( double * loc, Point * pta);

  double calcDistPtPt( Point * pta, Point * ptb);

  double calcDistPtAircraft( Point * pta);


  bool testDetail(
		  double zloc,
		  double yloc,
		  double xloc);


  long getLongRandom(
		     long vmin,
		     long vlim);

  bool isOkDouble( double val);


  bool isOkFloat( float val);


  void printRunTime(
		    const string& str,
		    struct timeval * ptva);


  void addDeltaTime(
		    struct timeval * ptva,
		    double * psum);


  void splitString(
		   const string& str,
		   const string& delimiters,
		   std::vector<string>& tokens);


  string formatTime( double dtm);

  long parseLong( string msg, string stg);

  double parseDouble( string msg, string stg);

  bool parseBool( string msg, string stg);

  void throwerr( const char * msg, ...);

  bool buildKdTree();
  bool findLimits();
  bool allocateCellMat();
  void freeCellMat();
  bool calcWinds();

private:

  // These are filled when parsing arguments

  bool preGridded;

  long bugs;
  int testMode;
  double * synWinds;
  long * radFiles;

  double zgridmin;
  double zgridmax;
  double zgridinc;

  double ygridmin;
  double ygridmax;
  double ygridinc;

  double xgridmin;
  double xgridmax;
  double xgridinc;

  std::string projName;
  double projLat0;
  double projLon0;

  double baseW;
  double epsilon;
  double maxDeltaAltKm;
  double maxAbsElevDeg;
  double minRadialDistKm;
  long numNbrMax;
  double maxDistBase;
  double maxDistFactor;
  bool forceOk;
  bool useEigen;
  //xxx also spec max dist of observation from cell center
  std::string inDir;
  std::string fileRegex;
  std::string fileList;
  std::string radialName;
  std::string dbzName;
  std::string ncpName;
  std::string outTxt;
  std::string outNc;
  double * detailSpec;

  double earthRadiusMeter;
  double flattening;

  double basex;
  double basey;

  double radarAlt;

  long nradx;
  long nrady;
  long nradz;

  int ndim;

  bool verbose;

  struct timeval timea;

  std::vector<Point *> *pointVec;
  KD_tree * radarKdTree;

  Bbox * aircraftBbox;   // overall bounding box of aircraft locs
  Bbox * pointBbox;      // overall bounding box of point locs
  double timeMin;
  double timeMax;

  GeographicLib::Geodesic *geodesic;
  static const GeographicLib::TransverseMercatorExact tranMerc;

  double minDbz;
  double minNcp;

  Cell *** cellMat;

}; // end class

#endif
