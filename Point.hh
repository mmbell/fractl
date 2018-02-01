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

#ifndef POINT_H
#define POINT_H

class Point {

public:

  long ifile;
  long iray;
  long ipt;
  double aircraftz;
  double aircrafty;
  double aircraftx;
  double coordz;
  double coordy;
  double coordx;
  double rayTime;      // UTC seconds since 1970 jan 1.
  double thetaRad;     // polar coords angle, radians
  double elevRad;      // angle above horiz plane, radians
  double vg;           // radial velocity
  double dbz;          // reflectivity, decibels
  double ncp;          // radar net coherent power

  double dist;         // distance from desired location; used for sorting


  Point(
	long ifile,
	long iray,
	long ipt,
	double aircraftz,
	double aircrafty,
	double aircraftx,
	double coordz,
	double coordy,
	double coordx,
	double rayTime,
	double thetaRad,
	double elevRad,
	double vg,
	double dbz,
	double ncp)
  {
    this->ifile = ifile;
    this->iray = iray;
    this->ipt = ipt;
    this->aircraftz = aircraftz;
    this->aircrafty = aircrafty;
    this->aircraftx = aircraftx;
    this->coordz = coordz;
    this->coordy = coordy;
    this->coordx = coordx;
    this->rayTime = rayTime;
    this->thetaRad = thetaRad;
    this->elevRad = elevRad;
    this->vg = vg;
    this->dbz = dbz;
    this->ncp = ncp;
  } // end constructor

  // Copy constructor
  Point(
	Point *pt)
  {
    this->ifile = pt->ifile;
    this->iray = pt->iray;
    this->ipt = pt->ipt;
    this->aircraftz = pt->aircraftz;
    this->aircrafty = pt->aircrafty;
    this->aircraftx = pt->aircraftx;
    this->coordz = pt->coordz;
    this->coordy = pt->coordy;
    this->coordx = pt->coordx;
    this->rayTime = pt->rayTime;
    this->thetaRad = pt->thetaRad;
    this->elevRad = pt->elevRad;
    this->vg = pt->vg;
    this->dbz = pt->dbz;
    this->ncp = pt->ncp;
  } // end constructor

  void print() {
    std::cout
      << "      ifile: " << ifile << std::endl
      << "      iray: " << iray << std::endl
      << "      ipt: " << ipt << std::endl
      << "      aircraft_zyx: " << aircraftz
      << "  " << aircrafty << "  " << aircraftx << std::endl
      << "      coord_zyx: " << coordz
      << "  " << coordy << "  " << coordx << std::endl
      << "      rayTime: " << rayTime << std::endl
      << "      theta deg: " << (thetaRad * 180 / M_PI) << std::endl
      << "      elev deg: " << (elevRad * 180 / M_PI) << std::endl
      << "      vg: " << vg << std::endl
      << "      dbz: " << dbz << std::endl
      << "      ncp: " << ncp << std::endl;
  }

}; // end class

#endif
