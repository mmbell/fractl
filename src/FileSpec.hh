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

#ifndef FILESPEC_H
#define FILESPEC_H

#include <string>
#include <limits>
#include <iostream>

class FileSpec {

public:

std::string fpath;         // file path
double altitudeKmMsl;
double latitudeDeg;
double longitudeDeg;



FileSpec( std::string ppath)
{
  this->fpath = ppath;
  this->altitudeKmMsl = std::numeric_limits<double>::quiet_NaN();
  this->latitudeDeg   = std::numeric_limits<double>::quiet_NaN();
  this->longitudeDeg  = std::numeric_limits<double>::quiet_NaN();
} // end constructor



FileSpec(
  std::string ppath,
  double altKmMsl,
  double latDeg,
  double lonDeg)
{
  this->fpath = ppath;
  this->altitudeKmMsl = altKmMsl;
  this->latitudeDeg = latDeg;
  this->longitudeDeg = lonDeg;
} // end constructor





void print() {
  std::cout
    << "      fpath: " << fpath << std::endl
    << "      altitudeKmMsl: " << altitudeKmMsl << std::endl
    << "      latitudeDeg: " << latitudeDeg << std::endl
    << "      longitudeDeg: " << longitudeDeg << std::endl;
}

}; // end class

#endif
