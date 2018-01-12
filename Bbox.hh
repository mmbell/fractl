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

#ifndef BBOX_H
#define BBOX_H

#include <limits>
#include <cmath>
class Bbox {

public:


long nob;
double coordzMin;
double coordzMax;
double coordyMin;
double coordyMax;
double coordxMin;
double coordxMax;



Bbox() {
  nob = 0;
  coordzMin = std::numeric_limits<double>::quiet_NaN();
  coordzMax = std::numeric_limits<double>::quiet_NaN();
  coordyMin = std::numeric_limits<double>::quiet_NaN();
  coordyMax = std::numeric_limits<double>::quiet_NaN();
  coordxMin = std::numeric_limits<double>::quiet_NaN();
  coordxMax = std::numeric_limits<double>::quiet_NaN();
} // end constructor




void print() {
  std::cout
    << "bbox: nob: " << nob << std::endl
    << "  coordz: min: " << coordzMin << "  max: " << coordzMax << std::endl
    << "  coordy: min: " << coordyMin << "  max: " << coordyMax << std::endl
    << "  coordx: min: " << coordxMin << "  max: " << coordxMax << std::endl;
  std::cout << std::endl;
}


void addOb(
  double coordz,
  double coordy,
  double coordx)
{
  nob++;
  if (std::isnan( coordzMin) || coordz < coordzMin) coordzMin = coordz;
  if (std::isnan( coordzMax) || coordz > coordzMax) coordzMax = coordz;
  if (std::isnan( coordyMin) || coordy < coordyMin) coordyMin = coordy;
  if (std::isnan( coordyMax) || coordy > coordyMax) coordyMax = coordy;
  if (std::isnan( coordxMin) || coordx < coordxMin) coordxMin = coordx;
  if (std::isnan( coordxMax) || coordx > coordxMax) coordxMax = coordx;
} // end addOb




void throwerr( const char * msg) {
  std::cout << "Bbox: throwerr: " << msg << std::endl;
  std::cout.flush();
  throw msg;
}


}; // end class Bbox

#endif


