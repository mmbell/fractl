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

#ifndef CELL_H
#define CELL_H
#include <limits>

class Cell {

public:

double ww;               // w wind
double vv;               // v wind
double uu;               // u wind

double meanNbrDbz;       // mean(dbz) over neighbors
double meanNbrNcp;       // mean(ncp) over neighbors
double meanNbrElevDeg;   // mean( elevation angle degrees) over neighbors
double meanNbrKeepDist;  // mean( distance from cell center) of keep nbrs
double meanNbrOmitDist;  // mean( distance from cell center) of omit nbrs
double conditionNumber;

double ustd;             // variance in U from geometry
double vstd;             // variance in V from geometry
double wstd;             // variance in W from geometry

Cell()
{
  this->ww = std::numeric_limits<double>::quiet_NaN();
  this->vv = std::numeric_limits<double>::quiet_NaN();
  this->uu = std::numeric_limits<double>::quiet_NaN();

  this->meanNbrDbz = std::numeric_limits<double>::quiet_NaN();
  this->meanNbrNcp = std::numeric_limits<double>::quiet_NaN();
  this->meanNbrElevDeg = std::numeric_limits<double>::quiet_NaN();
  this->meanNbrKeepDist = std::numeric_limits<double>::quiet_NaN();
  this->meanNbrOmitDist = std::numeric_limits<double>::quiet_NaN();
  this->conditionNumber = std::numeric_limits<double>::quiet_NaN();

  this->ustd = std::numeric_limits<double>::quiet_NaN();
  this->vstd = std::numeric_limits<double>::quiet_NaN();
  this->wstd = std::numeric_limits<double>::quiet_NaN();
} // end constructor


}; // end class

#endif
