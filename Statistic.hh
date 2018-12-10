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

#ifndef STATISTIC_H
#define STATISTIC_H

#include <iomanip>

class Statistic {

public:

long numTot;
long numNaN;
long numNegInf;
long numPosInf;
long numGood;

double dmin;
double dmax;
double dsum;
double dsumsq;


Statistic() {
  numTot = 0;
  numNaN = 0;
  numNegInf = 0;
  numPosInf = 0;
  numGood = 0;
  dmin   = std::numeric_limits<double>::quiet_NaN();
  dmax   = std::numeric_limits<double>::quiet_NaN();
  dsum   = 0;
  dsumsq = 0;
} // end constructor



void addOb(
  double val)
{
  numTot++;
  if (std::isnan(val)) numNaN++;
  else if (std::isinf(val) ) numNegInf++;
  else if (std::isinf(val) ) numPosInf++;
  else {
    numGood++;
    if (std::isnan( dmin) || val < dmin) dmin = val;
    if (std::isnan( dmax) || val > dmax) dmax = val;
    dsum += val;
    dsumsq += val * val;
  }
} // end addOb



void print( int prec) {
  int width = prec + 7;
  int svprec = std::cout.precision( 15);
  std::cout
    << "Statistic: "
    << "  numTot: " << numTot;

  std::cout << std::setprecision( prec);
  if (numGood > 0) {
    std::cout << "  min: " << std::setw( width) << dmin;
    std::cout << "  max: " << std::setw( width) << dmax;
    std::cout << "  mean: " << std::setw( width) << (dsum / numGood);
  }

  if (numGood > 1) {
    double variance = (dsumsq - dsum * dsum / numGood) / (numGood-1);
    double stddev = 0;
    if (variance > 0) stddev = sqrt( variance);
    std::cout << "  stddev: " << std::setw( width) << stddev;
  }

  std::cout << std::setprecision( 15);
  if (numGood != numTot) std::cout << "  numGood: " << numGood;
  if (numNaN != 0) std::cout << "  numNaN: " << numNaN;
  if (numNegInf != 0) std::cout << "  numNegInf: " << numNegInf;
  if (numPosInf != 0) std::cout << "  numPosInf: " << numPosInf;
  std::cout << std::endl;
  std::cout.precision( svprec);
}




///void throwerr( const char * msg) {
///  std::cout << "Statistic: throwerr: " << msg << std::endl;
///  std::cout.flush();
///  throw msg;
///}


}; // end class Statistic

#endif


