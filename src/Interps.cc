#include <iostream>
#include "Interps.hh"
#include "Utils.hh"

// The filters operate on an array of float
// So extract the right field into an array, apply the filter, then shove the updated
// values back in the right field

bool Interp::interpolate_U(Cell *** &cellMat, int xdim, int ydim, int zdim, double flag)
{
  
  std::cout << "Starting Interpolation on U" << std::endl;
  
  double *data = cell2array_U(cellMat, xdim, ydim, zdim);
  bool retval = interpolate(data, xdim, ydim, zdim, flag);
  array2cell_U(cellMat, xdim, ydim, zdim, data);
  delete[] data;
  
  std::cout << "Done with Interpolation on U" << std::endl;
  
  return retval;
}

bool Interp::interpolate_V(Cell *** &cellMat, int xdim, int ydim, int zdim, double flag)
{
  
  std::cout << "Starting Interpolation on V" << std::endl;
  
  double *data = cell2array_V(cellMat, xdim, ydim, zdim);
  bool retval = interpolate(data, xdim, ydim, zdim, flag);
  array2cell_V(cellMat, xdim, ydim, zdim, data);
  delete[] data;
  
  std::cout << "Done with Interpolation on V" << std::endl;
  
  return retval;
}

// -------------------- Interp Factory


Interp *InterpFactory::createInterp(Params::interp_t t)
{
  switch(t) {
  case Params::INTERP_LEISE:
    return new LeiseInterp();
  case Params::INTERP_RADAR_WIND:
    return new RadarWindInterp();
  case Params::INTERP_NONE:
    return NULL;
  default:
    std::cerr << "Unsupported Interpolation Type" << std::endl;
    return NULL;
  }
}
