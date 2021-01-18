#include <iostream>
#include "Filters.hh"
#include "Utils.hh"

// The filters operate on an array of float
// So extract the right field into an array, apply the filter, then shove the updated
// values back in the right field

// TODO DRY this code

bool Filter::filter_U(Cell *** &cellMat, int xdim, int ydim, int zdim, int nstep, int *ns)
{
  
  std::cout << "Starting Leise filtering on U" << std::endl;
  
  double *data = cell2array_U(cellMat, xdim, ydim, zdim);
  bool retval = filter(data, xdim, ydim, zdim, nstep, ns);
  array2cell_U(cellMat, xdim, ydim, zdim, data);
  delete[] data;
  
  std::cout << "Done with Leise filtering on U" << std::endl;
  
  return retval;
}

bool Filter::filter_V(Cell *** &cellMat, int xdim, int ydim, int zdim, int nstep, int *ns)
{
  
  std::cout << "Starting Leise filtering on V" << std::endl;
  
  double *data = cell2array_V(cellMat, xdim, ydim, zdim);
  bool retval = filter(data, xdim, ydim, zdim, nstep, ns);
  array2cell_V(cellMat, xdim, ydim, zdim, data);
  delete[] data;
  
  std::cout << "Done with Leise filtering on V" << std::endl;
  
  return retval;
}

bool Filter::filter_W(Cell *** &cellMat, int xdim, int ydim, int zdim, int nstep, int *ns)
{
  
  std::cout << "Starting Leise filtering on W" << std::endl;
  
  double *data = cell2array_W(cellMat, xdim, ydim, zdim);
  bool retval = filter(data, xdim, ydim, zdim, nstep, ns);
  array2cell_W(cellMat, xdim, ydim, zdim, data);
  delete[] data;
  
  std::cout << "Done with Leise filtering on W" << std::endl;
  
  return retval;
}

// -------------------- Filter Factory


Filter *FilterFactory::createFilter(Params::filter_t t)
{
  switch(t) {
  case Params::FILTER_LEISE:
    return new LeiseFilter();
  case Params::FILTER_NONE:
    return NULL;
  default:
    std::cerr << "Unsupported Filter Type" << std::endl;
    return NULL;
  }
}
