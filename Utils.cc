#include "Utils.hh"

double *cell2array_U(Cell *** &cellMat, int xdim, int ydim, int zdim)
{
  double *data = new double[xdim * ydim * zdim];
  int index = 0;

  // Put data into a one dimensional array
  
  for(int z = 0; z < zdim; ++z)  
    for(int y = 0; y < ydim; ++y)
      for(int x = 0; x < xdim; ++x)
	data[index++] = cellMat[z][y][x].uu;
  return data;
}

void array2cell_U(Cell *** &cellMat, int xdim, int ydim, int zdim, double *data)
{
  int index = 0;

  // Put data back into cell matrix
  
  for(int z = 0; z < z; ++z)  
    for(int y = 0; y < y; ++y)
      for(int x = 0; x < x; ++x)
	cellMat[z][y][x].uu = 	data[index++];
}

double *cell2array_V(Cell *** &cellMat, int xdim, int ydim, int zdim)
{
  double *data = new double[xdim * ydim * zdim];
  int index = 0;

  // Put data into a one dimensional array
  
  for(int z = 0; z < zdim; ++z)  
    for(int y = 0; y < ydim; ++y)
      for(int x = 0; x < xdim; ++x)  
	data[index++] = cellMat[z][y][x].vv;
  return data;
}

void array2cell_V(Cell *** &cellMat, int xdim, int ydim, int zdim, double *data)
{
  int index = 0;

  // Put data back into cell matrix
  for(int z = 0; z < zdim; ++z)  
    for(int y = 0; y < ydim; ++y)
      for(int x = 0; x < xdim; ++x)  
	cellMat[z][y][x].vv = 	data[index++];
}

double *cell2array_W(Cell *** &cellMat, int xdim, int ydim, int zdim)
{
  double *data = new double[xdim * ydim * zdim];
  int index = 0;

  // Put data into a one dimensional array
  
  for(int z = 0; z < zdim; ++z)
    for(int y = 0; y < ydim; ++y)
      for(int x = 0; x < xdim; ++x)
	data[index++] = cellMat[z][y][x].ww;
  return data;
}

void array2cell_W(Cell *** &cellMat, int xdim, int ydim, int zdim, double *data)
{
  int index = 0;

  // Put data back into cell matrix
  
  for(int z = 0; z < zdim; ++z)
    for(int y = 0; y < ydim; ++y)
      for(int x = 0; x < xdim; ++x)
	cellMat[z][y][x].ww = data[index++];
}
