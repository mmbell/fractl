#include "Utils.hh"

// These are used by the Cedric Leise filter 

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
  
  for(int z = 0; z < zdim; ++z)  
    for(int y = 0; y < ydim; ++y)
      for(int x = 0; x < xdim; ++x)
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

// These are used by the Cedric mass continuity integration

double *cell2plane_U(Cell *** &cellMat, int xdim, int ydim, int level)
{
  double *data = new double[xdim * ydim];
  int index = 0;

  // Put data into a Fortran 2 dimentional array
  
    for(int y = 0; y < ydim; ++y)
      for(int x = 0; x < xdim; ++x)
	data[index++] = cellMat[level][y][x].uu;
  return data;
}

void plane2cell_U(Cell *** &cellMat, int xdim, int ydim, int level, double *data)
{
  int index = 0;

  // Put data back into cell matrix
  
  for(int y = 0; y < ydim; ++y)
    for(int x = 0; x < xdim; ++x)
      cellMat[level][y][x].uu = data[index++];
}

double *cell2plane_V(Cell *** &cellMat, int xdim, int ydim, int level)
{
  double *data = new double[xdim * ydim];
  int index = 0;

  // Put data into a Fortran 2 dimentional array
  
    for(int y = 0; y < ydim; ++y)
      for(int x = 0; x < xdim; ++x)
	data[index++] = cellMat[level][y][x].vv;
  return data;
}

void plane2cell_V(Cell *** &cellMat, int xdim, int ydim, int level, double *data)
{
  int index = 0;

  // Put data back into cell matrix
  
  for(int y = 0; y < ydim; ++y)
    for(int x = 0; x < xdim; ++x)
      cellMat[level][y][x].vv = data[index++];
}

double *cell2plane_EU(Cell *** &cellMat, int xdim, int ydim, int level)
{
  double *data = new double[xdim * ydim];
  int index = 0;

  // Put data into a Fortran 2 dimentional array
  
    for(int y = 0; y < ydim; ++y)
      for(int x = 0; x < xdim; ++x)
	data[index++] = cellMat[level][y][x].ustd;
  return data;
}

void plane2cell_EU(Cell *** &cellMat, int xdim, int ydim, int level, double *data)
{
  int index = 0;

  // Put data back into cell matrix
  
  for(int y = 0; y < ydim; ++y)
    for(int x = 0; x < xdim; ++x)
      cellMat[level][y][x].ustd = data[index++];
}
double *cell2plane_EV(Cell *** &cellMat, int xdim, int ydim, int level)
{
  double *data = new double[xdim * ydim];
  int index = 0;

  // Put data into a Fortran 2 dimentional array
  
    for(int y = 0; y < ydim; ++y)
      for(int x = 0; x < xdim; ++x)
	data[index++] = cellMat[level][y][x].vstd;
  return data;
}

void plane2cell_EV(Cell *** &cellMat, int xdim, int ydim, int level, double *data)
{
  int index = 0;

  // Put data back into cell matrix
  
  for(int y = 0; y < ydim; ++y)
    for(int x = 0; x < xdim; ++x)
      cellMat[level][y][x].vstd = data[index++];
}

double *cell2plane_Dbz(Cell *** &cellMat, int xdim, int ydim, int level)
{
  double *data = new double[xdim * ydim];
  int index = 0;

  // Put data into a Fortran 2 dimentional array
  
    for(int y = 0; y < ydim; ++y)
      for(int x = 0; x < xdim; ++x)
	data[index++] = cellMat[level][y][x].meanNbrDbz;
  return data;
}

void plane2cell_W(Cell *** &cellMat, int xdim, int ydim, int level, double *data)
{
  int index = 0;

  // Put data back into cell matrix
  
  for(int y = 0; y < ydim; ++y)
    for(int x = 0; x < xdim; ++x)
      cellMat[level][y][x].ww = data[index++];
}

void fill_array(double *dest, std::size_t size, double val)
{
  for(std::size_t index = 0; index < size; index++)
    *dest++ = val;
}

void copy_array(double *dest, double *src, std::size_t size)
{
  for(std::size_t index = 0; index < size; index++)
    *dest++ = *src++;
}

