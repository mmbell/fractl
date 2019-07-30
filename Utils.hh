#include <cstddef>

#include "Cell.hh"

double *cell2array_U(Cell *** &cellMat, int xdim, int ydim, int zdim);
void array2cell_U(Cell *** &cellMat, int xdim, int ydim, int zdim, double *data);

double *cell2array_V(Cell *** &cellMat, int xdim, int ydim, int zdim);
void array2cell_V(Cell *** &cellMat, int xdim, int ydim, int zdim, double *data);

double *cell2array_W(Cell *** &cellMat, int xdim, int ydim, int zdim);
void array2cell_W(Cell *** &cellMat, int xdim, int ydim, int zdim, double *data);

double *cell2plane_U(Cell *** &cellMat, int xdim, int ydim, int level);
void plane2cell_U(Cell *** &cellMat, int xdim, int ydim, int level, double *data);
  
double *cell2plane_V(Cell *** &cellMat, int xdim, int ydim, int level);
void plane2cell_V(Cell *** &cellMat, int xdim, int ydim, int level, double *data);

double *cell2plane_W(Cell *** &cellMat, int xdim, int ydim, int level);
void plane2cell_W(Cell *** &cellMat, int xdim, int ydim, int level, double *data);

double *cell2plane_EU(Cell *** &cellMat, int xdim, int ydim, int level);
void plane2cell_EU(Cell *** &cellMat, int xdim, int ydim, int level, double *data);

double *cell2plane_EV(Cell *** &cellMat, int xdim, int ydim, int level);
void plane2cell_EV(Cell *** &cellMat, int xdim, int ydim, int level, double *data);

double *cell2plane_Dbz(Cell *** &cellMat, int xdim, int ydim, int level);

void fill_array(double *dest, std::size_t size, double val);
void copy_array(double *dest, double *src, std::size_t size);
