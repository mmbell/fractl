#include "Cell.hh"

double *cell2array_U(Cell *** &cellMat, int xdim, int ydim, int zdim);
void array2cell_U(Cell *** &cellMat, int xdim, int ydim, int zdim, double *data);

double *cell2array_V(Cell *** &cellMat, int xdim, int ydim, int zdim);
void array2cell_V(Cell *** &cellMat, int xdim, int ydim, int zdim, double *data);

double *cell2array_W(Cell *** &cellMat, int xdim, int ydim, int zdim);
void array2cell_W(Cell *** &cellMat, int xdim, int ydim, int zdim, double *data);
