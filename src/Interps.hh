#include <cstddef>
#include "Cell.hh"
#include "Params.hh"

class Interp {
  
public:
  
  Interp() {};
  virtual ~Interp() {};

  bool interpolate_U(Cell *** &cellMat, int xdim, int ydim, int zdim, double flag);
  bool interpolate_V(Cell *** &cellMat, int xdim, int ydim, int zdim, double flag);

  virtual bool interpolate(double *data, int xdim, int ydim, int zdim, double flag) = 0;
  
protected:
  
  int max(int i1, int i2) { return (i1 > i2) ? i1 : i2; };
  int min(int i1, int i2) { return (i1 < i2) ? i1 : i2; };
  
};

// Use the Interpolation method written by Jim Leise in Cedric

class LeiseInterp : public Interp {
  
public:
  LeiseInterp();
  ~LeiseInterp();
  bool interpolate(double *data, int xdim, int ydim, int zdim, double flag);
};

// Use the Interpolation method that came with the original lrose RadarWind

class RadarWindInterp : public Interp {
  
public:
  
  RadarWindInterp();
  ~RadarWindInterp();
  bool interpolate(double *data, int xdim, int ydim, int zdim, double flag);
};

// The Interp Factory
  
class InterpFactory {
  
public:
  
  static Interp *createInterp(Params::interp_t t);
};
