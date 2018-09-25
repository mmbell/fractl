#include <cstddef>
#include "Cell.hh"
#include "Params.hh"

class Filter {
  
public:
  
  Filter() {};
  virtual ~Filter() {};

  bool filter_U(Cell *** &cellMat, int xdim, int ydim, int zdim, int nstep, int *ns = NULL);
  bool filter_V(Cell *** &cellMat, int xdim, int ydim, int zdim, int nstep, int *ns = NULL);
  bool filter_W(Cell *** &cellMat, int xdim, int ydim, int zdim, int nstep, int *ns = NULL);  

  virtual bool filter(double *data, int xdim, int ydim, int zdim, int nstep, int *ns = NULL) = 0;
  
protected:
  
  int max(int i1, int i2) { return (i1 > i2) ? i1 : i2; };
  int min(int i1, int i2) { return (i1 < i2) ? i1 : i2; };
  
};

class LeiseFilter : public Filter {
  
public:
  
  LeiseFilter();
  ~LeiseFilter();
  bool filter(double *data, int xdim, int ydim, int zdim, int nstep, int *ns = NULL);
};


class FilterFactory {
  
public:
  
  static Filter *createFilter(Params::filter_t t);
};
