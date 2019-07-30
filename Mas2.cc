// Port of Cedric MS2DRV.F function

#include <iostream>

#include <math.h>
#include <string>

#include "Utils.hh"
#include "Fractl.hh"
#include "Utils.hh"

#define CSP_INDEX(a, i, j) (a)[j * 3 + i]
#define DATA_INDEX(a, i, j) (a)[j * nradx + i]

// #define RHOFUN(l, alpha) exp(-alpha * CSP_INDEX(csp, 0, 2) * CSP_INDEX(csp, 2, 2) * (l - 1))
// #define RHOFUN(l, alpha) exp(-alpha * CSP_INDEX(csp, 0, 2) * CSP_INDEX(csp, 2, 2) * (l))

#define RHOFUN(l, alpha) exp((-l) * (alpha))

#define HEIGHT(l) (CSP_INDEX(csp, 0, 2) + CSP_INDEX(csp, 2, 2) * (l - 1))

double rhofun(double *csp, int level, double alpha)
{
  double csp_1_3 = CSP_INDEX(csp, 0, 2);
  double csp_3_3 = CSP_INDEX(csp, 2, 2);
  double retval = exp(-alpha * csp_1_3 * csp_3_3 * (double) (level - 1));
  std::cout << "rhofun at level " << level << ": " << retval << std::endl;
  std::cout << "exp( -" << alpha << " * " << csp_1_3 << " * " << csp_3_3
	    << " * (" << level << " - 1))" << std::endl;
  return retval;
}

void dumpBuf(double *buf, size_t nradx, size_t nrady)
{
  int count = 0;
  double min = 10000;
  double max = -10000;
  
  for(size_t i = 0; i < nradx; i++)
    for(size_t j = 0; j < nrady; j++)
      if (! std::isnan(DATA_INDEX(buf, i, j)) ) {
	count++;
	if ( DATA_INDEX(buf, i, j) < min)
	  min = DATA_INDEX(buf, i, j);
	if ( DATA_INDEX(buf, i, j) > max)
	 max = DATA_INDEX(buf, i, j);
      }
			   
  std::cout << "values: " << count
	    << ", min: " << min
	    << ", max: " << max
	    << ", nan: " << nradx * nrady - count
	    << std::endl;
}

int Fractl::ms2drv() {
  
  // Initialization

  double c2rd[2] = { mas2ErrMax, mas2IterMax };

  size_t plane_size = nradx * nrady;
  double *xydeli;
  double *csp;
  int iactc = 0;
    
  initMas2(&csp, &xydeli);
  
  size_t start_level;
  size_t stop_level;
  int increment;

  double *w_current = new double[plane_size];
  double *w_previous = new double[plane_size];  
  double *conv_current = new double[plane_size];
  double *conv_previous = new double[plane_size];    
  double *ewU_previous = new double[plane_size];  
  double *ewV_previous = new double[plane_size];
  
  switch(mas2Dir) {
  case Params::DIR_DOWN:
    start_level = nradz;
    stop_level = 0;
    increment = -1;
  case Params::DIR_UP:
  default:
    start_level = 0;
    stop_level = nradz;
    increment = 1;
  }
  
  double dzh = (0.5 * CSP_INDEX(csp, 2, 2)) * increment;

  // params to dwitter()
  
  int iter;
  double dmn;
  int kst;
  
  // For each level

  for (size_t level = start_level; level != stop_level; level += increment) {

    std::cout << "---------------- level " << level << " ---------------" << std::endl;
    
    float alpha = 0.0;
    // TODO  alpha defaults to 0.1, but get sets according to ???
    if (alpha <= 0.0)
      alpha = 0.1;
    
    double *U = cell2plane_U(cellMat, nradx, nrady, level);
    double *V = cell2plane_V(cellMat, nradx, nrady, level);
    double *Dbz = cell2plane_Dbz(cellMat, nradx, nrady, level);

    double *ewU_current = cell2plane_EU(cellMat, nradx, nrady, level);
    double *ewV_current = cell2plane_EV(cellMat, nradx, nrady, level);
    
    // double zlev = HEIGHT(level);
    double zlev = level * zgridinc;
    double avg_den = vtA * exp(vtC * alpha * (zlev + zgridinc * increment));
    
    double rhop = RHOFUN(zlev, alpha);
    double rhoc = RHOFUN(zlev - zgridinc * increment, alpha);    

    std::cout << "---- Level " << level << ", zlev: " << zlev
	      << ", rhoc: " << rhoc << ", rhop: "<< rhop
	      << std::endl;

    std::cout << "-- U before manipulation:" << std::endl;
    dumpBuf(U, nradx, nrady);
    std::cout << "-- V before manipulation:" << std::endl;
    dumpBuf(V, nradx, nrady);
    
    std::cout << "*** ewU_current:" << std::endl;
    dumpBuf(ewU_current, nradx, nrady);
    std::cout << "*** ewV_current" << std::endl;
    dumpBuf(ewV_current, nradx, nrady);

    std::cout << "*** dbz" << std::endl;
    dumpBuf(Dbz, nradx, nrady);

    double vt_min = 1000;
    double vt_max = -1000;
    
    for(size_t index = 0; index < plane_size; index++) {
      if(std::isnan(ewU_current[index]))
	U[index] = numeric_limits<double>::quiet_NaN();
      if(std::isnan(ewV_current[index]))
	V[index] = numeric_limits<double>::quiet_NaN();

#if 0 // TODO
      if( ! std::isnan(Dbz[index]) ) {
	double vt = -avg_den * pow(10.0, (0.1 * vtB * Dbz[index]));
	if (vt < vt_min)
	  vt_min = vt;
	if (vt > vt_max)
	  vt_max = vt;
	
	if( ! std::isnan(U[index]))
	  U[index] += vt * ewU_current[index];
	if( ! std::isnan(V[index]))
	  V[index] += vt * ewV_current[index];
      } else {
	U[index] = numeric_limits<double>::quiet_NaN();
	V[index] = numeric_limits<double>::quiet_NaN();
      }
#endif
    } // index

    std::cout << "vt_min: " << vt_min << ", vt_max: " << vt_max << std::endl;
    
    std::cout << "++ U after manipulation:" << std::endl;
    dumpBuf(U, nradx, nrady);
    std::cout << "++ V after manipulation:" << std::endl;
    dumpBuf(V, nradx, nrady);
    
    pconvg(U, V, conv_current, nradx, nrady, 3, xydeli, csp, iactc);
    std::cout << std::endl << "conv_current after pconvg, level " << level << std::endl;
    dumpBuf(conv_current, nradx, nrady);
    
    fill_array(w_current, plane_size, numeric_limits<double>::quiet_NaN());

    if (level != start_level) {
      dwiter(ewU_current,        // rbuf(1, 3)
	     ewV_current,        // rbuf(1, 4)
	     conv_current,       // rbuf(1, 5)
	     w_current,          // rbuf(1, 6)
	     conv_previous,      // rbuf(1, 7)
	     w_previous,         // rbuf(1, 8)
	     ewU_previous,       // rbuf(1, 9)
	     ewV_previous,       // rbuf(1, 10)
	     nradx, nrady, xydeli,
	     dzh, rhoc, rhop, c2rd, &iter, &dmn, &kst, level, zlev); // TODO
      std::cout << "Retval at level " << level << ": " << kst << std::endl;
      std::cout << "// U after dwitter:" << std::endl;
      dumpBuf(U, nradx, nrady);
      std::cout << "// V after dwitter:" << std::endl;
      dumpBuf(V, nradx, nrady);
      
      std::cout << std::endl << "w_current after dwitter, level " << level << std::endl;
      dumpBuf(w_current, nradx, nrady);
    }

    if (level != stop_level) { // initialize boundary
      switch(lowerBoundaryInitMethod) {
      case Params::INIT_FRACT:
      case Params::INIT_FIELD:
	std::cout << "PARAMS::INIT_FRACT and PARAMS:INIT_FIELD not implemented yet" << std::endl;
      case Params::INIT_CONST:
	double val = std::stod(initStr);
	fill_array(w_previous, plane_size, val);
	break;
      }
      for(int i = 0; i < nradx; i++)
	for(int j = 0; j < nrady; j++) {
	  if ( std::isnan(DATA_INDEX(w_current, i, j)) &&
	       ! std::isnan(DATA_INDEX(conv_current, i, j)))
	    DATA_INDEX(w_current, i, j) = DATA_INDEX(w_previous, i, j);
	}
    }
    std::cout << "## w_current after boundary manip" << std::endl;
    dumpBuf(w_current, nradx, nrady);
    
    // Get everything ready for next level
    copy_array(conv_previous, conv_current, plane_size);
    copy_array(w_previous, w_current, plane_size);
    copy_array(ewU_previous, ewU_current, plane_size);
    copy_array(ewV_previous, ewV_current, plane_size);
  
    for(int i = 0; i < nradx; i++)
      for(int j = 0; j < nrady; j++) {
	if ( ! std::isnan(DATA_INDEX(w_current, i, j)) ) {
	  DATA_INDEX(U, i, j) += DATA_INDEX(w_current, i, j)
	    * DATA_INDEX(ewU_current, i, j);
	  DATA_INDEX(V, i, j) += DATA_INDEX(w_current, i, j)
	    * DATA_INDEX(ewV_current, i, j);
	} else {
	  DATA_INDEX(U, i, j) = numeric_limits<double>::quiet_NaN();
	  DATA_INDEX(V, i, j) = numeric_limits<double>::quiet_NaN();	
	}
      }
      std::cout << "\\ U before saving:" << std::endl;
      dumpBuf(U, nradx, nrady);
      std::cout << "\\ V before saving:" << std::endl;
      dumpBuf(V, nradx, nrady);
  
    // save U, V, w_current into the cell
    
    plane2cell_U(cellMat, nradx, nrady, level, U);
    plane2cell_V(cellMat, nradx, nrady, level, V);
    plane2cell_W(cellMat, nradx, nrady, level, w_current);
    
    delete[] U;
    delete[] V;
    delete[] ewU_current;
    delete[] ewV_current;    
    
  } // each level
  
  return 0; // TODO
}

bool Fractl::initMas2(double **csp, double **xydeli)
{
  *csp = new double[3 * 3];
  *xydeli = new double[2];
  
  // double sfi = 1.0 / 100;	// TODO
  double sfi = 1.0;

  CSP_INDEX(*csp, 0, 0) = xgridmin * sfi;
  CSP_INDEX(*csp, 1, 0) = xgridmax * sfi;
  CSP_INDEX(*csp, 2, 0) = xgridinc * 0.001;
    
  CSP_INDEX(*csp, 0, 1) = ygridmin * sfi;
  CSP_INDEX(*csp, 1, 1) = ygridmax * sfi;
  CSP_INDEX(*csp, 2, 1) = ygridinc * 0.001;

  if (zgridmin == 0.0)
    CSP_INDEX(*csp, 0, 2) = fabs(zgridinc) * sfi; // * 0.1;
  else
    CSP_INDEX(*csp, 0, 2) = zgridmin * sfi; // * 0.1;
  CSP_INDEX(*csp, 1, 2) = zgridmax * sfi; // * 0.1;
  CSP_INDEX(*csp, 2, 2) = zgridinc; // * 0.001;

  (*xydeli)[0] = 1.0 / CSP_INDEX(*csp, 2, 0);
  (*xydeli)[1] = 1.0 / CSP_INDEX(*csp, 2, 1);

  return true;	// TODO
}
