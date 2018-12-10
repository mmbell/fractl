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

//
// Read a pre-gridded set of files (the result of Radx2Grid for example)
//
// NOT operational yet. Depends on features not yet available in Radx2Grid

#include <iostream>

#include "Fractl.hh"

// NetCDF: See doc at:
//   netcdfC/tda/netcdf-4.1.3/examples/CXX/pres_temp_4D_wr.cpp

#include <netcdfcpp.h>

// Forward declaration of helper functions

bool getFillValue(NcVar *var, float &val);
bool getDimInfo(NcFile &file, int dim, const char *varName,
		float &spacing, float &min, float &max);
bool getGlobalAttribute(NcFile &file, std::string name, float &val);
bool getGridMapping(NcFile &file, float &radar_lat, float &radar_lon);
bool getOriginLatLon(NcFile &file, float &lat, float &lon);
bool getTime(NcFile &file, const char *varName, double &time);
  
// TODO
// Bbox      (overall bounding box of aircraft locs
// pointBox  (overall bounding box of point locs
// pointVec
// macAbsErr (Horiz and Vert)
// timeMin and timeMax
//

// TODO This isn't working yet (need to figure out what to do with
// needed variables that aren't available in a pre-gridded file.
// (Radar altitude, NCP, ...)

bool Fractl::readPreGriddedFile(
  long              ifile,
  FileSpec        * fspec,
  Bbox            * aircraftBbox,    // overall bounding box of aircraft locs
  Bbox            * pointBbox,       // overall bounding box of point locs
  Statistic       * statVg,          // statistics for radial velocity
  Statistic       * statDbz,         // statistics for DBZ
  Statistic       * statNcp,         // statistics for net coherent power
  vector<Point *> * pointVec,        // appended
  double          * maxAbsErrHoriz,
  double          * maxAbsErrVert,
  double          & timeMin,         // returned
  double          & timeMax)         // returned

{
  NcError ncError(NcError::verbose_nonfatal); // Prevent NertCDF error from exiting the program

  // Open the file
  NcFile file(fspec->fpath.c_str(), NcFile::ReadOnly);

  if (! file.is_valid() ) {
    std::cerr << "Can't read, (or not a netCDF file) " <<  fspec->fpath.c_str()
	      << std::endl;
    return false;
  }

  // Read the dimentions

  NcDim *x0 = file.get_dim("x0");
  NcDim *y0 = file.get_dim("y0");
  NcDim *z0 = file.get_dim("z0");

  // TODO set variables in Fractl object instead

  int iDim = x0->size();
  int jDim = y0->size();
  int kDim = z0->size();

  float iGridSp, xmin, xmax;
  float jGridSp, ymin, ymax;
  float kGridSp, zmin, zmax;

  float radarLat, radarLon, radarAlt;

  double time;

  if ( ! getTime(file, "start_time", time) ) {
    std::cerr << "Can't get start_time from file" << std::endl;
    return false;
  }
  
  // Get grid info. This will set *GridSp, *min, *max in all 3 dimentions
  
  if (! getDimInfo(file, iDim, "x0", iGridSp, xmin, xmax) ) {
    std::cerr << "Can't get x0 array from file" << std::endl;
    return false;
  }

  if (! getDimInfo(file, jDim, "y0", jGridSp, ymin, ymax) ) {
    std::cerr << "Can't get y0 array from file" << std::endl;
    return false;
  }

  if (! getDimInfo(file, kDim, "z0", kGridSp, zmin, zmax) ) {
    std::cerr << "Can't get z0 array from file" << std::endl;
    return false;
  }

#if 0
  // These are non-standard fields.
  // Get radar location

  if ( ! getGlobalAttribute(file, "radar_altitude", radarAlt)) {
    std::cerr << "Can't get radar_altitude from file" << std::endl;
    return false;
  }
  if ( ! getGlobalAttribute(file, "radar_latitude", radarLat)) {
    std::cerr << "Can't get radar_latitude from file" << std::endl;
    return false;
  }
  if ( ! getGlobalAttribute(file, "radar_longitude", radarLon)) {
    std::cerr << "Can't get radar_longitude from file" << std::endl;
    return false;
  }
#endif
  radarAlt = -1.0;
  
  // This gets the origin lat and lon (radar location, usually at 0, 0 kilometers)
  
  if (! getGridMapping(file, radarLat, radarLon) ) {
    std::cerr << "Can't get grid mapping from file" << std::endl;
    return false;
  }

  // Radar altitude must come from the config file as it isn't saved
  // by Radx2Grid, unless a couple of flags are turned on.

  if ( radarAlt < 0) {
    std::cerr << "radarAlt parameter was not set. Setting radar alt to 0.0" << std::endl;
    std::cerr << "   It is a required parameter when using preGridded mode" << std::endl;
    radarAlt = 0.0;
  }

  //TODO. Radar is at 0,0 km. Do we need distance or index here?
  double radarx = 0;
  double radary = 0;
  
  aircraftBbox->addOb(radarAlt, radary, radarx);
  
  float latReference, lonReference;
  // TODO what are these references?
  if (! getOriginLatLon(file, latReference, lonReference) ) {
    std::cerr << "Can't get origin Lat and Lon from file" << std::endl;
    return false;
  }

  // These are the lat and lon at index 0, 0
  
  fspec->latitudeDeg = latReference;
  fspec->longitudeDeg = lonReference;

  // Should this be radar altitude, or altitude at origin of grid?
  
  fspec->altitudeKmMsl = radarAlt;
  
  NcVar *reflectivity = file.get_var(dbzName.c_str());
  if (reflectivity == NULL) {
    std::cerr << "Can't get reflectivity '" << dbzName
	      << "' from " << fspec->fpath << std::endl;
    return false;
  }

  NcVar *velocity = file.get_var(radialName.c_str());
  if (velocity == NULL) {
    std::cerr << "Can't get velocity '" << dbzName
	      << "' from " << fspec->fpath << std::endl;
    return false;
  }

#if 0	// TODO no such field in pre-gridded files
  NcVar *netCoherentPower = file.get_var(ncpName.c_str());
  if (netCoherentPower == NULL) {
    std::cerr << "Can't get Net Coherent Power '" << ncpName
	      << "' from " << fspec->fpath << std::endl;
    return false;
  }
#endif
  
  float ref_fill;
  float vel_fill;
  // float ncp_fill;

  if (! getFillValue(reflectivity, ref_fill) ) {
    std::cerr << "Can't get reflectivity fill value from file" << std::endl;
    return false;
  }
  
  if (! getFillValue(velocity, vel_fill) ) {
    std::cerr << "Can't get velocity fill value from file" << std::endl;
    return false;
  }

#if 0
  if (! getFillValue(netCoherentPower, ncp_fill) ) {
    std::cerr << "Can't get NCP fill value from file"  << std::endl;
    return false;
  }
#endif
  
  // iDim, jDim, and kDim are float. That doesn't work very well for pointer arithmetic
  int xDim = (int) iDim;
  int yDim = (int) jDim;
  int zDim = (int) kDim;

  // Allocate buffers to hold the data

  float *ref = new float[xDim * yDim];
  float *vel = new float[xDim * yDim];
  // float *ncp;
    
  // Grab the data

  int timeIndex = 0;

  for(int k = 0; k < zDim; k++) {
    if (! reflectivity->set_cur(timeIndex, k, 0, 0, -1) ) {
	std::cerr << "Couldn't set reflectivity corner" << std::endl;
	break;
      }
    if (! velocity->set_cur(timeIndex, k, 0, 0, -1) ) {
	std::cerr << "Couldn't set velocity corner" << std::endl;
	break;
      }
#if 0
    if (! netCoherentPower->set_cur(timeIndex, k, 0, 0, -1) ) {
	std::cerr << "Couldn't set NCP width corner" << std::endl;
	break;
    }
#endif
    if (! reflectivity->get(ref, 1, 1, yDim, xDim) ) {
      std::cerr << "Couldn't get reflectivity value" << std::endl;
      break;
    }
    if (! velocity->get(vel, 1, 1, yDim, xDim) ) {
      std::cerr << "Couldn't get Velocity value" << std::endl;
      break;
    }
#if 0
    if (! netCoherentPower->get(ncp, 1, 1, yDim, xDim) ) {
      std::cerr << "Couldn't get NCP value" << std::endl;
      break;
    }
#endif
    
    //	std::cerr << "[" << i << ", " << j << ", " << k << "] " <<
    //	  "ref: " << ref << ", vel: " << vel << std::endl;
    float v;

    // Looks like x and y are swapped in the NetCDF file.

    long ipt = -1;

    int coordz = zmin + k * kGridSp;
    
    for(int i = 0; i < xDim; i++) {
      for(int j = 0; j < yDim; j++) {
	
	ipt++;	// number of points

#if 0
	int coordx = xmin + i * iGridSp;
	int coordy = ymin + j * jGridSp;
#endif
	int coordy = ymin + i * iGridSp;
	int coordx = xmin + j * jGridSp;
	
	pointBbox->addOb(coordz, coordy, coordx);

	v = *(ref + i * yDim + j);		// reflectivity (REF)
	if (v <= ref_fill)
	  v = -999;
	double valDbz = v;
	
	v = *(vel + i * yDim + j);		// dopler velocity magnitude (VU)
	if (v <= vel_fill)
	  v = -999;
	double valVg = v;
#if 0	
	v = *(ncp + i * yDim + j);		// spectral grid width (SW)
	if (v <= ncp_fill)
	  v = -999;
	double valNcp = v;
#endif
	double valNcp = 1.0;
	
	float aircraftAltKmMsl = 0;             // TODO radar alt?

	float thetaRad = atan2(coordy, coordx);  // azimuth?
	float elevRad = atan2(coordz, sqrt( coordx * coordx + coordy * coordy)); // elevation?
	
	pointVec->push_back( new Point(
            ifile,
	    0,	// assume only 1 ray
            ipt,
	    
            aircraftAltKmMsl,	// TODO radar alt?
            0, // aircrafty,		// Radar x and y?
            0, // aircraftx,
	    
	    coordz,
            coordy,
            coordx,
	    
            time,               // UTC seconds since 1970 jan 1.
            thetaRad,           // polar coords angle
            elevRad,            // angle from the aircraft between horizontal and pt
	    
            valVg,              // radial velocity
            valDbz,             // log reflectivity
            valNcp));           // radar net coherent power

	statVg->addOb(valVg);
	statDbz->addOb(valDbz);
	statNcp->addOb(valNcp);
      }
    }
  }
  
  delete[] ref;
  delete[] vel;

#if 0
  free(ncp);
#endif
  return true;
}

bool getFillValue(NcVar *var, float &val)
{
  NcAtt *fv = var->get_att("_FillValue");
  if (fv == NULL)
    return false;
  val = fv->as_float(0);
  return true;
}

bool getTime(NcFile &file, const char *varName, double &time){
  NcVar *var = file.get_var(varName);
  if (var == NULL) return false;

  return var->get(&time, 1);
}
  
bool getDimInfo(NcFile &file, int dim, const char *varName,
			   float &spacing, float &min, float &max)
{
  NcVar *var = file.get_var(varName);
  if (var == NULL) return false;
    
  float *vals = new float[dim];
  bool retVal = var->get(vals, dim);
  
  if(retVal) {
    spacing = vals[1] - vals[0];
    min = vals[0];
    max = vals[dim - 1];
  }
  delete[] vals;
  return retVal;
}

// Radar location. usually at (0, 0) kilometers

bool getGridMapping(NcFile &file, float &radar_lat, float &radar_lon)
{
  NcVar *grid_mapping = file.get_var("grid_mapping_0");
  if (grid_mapping == NULL)
    return false;
  NcAtt *olat = grid_mapping->get_att("latitude_of_projection_origin");
  if (olat == NULL)
    return false;
  radar_lat = olat->as_float(0);
  
  NcAtt *olon = grid_mapping->get_att("longitude_of_projection_origin");
  if (olon == NULL)
    return false;
  radar_lon = olon->as_float(0);
  return true;
}

// Get Lat and Lon at lower left corner of grid.

bool getOriginLatLon(NcFile &file, float &lat, float &lon)
{
  NcVar *lat0 = file.get_var("lat0");
  if ( lat0 == NULL )
    return false;
  
  NcVar *lon0 = file.get_var("lon0");
  if ( lon0 == NULL )
    return false;
  
  if (! lat0->set_cur(0, 0) ) {
    std::cerr << "Couldn't set lat corner" << std::endl;
    return false;
  }
  if (! lon0->set_cur(0, 0) ) {
    std::cerr << "Couldn't set lat corner" << std::endl;
    return false;
  }

  if (! lat0->get(&lat, 1, 1) ) {
    std::cerr << "Couldn't get lat at (0, 0)" << std::endl;
    return false;
  }
  
  if (! lon0->get(&lon, 1, 1) ) {
    std::cerr << "Couldn't get lon at (0, 0)" << std::endl;
    return false;
  }
#if 0
  std::cout << "** lat(0, 0): " << lat << ", lon(0, 0): " << lon << std::endl;
#endif
  
  return true;
}

bool getGlobalAttribute(NcFile &file, std::string name, float &val) {
  NcAtt *attr = file.get_att(name.c_str());
  if ( attr == NULL )
    return false;
  val = attr->as_float(0);
  return true;
}
