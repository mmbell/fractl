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

#include <sys/stat.h>
#include <cerrno>
#include <cstring>

#include "Fractl.hh"
#include "Params.hh"

// NetCDF: See doc at:
//   netcdfC/tda/netcdf-4.1.3/examples/CXX/pres_temp_4D_wr.cpp
//   /usr/local/netcdf4/include/netcdfcpp.h

#include <Ncxx/Nc3File.hh>

// Write the variables in the cellMat to the Netcdf file outNc.

void Fractl::writeNetcdf()
{

  long num_x, num_y, num_z;

  if (gridType == Params::GRID_MISH) { // TODO Verify the + 4
    num_x = (xgridmax - xgridmin) * 2 / xgridinc + 4;
    num_y = (ygridmax - ygridmin) * 2 / ygridinc + 4;  
    num_z = (zgridmax - zgridmin) * 2 / zgridinc + 4;      
  }
  else {
    num_x = nradx;
    num_y = nrady;
    num_z = nradz;
  }
  // Create linear arrays
  //  long totalLen = nradz * nrady * nradx;
  long totalLen = num_z * num_y * num_x;
  float * wLinear = new float[ totalLen];
  float * vLinear = new float[ totalLen];
  float * uLinear = new float[ totalLen];
  float * dbzLinear = new float[ totalLen];
  float * ncpLinear = new float[ totalLen];
  float * condNumLinear = new float[ totalLen];
  float * ustdLinear = new float[ totalLen];
  float * vstdLinear = new float[ totalLen];
  float * wstdLinear = new float[ totalLen];
  for (long iz = 0; iz < num_z; iz++) {
    for (long iy = 0; iy < num_y; iy++) {
      for (long ix = 0; ix < num_x; ix++) {
	
        // long kk = iz * nrady * nradx + iy * nradx + ix;    // linear index
	long kk = iz * num_y * num_x + iy * num_x + ix;    // linear index
        wLinear[kk] = cellMat[iz][iy][ix].ww;
        vLinear[kk] = cellMat[iz][iy][ix].vv;
        uLinear[kk] = cellMat[iz][iy][ix].uu;
        dbzLinear[kk] = cellMat[iz][iy][ix].meanNbrDbz;
        ncpLinear[kk] = cellMat[iz][iy][ix].meanNbrNcp;
        condNumLinear[kk] = cellMat[iz][iy][ix].conditionNumber;
        ustdLinear[kk] = cellMat[iz][iy][ix].ustd;
        vstdLinear[kk] = cellMat[iz][iy][ix].vstd;
        wstdLinear[kk] = cellMat[iz][iy][ix].wstd;

      } // for ix
    } // for iy
  } // for iz

  string fname;
  if (outNc[outNc.length()-1] == '/') {
    // file name format: outNc/yyyymmdd/ncf_yyyymmdd_hhmmss.nc
    time_t tm = timeMax;
    struct tm tmstr;
    gmtime_r( &tm, &tmstr);

    char tbuf[1000];
    strftime( tbuf, 1000, "%Y%m%d%H%M%S", &tmstr);
    string fulldate(tbuf, 0, 14);
    string subdirname = outNc + fulldate.substr( 0, 8);
    if (bugs >= Params::DEBUG_NORM) {
      cout << "fulldate: \"" << fulldate << "\"" << endl;
      cout << "subdirname: \"" << subdirname << "\"" << endl;
    }

    if ( (mkdir( subdirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) &&
	 (errno != EEXIST))  {
      std::cerr <<  "Couldn't mkdir" <<  subdirname << ": "
		<< std::strerror(errno) << std::endl
		<< "                 Writing the output to /tmp" << std::endl;
      subdirname = "/tmp";
    }

    fname = subdirname + "/ncf_" + fulldate.substr( 0, 8) + "_"
      + fulldate.substr( 8, 6) + ".nc";
  }
  else {
    fname = outNc;
  }
  if (bugs >= Params::DEBUG_NORM) cout << "fname: \"" << fname << "\"" << endl;
  Nc3File ncout( fname.c_str(), Nc3File::Replace);
  if ( ! ncout.is_valid() ) {
    std::cerr << "Couldn't create output NetCDF file " << fname << std::endl;
    return;
  }

  ncout.add_att("Conventions", "CF-1.5");

  // Dimension arrays
  int ntime = 1;
  double * tvals = new double[ ntime];
  tvals[0] = timeMax;

  
  double * zvals = new double[num_z];
  double * yvals = new double[num_y];
  double * xvals = new double[num_x];

  if (gridType == Params::GRID_MESH) {
    for (long ii = 0; ii < num_z; ii++) {
      zvals[ii] = 1000 * (zgridmin + ii * zgridinc);  // convert km to m
    }

    for (long ii = 0; ii < nrady; ii++) {
      yvals[ii] = 1000 * (ygridmin + ii * ygridinc);
    }

    for (long ii = 0; ii < nradx; ii++) {
      xvals[ii] = 1000 * (xgridmin + ii * xgridinc);
    }
  } else {
    fillMishAxisCoordinates(zvals, nradz, zgridmin, zgridinc);
    fillMishAxisCoordinates(yvals, nrady, ygridmin, ygridinc);
    fillMishAxisCoordinates(xvals, nradx, xgridmin, xgridinc);    
  }
  
  // Dimensions
  Nc3Dim * tdim = ncout.add_dim("time", ntime);
  Nc3Dim * zdim = ncout.add_dim("z0", num_z);
  Nc3Dim * ydim = ncout.add_dim("y0", num_y);
  Nc3Dim * xdim = ncout.add_dim("x0", num_x);

  // Only if we output a Samurai mish
#if 0
  Nc3Dim *sam_idim, *sam_jdim, *sam_kdim;
  
  Nc3Dim *sam_imin, *sam_imax, *sam_iincr;
  Nc3Dim *sam_jmin, *sam_jmax, *sam_jincr;
  Nc3Dim *sam_kmin, *sam_kmax, *sam_kincr; 
#endif
  
  // If writing out a Samurai Mish, add variables so that Samurai can reconstruct the grid
  
  if (gridType == Params::GRID_MISH) {
    ncout.add_att("sam_idim",  nradx);
    ncout.add_att("sam_jdim", nrady);
    ncout.add_att("sam_kdim", nradz);

    ncout.add_att("sam_imin", xgridmin);
    ncout.add_att("sam_imax", xgridmax);
    ncout.add_att("sam_iincr", xgridinc);    

    ncout.add_att("sam_jmin", ygridmin);
    ncout.add_att("sam_jmax", ygridmax);
    ncout.add_att("sam_jincr", ygridinc);    

    ncout.add_att("sam_kmin", zgridmin);
    ncout.add_att("sam_kmax", zgridmax);
    ncout.add_att("sam_kincr", zgridinc);    
  }
  
  const Nc3Dim * dataDims4[] = { tdim, zdim, ydim, xdim};
  const Nc3Dim * dataDims2[] = { ydim, xdim};
  const Nc3Dim * dataDims1[] = { tdim};

  // Dimension vars
  Nc3Var * tvar = ncout.add_var("time", nc3Double, tdim);
  tvar->add_att("standard_name", "time");
  tvar->add_att("units", "seconds since 1970-01-01T00:00:00Z");
  tvar->add_att("calendar", "gregorian");

  Nc3Var * zvar = ncout.add_var("z0", nc3Double, zdim);
  zvar->add_att("standard_name", "height");
  zvar->add_att("units", "m");
  zvar->add_att("positive", "up");

  Nc3Var * yvar = ncout.add_var("y0", nc3Double, ydim);
  yvar->add_att("standard_name", "projection_y_coordinate");
  yvar->add_att("units", "m");

  Nc3Var * xvar = ncout.add_var("x0", nc3Double, xdim);
  xvar->add_att("standard_name", "projection_x_coordinate");
  xvar->add_att("units", "m");

  // Metadata vars
  Nc3Var * gridMapVar = ncout.add_var("grid_mapping_0", nc3Int);
  gridMapVar->add_att("grid_mapping_name", "transverse_mercator");
  gridMapVar->add_att("semi_major_axis", earthRadiusMeter);
  gridMapVar->add_att("semi_minor_axis", earthRadiusMeter);
  gridMapVar->add_att("inverse_flattening", 1.0 / flattening);
  gridMapVar->add_att("latitude_of_projection_origin", projLat0);
  gridMapVar->add_att("longitude_of_projection_origin", projLon0);
  gridMapVar->add_att("false_easting", 0);
  gridMapVar->add_att("false_northing", 0);

  Nc3Var * latVar = ncout.add_var("lat0", nc3Double, 2, dataDims2);
  latVar->add_att("standard_name", "latitude");
  latVar->add_att("units", "degrees_north");

  Nc3Var * lonVar = ncout.add_var("lon0", nc3Double, 2, dataDims2);
  lonVar->add_att("standard_name", "longitude");
  lonVar->add_att("units", "degrees_east");

  // Data variables

  Nc3Var * timeStartVar = ncout.add_var("start_time", nc3Double, 1, dataDims1);
  timeStartVar->add_att("standard_name", "start_time");
  timeStartVar->add_att("units", "s");

  Nc3Var * timeStopVar = ncout.add_var("stop_time", nc3Double, 1, dataDims1);
  timeStopVar->add_att("standard_name", "stop_time");
  timeStopVar->add_att("units", "s");

  Nc3Var * wVar = ncout.add_var("W", nc3Float, 4, dataDims4);
  wVar->add_att("standard_name", "upward_air_velocity");
  wVar->add_att("units", "m s-1");
  wVar->add_att("grid_mapping", "grid_mapping_0");

  Nc3Var * vVar = ncout.add_var("V", nc3Float, 4, dataDims4);
  vVar->add_att("standard_name", "northward_wind");
  vVar->add_att("units", "m s-1");
  vVar->add_att("grid_mapping", "grid_mapping_0");

  Nc3Var * uVar = ncout.add_var("U", nc3Float, 4, dataDims4);
  uVar->add_att("standard_name", "eastward_wind");
  uVar->add_att("units", "m s-1");
  uVar->add_att("grid_mapping", "grid_mapping_0");

  Nc3Var * dbzVar = ncout.add_var("DBZ", nc3Float, 4, dataDims4);
  dbzVar->add_att("standard_name", "mean_neighbor_dbz");
  dbzVar->add_att("units", "1");  // Apparently CF uses "1" to mean ratio
  dbzVar->add_att("grid_mapping", "grid_mapping_0");

  Nc3Var * ncpVar = ncout.add_var("NCP", nc3Float, 4, dataDims4);
  ncpVar->add_att("standard_name", "mean_neighbor_ncp");
  ncpVar->add_att("units", "1");  // Apparently CF uses "1" to mean ratio
  ncpVar->add_att("grid_mapping", "grid_mapping_0");

  Nc3Var * condNumVar = ncout.add_var("conditionNumber", nc3Float, 4, dataDims4);
  condNumVar->add_att("standard_name", "matrix_condition_number");
  condNumVar->add_att("units", "1");
  condNumVar->add_att("grid_mapping", "grid_mapping_0");

  Nc3Var * wstdVar = ncout.add_var("W_std", nc3Float, 4, dataDims4);
  wstdVar->add_att("standard_name", "w_std_deviation");
  wstdVar->add_att("units", "m s-1");
  wstdVar->add_att("grid_mapping", "grid_mapping_0");

  Nc3Var * vstdVar = ncout.add_var("V_std", nc3Float, 4, dataDims4);
  vstdVar->add_att("standard_name", "v_std_deviation");
  vstdVar->add_att("units", "m s-1");
  vstdVar->add_att("grid_mapping", "grid_mapping_0");

  Nc3Var * ustdVar = ncout.add_var("U_std", nc3Float, 4, dataDims4);
  ustdVar->add_att("standard_name", "u_std_deviation");
  ustdVar->add_att("units", "m s-1");
  ustdVar->add_att("grid_mapping", "grid_mapping_0");

  // Write coord vars
  tvar->put( tvals, ntime);
  zvar->put( zvals, num_z);
  yvar->put( yvals, num_y);
  xvar->put( xvals, num_x);

  // Write metadata vars
  int gridMapValue = 0;
  gridMapVar->put( &gridMapValue, 1);

  double * latLinear = new double[num_y * num_x];
  double * lonLinear = new double[num_y * num_x];

  if (gridType == Params::GRID_MESH) {  
    for (long iy = 0; iy < num_y; iy++) {
      for (long ix = 0; ix < num_x; ix++) {
	double latDeg;
	double lonDeg;
	// TODO These are not right for the MISH
	double coordy = ygridmin + iy * ygridinc;
	double coordx = xgridmin + ix * xgridinc;
	yxToLatLon( // tranMerc,
		   projLon0,
		   // basey,
		   // basex,
		   coordy,
		   coordx,
		   latDeg,             // output value
		   lonDeg);            // output value
	long kk = iy * num_x + ix;    // linear index
	latLinear[kk] = latDeg;
	lonLinear[kk] = lonDeg;
      } // for ix
    } // for iy
  } else {
    fillMishLatLon(latLinear, lonLinear);
  }
  
  long counts2[] = { num_y, num_x};
  latVar->put( latLinear, counts2);
  lonVar->put( lonLinear, counts2);

  delete[] latLinear;
  delete[] lonLinear;

  // Write data vars
  long counts1[] = { ntime};
  double * timeStarts = new double[ntime];
  double * timeStops = new double[ntime];
  for (int ii = 0; ii < ntime; ii++) {
    timeStarts[ii] = timeMin;
    timeStops[ii] = timeMax;
  }
  timeStartVar->put( timeStarts, counts1);
  timeStopVar->put( timeStops, counts1);
  delete[] timeStarts;
  delete[] timeStops;


  long counts4[] = { ntime, num_z, num_y, num_x};
  wVar->put( wLinear, counts4);
  vVar->put( vLinear, counts4);
  uVar->put( uLinear, counts4);
  dbzVar->put( dbzLinear, counts4);
  ncpVar->put( ncpLinear, counts4);
  condNumVar->put( condNumLinear, counts4);
  wstdVar->put( wstdLinear, counts4);
  vstdVar->put( vstdLinear, counts4);
  ustdVar->put( ustdLinear, counts4);

  // The Nc3File destructor automatically closes the file.

  delete[] tvals;
  delete[] zvals;
  delete[] yvals;
  delete[] xvals;

  delete[] wLinear;
  delete[] vLinear;
  delete[] uLinear;
  delete[] dbzLinear;
  delete[] ncpLinear;
  delete[] condNumLinear;
  delete[] wstdLinear;
  delete[] vstdLinear;
  delete[] ustdLinear;

} // end writeNetcdf

// Fill a Samurai mish coordinate axis

void Fractl::fillMishAxisCoordinates(double *vals, long nrad, double gridmin, double gridincr)
{
  for (int i = -1; i < (nrad); i++) {
    for (int mu = -1; mu <= 1; mu += 2) {
      double pos = gridmin + gridincr * (i + (0.5 * sqrt(1.0 / 3.0) * mu + 0.5));
      long index = (i + 1) * 2 + (mu + 1) / 2;
      vals[index] = 1000 * pos; // convert km to m
    }
  }
}

void Fractl::fillMishLatLon(double *latLinear, double *lonLinear)
{
  // for (long iy = 0; iy < num_y; iy++)
  // for (long ix = 0; ix < num_x; ix++)

  // TODO verify the + 4
  
  long num_x = (xgridmax - xgridmin) * 2 / xgridinc + 4;
  
  for (int ji = -1; ji < (nrady); ji++) {
    for (int jmu = -1; jmu <= 1; jmu += 2) {
      double coordy = ygridmin + ygridinc * (ji + (0.5 * sqrt(1. / 3.) * jmu + 0.5));

      for (int ii = -1; ii < (nradx); ii++) {
	for (int imu = -1; imu <= 1; imu += 2) {
	  double coordx = xgridmin + xgridinc * (ii + (0.5 * sqrt(1. / 3.) * imu + 0.5));

	  long ix = (ii + 1) * 2 + (imu + 1) / 2;
	  long iy = (ji + 1) * 2 + (jmu + 1) / 2;
	  
	  double latDeg;
	  double lonDeg;
      
	  yxToLatLon( // tranMerc,
		     projLon0,
		     coordy,
		     coordx,
		     latDeg,             // output value
		     lonDeg);            // output value
	  long kk = iy * num_x + ix;    // linear index
	  latLinear[kk] = latDeg;
	  lonLinear[kk] = lonDeg;
	}
      }
    }
  }
}
