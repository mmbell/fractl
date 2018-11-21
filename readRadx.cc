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

#include <iostream>
#include "Fractl.hh"
#include "Params.hh"

// See Radx doc at:
// http://www.ral.ucar.edu/projects/titan/docs/radial_formats/

#include <Radx/Radx.hh>
#include <Radx/RadxField.hh>
#include <Radx/RadxGeoref.hh>
#include <Radx/RadxRangeGeom.hh>
#include <Radx/RadxRay.hh>
#include <Radx/RadxVol.hh>
#include <Radx/NcfRadxFile.hh>
#include <Radx/DoradeRadxFile.hh>
#include <Radx/RadxTime.hh>
#include <Radx/RadxTimeList.hh>

// Read a radar file and append the Points to pointVec.

// Refraction calculations
// See the accompanying inkscape diagram in file refraction.svg.
//
// Let:
//   re = radius of earth, assumed spherical
//   rv = radius of curvature of the radar beam
//   elev = beam elevation angle
//         = angle of the radar beam with respect to the horizontal
//           plane through the aircraft.
//   k = constant of rv / re when elev == 0.  Typically k = 4/3.
//   ha = height of aircraft above MSL
//
// rv = k re / cos(elev)                      eqn 1
//
// Consider a 2 dim projection of the earth and aircraft with
// the earth center at (0,0), the ground beneath the aircraft at (0,re),
// and the aircraft at (0, re+hp).
// Assume the radar beam lies in the 2 dim projection
// and has angle elev with respect to the x axis.
//
// Let
//   cx = x coord of the center of the arc created by the radar beam.
//     This arc has radius rv.
//   cy = y coord of the center of the arc created by the radar beam.
//
// cx = rv * sin(elev)
//
// re + ha - cy = rv * cos(elev) = k * re          // from eqn 1
// cy = re + hp - k * re = re * (1 - k) + hp        // generally cy < 0
//
// Let
//   px,py = coords of the reflecting hydrometeoroid.
//   sp = distance along the radar beam from the aircraft to
//     the reflecting hydrometeoroid.
//   phi = angle portended at cx,cy by the arc of the radar beam
//     from the aircraft to px,py.
//   gamma = angle, at cx,cy, from the parallel to the x axis,
//     to px,py.
//
// phi = sp / rv
// gamma = pi/2 + elev - phi
//
// px = cx + rv * cos(gamma)
// py = cy + rv * sin(gamma)
//
// Let
//   hp = height of the hydrometeoroid above MSL
//   se = distance along the earth surface from the point directly under
//     the aircraft to the point directly under the reflecting hydrometeoroid.
//
// (re + hp)^2 = px^2 + py^2
// hp = sqrt( px^2 + py^2) - re
//
//
// Rough checks:
// angle zeta = arctan( py/px)
// se/re = pi/2 - zeta
// se = re * (pi/2 - zeta)
//
// For sp << re, the path is essentially a straight line
// at angle elev, and we should have:
//   se =about sp * cos(elev)
//   hp =about ha + sp * sin(elev)

// TODO split this (500+ lines)

void Fractl::readRadarFile(
  long forceOk,                     // if != 0, force all point ncp and dbz to ok
  long ifile,
  FileSpec * fspec,
  double maxAbsElevDeg,
  Bbox * aircraftBbox,              // overall bounding box of aircraft locs
  Bbox * pointBbox,                 // overall bounding box of point locs
  Statistic * statVg,               // statistics for radial velocity
  Statistic * statDbz,              // statistics for DBZ
  Statistic * statNcp,              // statistics for net coherent power
  vector<Point *> *pointVec,        // appended
  double *maxAbsErrHoriz,
  double *maxAbsErrVert,
  double &timeMin,                  // returned
  double &timeMax)                  // returned
{
  RadxFile rxfile;
  RadxVol rxvol;

  if (0 != rxfile.readFromPath( fspec->fpath, rxvol))
    badparms("cannot read file: \"%s\"", fspec->fpath.c_str());

  double tmStart = rxvol.getStartTimeSecs() + 1.e-9 * rxvol.getStartNanoSecs();
  double tmEnd = rxvol.getEndTimeSecs() + 1.e-9 * rxvol.getEndNanoSecs();

  std::cout << "\n\nreadRadarFile:" << std::endl
    << "    fpath:         " << fspec->fpath << std::endl
    << "    startTime:     " << formatTime( tmStart) << std::endl
    << "    endTime:       " << formatTime( tmEnd) << std::endl
    << "    latitudeDeg:   " << rxvol.getLatitudeDeg() << std::endl
    << "    longitudeDeg:  " << rxvol.getLongitudeDeg() << std::endl
    << "    altitudeKm:    " << rxvol.getAltitudeKm() << std::endl;

  if (bugs >= Params::DEBUG_VERBOSE) {
    std::cout << std::endl << std::endl << "readRadarFile: RadxVol: " << std::endl;
    rxvol.print( std::cout);
    std::cout << "end RadxVol" << std::endl;
  }

  // TODO: This should be set with method calls
  fspec->altitudeKmMsl = rxvol.getAltitudeKm();
  fspec->latitudeDeg = rxvol.getLatitudeDeg();
  fspec->longitudeDeg = rxvol.getLongitudeDeg();

  // Get coords for debug only
  double tmpy;            // northing in meters
  double tmpx;            // easting in meters

  latLonToYX(
    projLon0,             // central meridian of projection
    basey,                // coord y base, km
    basex,                // coord x base, km
    fspec->latitudeDeg,
    fspec->longitudeDeg,
    tmpy,                 // output value
    tmpx);                // output value

  if (bugs >= Params::DEBUG_VERBOSE) {
    std::cout << std::endl << "readRadarFile: fspec:" << std::endl;
    fspec->print();
    std::cout << "readRadarFile: x (km): " << tmpx << std::endl;
    std::cout << "readRadarFile: y (km): " << tmpy << std::endl;
    std::cout << "set arrow " << (9000000 + ifile)
	      << " from " << tmpx << "," << tmpy
	      << " to " << tmpx << "," << tmpy
	      << std::endl;
  }
  vector<RadxRay *> rays = rxvol.getRays();
  long itotpt = 0;
  long numGoodPoint = 0;
  long numInvalidPoint = 0;
  long numMissLocPoint = 0;

  struct timeval timeh;
  addDeltaTime( &timeh, NULL);

  for (size_t iray = 0; iray < rays.size(); iray++) {
    RadxRay * ray = rays[iray];
    if (ray == NULL)
      badparms("cannot read file: \"%s\"", fspec->fpath.c_str());

    // RadxRay extends RadxRangeGeom
    double aircraftLatDeg = numeric_limits<double>::quiet_NaN();
    double aircraftLonDeg = numeric_limits<double>::quiet_NaN();
    double aircraftAltKmMsl = numeric_limits<double>::quiet_NaN();

    // gref is NULL for some ground-based stations.
    const RadxGeoref * gref = ray->getGeoreference();
    if (gref == NULL) {
      aircraftAltKmMsl = fspec->altitudeKmMsl;
      aircraftLatDeg = fspec->latitudeDeg;
      aircraftLonDeg = fspec->longitudeDeg;
    }
    else {
      aircraftLatDeg = gref->getLatitude();
      aircraftLonDeg = gref->getLongitude();
      aircraftAltKmMsl = gref->getAltitudeKmMsl();
    }
    if (std::isnan( aircraftLatDeg)
      || std::isnan( aircraftLonDeg)
      || std::isnan( aircraftAltKmMsl))
    {
      std::cout << "Error: incomplete file spec or radx Georeference" << std::endl;
      fspec->print();
      std::cout.flush();
      throwerr("incomplete file spec or radx Georeference");
    }

    double aircrafty;       // northing in meters
    double aircraftx;       // easting in meters

    latLonToYX(
      projLon0,             // central meridian of projection
      basey,                // coord y base, km
      basex,                // coord x base, km
      aircraftLatDeg,
      aircraftLonDeg,
      aircrafty,            // output value
      aircraftx);           // output value

    aircraftBbox->addOb( aircraftAltKmMsl, aircrafty, aircraftx);

    double rayTime = ray->getTimeDouble();
    if (std::isnan(timeMin) || rayTime < timeMin) timeMin = rayTime;
    if (std::isnan(timeMax) || rayTime > timeMax) timeMax = rayTime;

    // Print info for first ray only
    if (bugs >= Params::DEBUG_NORM && iray == 0) {
      std::cout << setprecision(7);
      std::cout << "first_ray_loc:"
        << "  fpath: " << fspec->fpath
        << "  time: " << setprecision(16) << ray->getTimeDouble()
        << "  " << formatTime( ray->getTimeDouble()) << std::endl
        << "  aircraftAltKmMsl: " << aircraftAltKmMsl
        << "  aircraftLatDeg: " << aircraftLatDeg
        << "  aircraftLonDeg: " << aircraftLonDeg << std::endl
        << "  aircrafty: " << aircrafty
        << "  aircraftx: " << aircraftx << std::endl;
      std::cout << "  RadxRangeGeom:" << std::endl;
      std::cout << "  getStartRangeKm(): " << ray->getStartRangeKm() << std::endl;
      std::cout << "  getGateSpacingKm(): " << ray->getGateSpacingKm() << std::endl;

      std::cout << "  RadxRay:" << std::endl;
      std::cout << "  getAzimuthDeg(): " << ray->getAzimuthDeg() << std::endl;
      std::cout << "  getElevationDeg(): " << ray->getElevationDeg() << std::endl;

      if (gref != NULL) {
        std::cout << "  RadxGeoref:" << std::endl;
        std::cout << "  getRotation(): " << gref->getRotation() << std::endl;
        std::cout << "  getEwVelocity(): " << gref->getEwVelocity() << std::endl;
        std::cout << "  getNsVelocity(): " << gref->getNsVelocity() << std::endl;
        std::cout << "  getVertVelocity(): " << gref->getVertVelocity() << std::endl;
        std::cout << "  getHeading(): " << gref->getHeading() << std::endl;
        std::cout << "  getRoll(): " << gref->getRoll() << std::endl;
        std::cout << "  getPitch(): " << gref->getPitch() << std::endl;
        std::cout << "  getDrift(): " << gref->getDrift() << std::endl;
        std::cout << "  getRotation(): " << gref->getRotation() << std::endl;
        std::cout << "  getTilt(): " << gref->getTilt() << std::endl;
      }

      // http://www.ral.ucar.edu/projects/titan/docs/radial_formats/
      // Standard names for moments variables:
      //   radial_velocity_of_scatterers_away_from_instrument VEL m/s

      // Eldora: "VG".  CSU-Chill: "VE". Michael's synthetic data: VR
      // RadxField * fieldVg = ray->getField( *radialName);
      RadxField * fieldVg = ray->getField(radialName);
      if (fieldVg == NULL) throwerr("radialName not found");
      std::cout << std::endl << std::endl
        << "readRadarFile: first ray fieldVg print: " << std::endl;
      fieldVg->print( std::cout);
      std::cout << std::endl;

      // Eldora: "DBZ".  CSU-Chill: "DZ".
      RadxField * fieldDbz = ray->getField( dbzName);
      if (fieldDbz == NULL) throwerr("dbzName not found");
      std::cout << std::endl << std::endl
        << "readRadarFile: first ray fieldDbz print: " << std::endl;
      fieldDbz->print( std::cout);
      std::cout << std::endl;

      // Eldora: "NCP".  CSU-Chill: "NC".
      RadxField * fieldNcp = ray->getField( ncpName);
      if (fieldNcp == NULL) throwerr("ncpName not found");
      std::cout << std::endl << std::endl
        << "readRadarFile: first ray fieldNcp print: " << std::endl;
      fieldNcp->print( std::cout);
      std::cout << std::endl;
    } // if bugs >= 0 && iray == 0

    if (bugs >= Params::DEBUG_VERBOSE) {
      std::cout << setprecision(15);
      std::cout << std::endl << std::endl << "readRadarFile: RadxRay: " << std::endl;
      ray->print( std::cout);
      std::cout << std::endl;
      std::cout << "end RadxRay" << std::endl;
      std::cout << "readRadarFile: iray: " << iray
        //<< "  rotation: " << ray->getGeoreference()->getRotation()
        << "  ray time: " << formatTime( ray->getTimeDouble()) << std::endl;
    }

    RadxField * fieldVg = ray->getField( radialName);
    RadxField * fieldDbz = ray->getField( dbzName);
    RadxField * fieldNcp = ray->getField( ncpName);

    if (fieldVg == NULL) throwerr("radialName not found");
    if (fieldDbz == NULL) throwerr("dbzName not found");
    if (fieldNcp == NULL) throwerr("ncpName not found");

    if (bugs >= Params::DEBUG_NORM && iray == 0) {
      std::cout << setprecision(15);
      std::cout << "readRadarFile: fieldVg getNRays:    "
        << fieldVg->getNRays() << std::endl;
      std::cout << "readRadarFile: fieldDbz getNRays:   "
        << fieldDbz->getNRays() << std::endl;
      std::cout << "readRadarFile: fieldNcp getNRays:   "
        << fieldNcp->getNRays() << std::endl;

      std::cout << "readRadarFile: fieldVg getNPoints:  "
        << fieldVg->getNPoints() << std::endl;
      std::cout << "readRadarFile: fieldDbz getNPoints: "
        << fieldDbz->getNPoints() << std::endl;
      std::cout << "readRadarFile: fieldNcp getNPoints: "
        << fieldNcp->getNPoints() << std::endl;
    }

    double aziDeg = ray->getAzimuthDeg();
    double elevDeg = ray->getElevationDeg();

    double startDistKm = ray->getStartRangeKm();
    double spacingKm = ray->getGateSpacingKm();

    // Missing values in Radx:
    //
    // getDataType()      returns 1 == SI16
    // getDoubleValue(i)  returns -9999.00000
    // getStoredValue(i)  returns -9999.00000
    //
    // getMissing()       returns -32768
    // getMissingFl64()   returns -9999
    // getMissingSi16()   returns -32768
    //
    // If a program retrieves via getDoubleValue(),
    // it must compare with getMissingFl64(), not getMissing().
    //
    // For missing values getStoredValue has the same behavior
    // as getDoubleValue - both return _missingFl64.

    double missVg = fieldVg->getMissingFl64();
    double missDbz = fieldDbz->getMissingFl64();
    double missNcp = fieldNcp->getMissingFl64();

    if (bugs >= Params::DEBUG_NORM && iray == 0) {
      Radx::DataType_t typeVg = fieldVg->getDataType();
      Radx::DataType_t typeDbz = fieldDbz->getDataType();
      Radx::DataType_t typeNcp = fieldNcp->getDataType();

      std::cout << setprecision(15);
      std::cout << "  missVg: " << missVg << std::endl;
      std::cout << "  missDbz: " << missDbz << std::endl;
      std::cout << "  missNcp: " << missNcp << std::endl;
      std::cout << "  typeVg: " << typeVg << std::endl;
      std::cout << "  typeDbz: " << typeDbz << std::endl;
      std::cout << "  typeNcp: " << typeNcp << std::endl;
      std::cout << "  vg getMissingFl64(): " << fieldVg->getMissingFl64() << std::endl;
      std::cout << "  vg getMissingSi16(): " << fieldVg->getMissingSi16() << std::endl;
    }

    struct timeval timei;
    addDeltaTime( &timei, NULL);

    // TODO bring out loop invariants

    for (size_t ipt = 0; ipt < fieldVg->getNPoints(); ipt++) {
      struct timeval timej;
      addDeltaTime( &timej, NULL);

      // Calcs that assume a locally flat earth and no refraction
      double elevRad = elevDeg * M_PI / 180;
      double slantDistKm = startDistKm + ipt * spacingKm;
      double flatHorizDistKm = cos( elevRad) * slantDistKm;
      double flatVertDistKm = sin( elevRad) * slantDistKm;
      double flatAltKmMsl = aircraftAltKmMsl + flatVertDistKm;

      // Calcs that assume a spherical earth with refraction.
      // See notes at the top.
      double krefract = 4.0 / 3.0;
      double earthRadiusKm = 0.001 * geodesic->MajorRadius();
      double refractRadiusKm = krefract * earthRadiusKm / cos( elevRad);
      double cx = refractRadiusKm * sin( elevRad);
      double cy = earthRadiusKm * (1 - krefract) + aircraftAltKmMsl;
      double phi = slantDistKm / refractRadiusKm;
      double gamma = M_PI / 2 + elevRad - phi;
      double px = cx + refractRadiusKm * cos( gamma);
      double py = cy + refractRadiusKm * sin( gamma);
      double altKmMsl = sqrt( px*px + py*py) - earthRadiusKm;  // hp
      double zeta = atan2( py, px);
      double horizDistKm = earthRadiusKm * (M_PI/2 - zeta);    // se

      // Calc error of the simple version
      double errHoriz = flatHorizDistKm - horizDistKm;
      double errVert = flatAltKmMsl - altKmMsl;

      if (fabs(errHoriz) > *maxAbsErrHoriz)
        *maxAbsErrHoriz = fabs(errHoriz);
      if (fabs(errVert) > *maxAbsErrVert)
        *maxAbsErrVert = fabs(errVert);

      if (bugs >= Params::DEBUG_VERBOSE) {
        std::cout << setprecision(15);
        std::cout << "readRadarFile calcs:" << std::endl;
        std::cout << "  Slant dist: " << slantDistKm
          << "  elevDeg: " << elevDeg << std::endl;
        std::cout << "  cx: " << cx << "  cy: " << cy << std::endl;
        std::cout << "  px: " << px << "  py: " << py << std::endl;
        std::cout << "  elevRad: " << elevRad << std::endl;
        std::cout << "  phi: " << phi << std::endl;
        std::cout << "  zeta: " << zeta << std::endl;
        std::cout << "  Horiz dist: flatHorizDistKm: " << flatHorizDistKm
          << "  with refract: horizDistKm: " << horizDistKm
          << "  err: " << errHoriz
          << "  maxAbs: " << *maxAbsErrHoriz << std::endl;
        std::cout << "  flatAltKmMsl: " << flatAltKmMsl
          << "  with refract: altKmMsl: " << altKmMsl
          << "  err: " << errVert
          << "  maxAbs: " << *maxAbsErrVert << std::endl;
      }

      double valVg = fieldVg->getDoubleValue(ipt);
      double valDbz = fieldDbz->getDoubleValue(ipt);
      double valNcp = fieldNcp->getDoubleValue(ipt);

      if (bugs >= Params::DEBUG_VERBOSE)
	std::cout << "valVg: " << valVg << ", valDbz: " << valDbz << ", valNcp: " << valNcp << std::endl;

      cntTimej++;
      addDeltaTime( &timej, &sumTimej);
      struct timeval timek;
      addDeltaTime( &timek, NULL);

      // In one test of 237 files, this is called 28766400 times
      // with a total time of 31.8 seconds, or 1.1e-6 seconds per call.
      double latDeg, lonDeg, azi2Deg, m12;
      double a12 = geodesic->Direct(
        aircraftLatDeg, aircraftLonDeg,
        aziDeg, 1000 * horizDistKm,
        latDeg, lonDeg, azi2Deg, m12);

      cntTimek++;
      addDeltaTime( &timek, &sumTimek);
      struct timeval timem;
      addDeltaTime( &timem, NULL);

      // Make sure the data are valid ...
      // before going to the expense of finding the y,x coords.
      // If not missing and ...

      bool usePoint = false;
      if ( valVg != missVg
        && valDbz != missDbz && (forceOk || valDbz >= minDbz)
        && valNcp != missNcp && (forceOk || valNcp >= minNcp)
        && (maxDeltaAltKm == 0
          || fabs( altKmMsl - aircraftAltKmMsl) <= maxDeltaAltKm)
        && (maxAbsElevDeg == 0 || fabs( elevDeg) <= maxAbsElevDeg)
        && (minRadialDistKm == 0 || slantDistKm >= minRadialDistKm))
      {

        struct timeval timel;
        addDeltaTime( &timel, NULL);

        double coordy;       // northing in meters
        double coordx;       // easting in meters

        // Caution: the following call is slow!
        // In one test of 237 files, it is called 28766400 times
        // with a total time of 180.9 seconds,
        // or 6.3e-6 seconds per call.

        latLonToYX(
          projLon0,          // central meridian of projection
          basey,             // coord y base, km
          basex,             // coord x base, km
          latDeg,
          lonDeg,
          coordy,            // output value
          coordx);           // output value

        if (bugs >= Params::DEBUG_VERBOSE) {
          std::cout << setprecision(15);
          std::cout << "    add ob: altKmMsl: " << altKmMsl
            << "  coordx: " << coordx
            << "  coordy: " << coordy << std::endl;
        }

        // Convert azimuth degrees to polar coordinates radians
        double thetaRad = 0.5 * M_PI - aziDeg * M_PI / 180;
        while (thetaRad < 0) thetaRad += 2*M_PI;
        while (thetaRad > 2*M_PI) thetaRad -= 2*M_PI;

        cntTimel++;
        addDeltaTime( &timel, &sumTimel);

        // Make sure coordz,y,x are within the grid
        double coordz = altKmMsl;
        if (
             (zgridmin==MISS_PARM || coordz >= zgridmin && coordz < zgridmax)
          && (ygridmin==MISS_PARM || coordy >= ygridmin && coordy < ygridmax)
          && (xgridmin==MISS_PARM || coordx >= xgridmin && coordx < xgridmax))
        {
          usePoint = true;
          numGoodPoint++;
          bool showDetail = testDetail(
            coordz,            // z
            coordy,            // y
            coordx             // x
				       );
          if (showDetail && (bugs >= Params::DEBUG_NORM)) {
            std::cout << "readRadarFile: showDetail: coordz: " << coordz << std::endl;
            std::cout << "readRadarFile: showDetail: coordy: " << coordy << std::endl;
            std::cout << "readRadarFile: showDetail: coordx: " << coordx << std::endl;
          }

          pointBbox->addOb( coordz, coordy, coordx);

	  // testModes explained in Args.cc

          if ((testMode != Params::MODE_ZETA) && (testMode != Params::MODE_GAMMA)) {

	    // Use synthetic wind data

            double synVels[3];    // W, V, U

            if (testMode == Params::MODE_BETA) {
              calcSyntheticWinds(
                showDetail,
                coordz,             // locz
                coordy,             // locy
                coordx,             // locx
                synVels);           // returned W, V, U
            }
            else if (testMode == Params::MODE_ZETA_BELTRAMI) {
              calcBeltramiFlow(
                showDetail,
                coordz,             // locz
                coordy,             // locy
                coordx,             // locx
                synVels);           // returned W, V, U
            }

            valVg = calcRadialVelocity(
              showDetail,
              thetaRad,           // polar coord angle from observer, radians
              elevRad,            // elevation angle from observer, radians
              synVels);           // W, V, U

            valDbz = 30;
            valNcp = 1;
          } // if ! ZETA

          if (fabs(valVg) > 200) throwerr("invalid valVg: %g", valVg);
          if (fabs(valDbz) > 200) throwerr("invalid valDbz: %g", valDbz);
#if 0
          if (valNcp < 0 || valNcp > 1) throwerr("invalid valNcp: %g", valNcp);
#endif
          pointVec->push_back( new Point(
            ifile,
            iray,
            ipt,
            aircraftAltKmMsl,
            aircrafty,
            aircraftx,
            altKmMsl,
            coordy,
            coordx,
            ray->getTimeDouble(),      // UTC seconds since 1970 jan 1.
            thetaRad,                  // polar coords angle
            elevRad,                   // angle from the aircraft between horizontal and pt
            valVg,                     // radial velocity
            valDbz,                    // log reflectivity
            valNcp));                  // radar net coherent power

          statVg->addOb( valVg);
          statDbz->addOb( valDbz);
          statNcp->addOb( valNcp);

          if (bugs >= Params::DEBUG_EXTRA && itotpt % 10000 == 0) {
            std::cout << setprecision(5);
            std::cout << "readRadarFile: ok:"
              << "  iray: " << iray
              << "  ipt: " << ipt
              << "  itotpt: " << itotpt << std::endl
              << "  elevDeg: " << elevDeg
              << "  altKmMsl: " << altKmMsl
              << "  coordy: " << coordy
              << "  coordx: " << coordx << std::endl
              << "  valVg: " << valVg
              << "  valDbz: " << valDbz
              << "  valNcp: " << valNcp << std::endl;
          }
        } // if point is within the grid
        else numMissLocPoint++;
      } // if point is valid
      else numInvalidPoint++;

      cntTimem++;
      addDeltaTime( &timem, &sumTimem);
      struct timeval timen;
      addDeltaTime( &timen, NULL);

      cntTimen++;
      addDeltaTime( &timen, &sumTimen);

      itotpt++;
    } // for ipt
    cntTimei++;
    addDeltaTime( &timei, &sumTimei);

  } // for iray

  cntTimeh++;
  addDeltaTime( &timeh, &sumTimeh);

  std::cout << "readRadarFile: numGoodPoint: " << numGoodPoint << std::endl;
  std::cout << "readRadarFile: numInvalidPoint: " << numInvalidPoint << std::endl;
  std::cout << "readRadarFile: numMissLocPoint: " << numMissLocPoint << std::endl;
} // end readRadarFile
