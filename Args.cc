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

#include "FileSpec.hh"
#include "Fractl.hh"
#include "Params.hh"

#include <cstdarg>
#include <iostream>
#include <fstream>
#include <regex>
#include <unordered_set>


bool Fractl::parseArgs(int argc, char *argv[])
{
  std::unordered_set<std::string> keywordSet( {
      "minDbz", "minNcp", "testMode", "synWinds",
	"conditionNumberCutoff", "maxU", "maxV", "maxW",
	"radFiles", "zGrid", "yGrid", "xGrid", "gridType", "projName",
	"projLat0", "projLon0", "baseW", "epsilon",
	"maxDeltaAltKm", "maxAbsElevDeg", "minRadialDistKm",
	"mumNbrMax", "maxDistBase", "maxDistFactor", "forceOk",
	"useEigen", "preGridded", "inDir", "fileRegex", "fileList",
	"radialName", "dbzName", "ncpName", "outText", "outNc",
	"detailSpec", "radarAlt", "gridType",
	"uvFilter", "wFilter", "uvSteps", "uvMultiStep",
	"wStep", "wMultiStep" "uvInterp", "windComputationMethod",
	// MASS2
	"direction", "alpha", "integrationType", "lowerBoundaryInitMethod",
	"initVal", "errMax", "iterMax",
	"zLimitUnits", "zLowerLevel", "zUpperLevel", "zLastLevel",
	"constA", "constB", "consC"
	} );

  tdrp_override_t override;
  TDRP_init_override(&override);

  char tmp_str[BUFSIZ];

  // fill override from command line arguments
  
  for (int i = 1; i < argc; i++) {
    const char *arg = argv[i] + 1;
    std::unordered_set<std::string>::const_iterator found = keywordSet.find(arg);
    if ( found != keywordSet.end() ) {
      sprintf(tmp_str, "%s = %s;", arg, argv[++i]);
      TDRP_add_override(&override, tmp_str);      
    } else {
      // This one has specific meaning to TDRP.
      if ( ! strcmp(argv[i], "-params" ) )
	i++;
    }
  }

  // Now parse the config file
  
  Params params;
  char *paramsPath;
  
  if (params.loadFromArgs(argc, argv, override.list, &paramsPath))
    badparms("Problem with command line args");

  TDRP_free_override(&override);

  // and finally process the arguments

  bugs = params.debug;
  minDbz = params.minDbz;
  minNcp = params.minNcp;
  testMode = params.testMode;
  gridType = params.gridType;
  
  synWinds = new double[3];
  if (! parseSynWinds(params.synWinds, synWinds))
    badparms("Invalid synWinds");

  if ( strcmp(params.radFiles, "not_set")) {
    radFiles = new long[2];
    if (! parseRadFiles(params.radFiles, radFiles))
      badparms("Invalid radFiles");
  }

  if (! parseGridSpec(params.zGrid, zgridmin, zgridmax, zgridinc) )
    badparms("Invalid zgrid");

  if (! parseGridSpec(params.yGrid, ygridmin, ygridmax, ygridinc) )
    badparms("Invalid ygrid");

  if (! parseGridSpec(params.xGrid, xgridmin, xgridmax, xgridinc) )
    badparms("Invalid xgrid");

  projName = params.projName;
  projLat0 = params.projLat0;
  projLon0 = params.projLon0;
  radarAlt = params.radarAlt;
  
  baseW = params.baseW;
  epsilon = params.epsilon;

  maxDeltaAltKm = params.maxDeltaAltKm;
  maxAbsElevDeg = params.maxAbsElevDeg;
  minRadialDistKm = params.minRadialDistKm;
  numNbrMax = params.numNbrMax;
  maxDistBase = params.maxDistBase;
  maxDistFactor = params.maxDistFactor;

  conditionNumberCutoff = params.conditionNumberCutoff;
  maxU = params.maxU;
  maxV = params.maxV;
  maxW = params.maxW;

  forceOk = parseBool("forceOk", params.forceOk);
  useEigen = parseBool("useEigen", params.useEigen);
  preGridded = parseBool("preGridded", params.preGridded);

  if (strcmp(params.inDir, "not_set"))
    inDir = params.inDir;
  if (strcmp(params.fileList, "not_set"))  
    fileList = params.fileList;
  if (strcmp(params.fileRegex, "not_set"))    
    fileRegex = params.fileRegex;
  
  radialName = params.radialName;
  dbzName = params.dbzName;
  ncpName = params.ncpName;

  if (strcmp(params.outTxt, "not_set"))      
    outTxt = params.outTxt;
  if (strcmp(params.outNc, "not_set"))        
    outNc = params.outNc;

  uvFilter = params.uvFilter;
  wFilter = params.wFilter;
  uvSteps = params.uvSteps;
  wSteps = params.wSteps;
  uvInterp = params.uvInterp;

  uvMultiStep = parseMultiStep(params.uvMultiStep);
  wMultiStep = parseMultiStep(params.wMultiStep);  
  
  if (strcmp(params.detailSpec, "not_set")) {
    detailSpec = new double[4];
    if (! parseDetailSpec(params.detailSpec, detailSpec))
      badparms("Invalid detailSpec");
  }

  if(verbose)
    dumpArgs();

  windComputationMethod = params.windComputationMethod;
  
  // Cedric MASS2

  mas2ErrMax = params.errMax;
  mas2IterMax = params.iterMax;
  mas2Dir = params.direction;
  lowerBoundaryInitMethod = params.lowerBoundaryInitMethod;
  initStr = params.initVal;
  
  vtA = params.constA;
  vtB = params.constB;
  vtC = params.constC;

  return checkArgs();
}

bool Fractl::parseFile(const char *fileName, int &argc, char ***argv)
{

  ifstream myFile(fileName);
  
  if (myFile.is_open()) {
  } else {
    std::cerr << "Failed to open config file " << fileName << std::endl;
    return false;
  }

  argc = 0;
  string line;
  std::vector<std::string> myArgv;

  myArgv.push_back("Fractl::parseFile");
  
  std::regex empty("^\\s*$");
  std::regex kv("^\\s*([^\\s]+)\\s*([^\\s]+)");
  std::smatch match;
  
  while (std::getline(myFile, line) ) {
    if (std::regex_match(line, empty))
      continue;
    if (std::regex_search(line, match, kv)) {
      argc += 2;
      myArgv.push_back(match.str(1));
      myArgv.push_back(match.str(2));
    }
  }
  if (argc == 0)
    return false;

  argc += 1;
  char **cp = (char **) malloc( sizeof (char *) * argc);
  *argv = cp;
  for(int i = 0; i < argc; i++, cp++)
    *cp = strdup(myArgv[i].c_str());
  
  // convert vector of strings to array of strings
  return true;
}

// Print error message and usage info, and exit 1.

void Fractl::badparms( const string msg, ...) {
  va_list vaList;
  va_start( vaList, msg);
  long bufLen = 10000;
  char * buf = new char[bufLen];
  vsnprintf( buf, bufLen, msg.c_str(), vaList);
  va_end( vaList);

  cout << "Fractl: error:" << endl;
  cout << buf << endl;
  cout << "Parameters:\n" << endl;
  cout << "  -bugs           debug level" << endl;
  cout << endl;
  cout << "  -testMode       test mode:" << endl;
  cout << endl;
  cout << "     MODE_ALPHA:" << endl;
  cout << "       Use simple synthetic wind data." << endl;
  cout << "       Skip the radar files." << endl;
  cout << "       Must specify synWinds: W,V,U" << endl;
  cout << "       Obs are located at the verif cell centers." << endl;
  cout << endl;
  cout << "     MODE_BETA:" << endl;
  cout << "       Use simple synthetic wind data." << endl;
  cout << "       Read swp files and use swp obs locations," << endl;
  cout << "       but replace radial vel with synthetic." << endl;
  cout << "       Must specify synWinds: W,V,U" << endl;
  cout << "       radFiles beg,lim" << endl;
  cout << "       Use radar files i: beg <= i < lim." << endl;
  cout << "       If beg==lim==0, read all files" << endl;
  cout << endl;
  cout << "     MODE_ZETA:" << endl;
  cout << "       Operational: use specified radar data." << endl;
  cout << "       Specify synWinds: 0,0,0" << endl;
  cout << "       radFiles beg,lim" << endl;
  cout << "       Use radar files i: beg <= i < lim." << endl;
  cout << "       If beg==lim==0, read all files" << endl;
  cout << endl;
  cout << "     MODE_ZETA_BELTRAMI:" << endl;
  cout << "       Like zeta, but compare the results with" << endl;
  cout << "       the test Beltrami flow and print statistics." << endl;
  cout << "       Specify synWinds: 0,0,0" << endl;
  cout << "       This is only useful for swp files generated" << endl;
  cout << "       by Michael Bell's Beltrami flow simulation." << endl;
  cout << endl;
  cout << endl;
  cout << "  -synWinds       Specify wind field for any testMode" << endl;
  cout << "                    other than zeta." << endl;
  cout << "                    Specify a comma sep triplet: W,V,U" << endl;
  cout << "                    If a wind spec is numeric," << endl;
  cout << "                      use that constant uniform value." << endl;
  cout << "                    Else the wind spec is the name for" << endl;
  cout << "                      one of SYNFUNC_* values:" << endl;
  cout << "                      \"sinx\": wind component = sin(locx)" << endl;
  cout << "                      \"siny\": wind component = sin(locy)" << endl;
  cout << "                      \"sinz\": wind component = sin(locz)" << endl;
  cout << "                    Example:" << endl;
  cout << "                      -synWinds 3,sinx,sinx" << endl;
  cout << "                      Means uniform Z wind at 3 m/s." << endl;
  cout << "                      Both V and U winds = sin(locx)," << endl;
  cout << "                        so the wind would be from the SW" << endl;
  cout << "                        to the NW, varying as sin(locx)." << endl;
  cout << endl;
  cout << "  -zGrid          z grid spec: min,max,inc" << endl;
  cout << "                    To get min,max from the data," << endl;
  cout << "                    specify only one value: increment" << endl;
  cout << "  -yGrid          like zgrid" << endl;
  cout << "  -xGrid          like zgrid" << endl;
  cout << endl;
  cout << "  -projName       projection name: must be transverseMercator" << endl;
  cout << "  -projLat0       projection lat0.  For example 16.5" << endl;
  cout << "  -projLon0       projection lon0.  For example 148.0" << endl;
  cout << endl;
  cout << "  -baseW          W wind into the base of the lowest layer" << endl;
  cout << "  -epsilon        epsilon" << endl;
  cout << endl;
  cout << "  -maxDeltaAltKm  max abs delta altitude of observations" << endl;
  cout << "                  from aircraft, km, or 0";
  cout << endl;
  cout << "  -maxAbsElevDeg  max abs elevation angle observations," << endl;
  cout << "                  degrees, or 0";
  cout << endl;
  cout << "  -minRadialDistKm  min radial distance of observations" << endl;
  cout << "                  from aircraft, km, or 0";
  cout << endl;
  cout << "  -maxU           Maximum value for U wind component" << endl;
  cout << "  -maxV           Maximum value for V wind component" << endl;
  cout << "  -maxU           Maximum value for W wind component" << endl;
  cout << "  -conditionNumberCutoff    Maximum value for a cell ConditionNumber" << endl;
  cout << endl;  
  cout << "  -numNbrMax      max num nearest nbrs" << endl;
  cout << "  -maxDistBase    max pt dist = base + factor*aircraftDist" << endl;
  cout << "  -maxDistFactor  max pt dist = base + factor*aircraftDist" << endl;
  cout << endl;
  cout << "  -forceOk        y/n: force ncp and dbz to ok on all points" << endl;
  cout << "  -useEigen       y/n: y: use Eigen.  n: use Cramer" << endl;
  cout << endl;
  cout << "  -preGridded     y/n: Data is pre gridded and doesn't need projection" << endl;
  cout << "  -inDir          dir with radx input files" << endl;
  cout << "                  Used airborne radars." << endl;
  cout << "                  Mutually exclusive with -fileList." << endl;
  cout << "  -fileRegex      regex for files in inDir, like '^swp'" << endl;
  cout << "                  Mutually exclusive with -fileList." << endl;

  cout << endl;
  cout << "  -fileList       file containing list of input files." << endl;
  cout << "                  Used for fixed location radars." << endl;
  cout << "                  Mutually exclusive with -inDir, -fileRegex." << endl;
  cout << "                  Format is one entry per line." << endl;
  cout << "                  # in the first column starts a comment." << endl;
  cout << "                  Each line has the format:" << endl;
  cout << "                    fileName altKmMsl latDeg lonDeg" << endl;
  cout << "                   where altKmMsl latDeg lonDeg" << endl;
  cout << "                   give the location of the fixed radar." << endl;
  cout << endl;
  cout << "  -radialName     Name of the radial velocity field" << endl;
  cout << "                  in the input files." << endl;
  cout << "                  For Eldora, this is VG." << endl;
  cout << "                  For CSU-Chill, this is VE." << endl;
  cout << endl;
  cout << "  -dbzName        Name of the reflectivity field" << endl;
  cout << "                  in the input files." << endl;
  cout << "                  For Eldora, this is DBZ." << endl;
  cout << "                  For CSU-Chill, this is DZ." << endl;
  cout << endl;
  cout << "  -ncpName        Name of the net coherent power field" << endl;
  cout << "                  in the input files." << endl;
  cout << "                  For Eldora, this is NCP." << endl;
  cout << "                  For CSU-Chill, this is NC." << endl;
  cout << endl;
  cout << "  -outTxt         output text file:" << endl;
  cout << "                    has verif or grid results." << endl;
  cout << endl;
  cout << "  -outNc          output netcdf file" << endl;
  cout << "                    If outNc ends in a slash, \"x/\"," << endl;
  cout << "                    it is a directory." << endl;
  cout << "                    We make a subdir and write to:" << endl;
  cout << "                    x/yyyymmdd/ncf_yyyymmdd_hhmmss.nc" << endl;
  cout << endl;
  cout << "  -detailSpec     z,y,x,radius     (optional)" << endl;
  cout << "                  If none, it's the same as omitting it." << endl;
  cout << "                  Specifies a region for detailed logging." << endl;
  cout << "                  Used for debugging only." << endl;
  cout << endl;
  cout << "  -uvFilter       Filtering method to use on U and V" << endl;
  cout << "                     Default FILTER_OFF" << endl;
  cout << endl;
  cout << "  -wFilter        Filtering method to use on W" << endl;
  cout << "                     Default FILTER_OFF" << endl;
  cout << endl;
  cout << "  -uvStep         How many steps to use in the filtering" << endl;
  cout << "                     This is used for all dimensions" << endl;
  cout << "                     Use uvMultiStep to specify individual steps instead" << endl;
  cout << endl;
  cout << "  -wStep         How many steps to use in the filtering" << endl;
  cout << "                     This is used for all dimensions" << endl;
  cout << "                     Use wMultiStep to specify individual steps instead" << endl;
  cout << endl;
  cout << "  -uvMultiStep   How many steps to use in the filtering" << endl;
  cout << "                     Enter one value per dimension for example \"1,2,3\"" << endl;
  cout << endl;
  cout << "  -wMultiStep    How many steps to use in the filtering" << endl;
  cout << "                     Enter one value per dimension for example \"1,2,3\"" << endl;
  cout << endl;
  cout << " -uvInterp       What method to use to interpolate missing U and V values" << endl;
  cout << "                     Default INTERP_NONE" << endl;
  cout << endl;

  exit(1);
}

bool Fractl::parseSynWinds(char *param, double *synWinds)
{
  vector<std::string> tokvec;
  splitString( param, ",", tokvec);
  if (tokvec.size() != 3)
    return false;
  
  for (int ii = 0; ii < 3; ii++) {
    std::string tok = tokvec.at(ii);
    if (isdigit(tok.at(0)) || tok.at(0) == '-')
      synWinds[ii] = parseDouble("synWinds", tokvec.at(ii));
    else if (tokvec.at(ii) == "sinz") synWinds[ii] = SYNFUNC_SINZ;
    else if (tokvec.at(ii) == "siny") synWinds[ii] = SYNFUNC_SINY;
    else if (tokvec.at(ii) == "sinx") synWinds[ii] = SYNFUNC_SINX;
    else throwerr("unknown synWinds function");
  }
  return true;
}

bool Fractl::parseRadFiles(char *param, long *radFiles)
{
  vector<std::string> tokvec;
  splitString(param, ",", tokvec);
  if (tokvec.size() != 2)
    return false;
  for (int ii = 0; ii < 2; ii++) {
    std::string tok = tokvec.at(ii);
    radFiles[ii] = parseLong("radFiles", tokvec.at(ii));
  }
  return true;
}

bool Fractl::parseGridSpec(char *param, double &min, double &max, double &step)
{
  vector<std::string> tokvec;
  if ( ! strcmp(param, "not_set")) {
    return false;
  }
  splitString(param, ",", tokvec);
  if (tokvec.size() == 1) {
    step = parseDouble("step", tokvec.at(0));
  }
  else if (tokvec.size() == 3) {
    min = parseDouble("min", tokvec.at(0));
    max = parseDouble("max", tokvec.at(1));
    step = parseDouble("step", tokvec.at(2));
  }
  else
    return false;
  return true;
}

bool Fractl::parseDetailSpec(char *param, double *spec)
{
  vector<std::string> tokvec;
  splitString(param, ",", tokvec);
  if (tokvec.size() != 4)
    return false;

  for (int ii = 0; ii < 4; ii++) {
    spec[ii] = parseDouble("detailSpec", tokvec.at(ii));
  }
  return true;
}
  
int *Fractl::parseMultiStep(char *param)
{
  if ( ! strcmp(param, "not_set"))
    return NULL;
  vector<std::string> tokvec;
  splitString(param, ",", tokvec);
  int nitems = tokvec.size();
  if(nitems > 3) {
    std::cerr << "A maximum of 3 item is expected for multiSpec. Spec ignored" << std::endl;
    return NULL;
  }
  
  int *retval = new int[3]();
  
  for (int ii = 0; ii < nitems; ii++) {
    retval[ii] = parseLong("multiStep", tokvec.at(ii));
  }
  return retval;
}
