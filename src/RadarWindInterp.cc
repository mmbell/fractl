#include "Interps.hh"

RadarWindInterp::RadarWindInterp()
{
}

RadarWindInterp::~RadarWindInterp()
{
}

bool RadarWindInterp::interpolate(double *y, int n1, int n2, int n3, double flag)
{
  std::cout << "RadarWind original interpolation is not supported yet" << std::endl;
  return false;
}

// TODO
// Integrate this original code into the Interp class hierarchy

// NO LONGER USED:
// Interpolate missing values from neighbors.
//
// The problem is that some cells don't have enough valid
// neighbors to calculate their winds, so their values are
// left as NaN.
//
// The idea was to use fill the NaN cells
// by interpolating the wind values from nearby cells.
//
// But what happens is that generally if a cell doesn't have
// enough good neighbors for a valid wind calculation,
// any neighbor cells we might use for interpolation are
// pretty flakey.

#if 0

void RadarWind::interpMissing(
			      
  long bugs,
  long nradz,                  // grid z dim
  long nrady,                  // grid y dim
  long nradx,                  // grid x dim
  double ***& radmat,          // interpolate to fill NaN values
  double ***& tmpmat)          // temp work
  
  //xxx make this a parm
  char nbrType = 'f';          // f: faces, c: corners + faces

  long nnbr = 0;
  // Parallel arrays:
  long * izs;         // offsets in the z direction
  long * iys;         // offsets in the y direction
  long * ixs;         // offsets in the x direction

  if (nbrType == 'f') {     // nbrs at faces only
    nnbr = 6;
    izs = new long[ nnbr];
    iys = new long[ nnbr];
    ixs = new long[ nnbr];
    long ii = 0;
    // Each row gives the offsets to a neighbor.
    izs[ii] = -1;  iys[ii] =  0;  ixs[ii] =  0;  ii++;   // nbr below
    izs[ii] =  0;  iys[ii] = -1;  ixs[ii] =  0;  ii++;   // nbrs at same z
    izs[ii] =  0;  iys[ii] =  0;  ixs[ii] = -1;  ii++;
    izs[ii] =  0;  iys[ii] =  0;  ixs[ii] =  1;  ii++;
    izs[ii] =  0;  iys[ii] =  1;  ixs[ii] =  0;  ii++;
    izs[ii] =  1;  iys[ii] =  0;  ixs[ii] =  0;  ii++;   // nbr above
  }

  else if (nbrType == 'c') {     // nbrs at corners + faces
    nnbr = 26;
    izs = new long[ nnbr];
    iys = new long[ nnbr];
    ixs = new long[ nnbr];
    long ii = 0;
    // Each row gives the offsets to a neighbor.
    izs[ii] = -1;  iys[ii] = -1;  ixs[ii] = -1;  ii++;   // nbrs below
    izs[ii] = -1;  iys[ii] = -1;  ixs[ii] =  0;  ii++;
    izs[ii] = -1;  iys[ii] = -1;  ixs[ii] =  1;  ii++;

    izs[ii] = -1;  iys[ii] =  0;  ixs[ii] = -1;  ii++;
    izs[ii] = -1;  iys[ii] =  0;  ixs[ii] =  0;  ii++;
    izs[ii] = -1;  iys[ii] =  0;  ixs[ii] =  1;  ii++;

    izs[ii] = -1;  iys[ii] =  1;  ixs[ii] = -1;  ii++;
    izs[ii] = -1;  iys[ii] =  1;  ixs[ii] =  0;  ii++;
    izs[ii] = -1;  iys[ii] =  1;  ixs[ii] =  1;  ii++;

    izs[ii] =  0;  iys[ii] = -1;  ixs[ii] = -1;  ii++;   // nbrs at same z
    izs[ii] =  0;  iys[ii] = -1;  ixs[ii] =  0;  ii++;
    izs[ii] =  0;  iys[ii] = -1;  ixs[ii] =  1;  ii++;

    izs[ii] =  0;  iys[ii] =  0;  ixs[ii] = -1;  ii++;   // nbrs at same z,y
    izs[ii] =  0;  iys[ii] =  0;  ixs[ii] =  1;  ii++;

    izs[ii] =  0;  iys[ii] =  1;  ixs[ii] = -1;  ii++;
    izs[ii] =  0;  iys[ii] =  1;  ixs[ii] =  0;  ii++;
    izs[ii] =  0;  iys[ii] =  1;  ixs[ii] =  1;  ii++;

    izs[ii] =  1;  iys[ii] = -1;  ixs[ii] = -1;  ii++;   // nbrs above
    izs[ii] =  1;  iys[ii] = -1;  ixs[ii] =  0;  ii++;
    izs[ii] =  1;  iys[ii] = -1;  ixs[ii] =  1;  ii++;

    izs[ii] =  1;  iys[ii] =  0;  ixs[ii] = -1;  ii++;
    izs[ii] =  1;  iys[ii] =  0;  ixs[ii] =  0;  ii++;
    izs[ii] =  1;  iys[ii] =  0;  ixs[ii] =  1;  ii++;

    izs[ii] =  1;  iys[ii] =  1;  ixs[ii] = -1;  ii++;
    izs[ii] =  1;  iys[ii] =  1;  ixs[ii] =  0;  ii++;
    izs[ii] =  1;  iys[ii] =  1;  ixs[ii] =  1;  ii++;
  }
  else throwerr("invalid nbrType");


  // Copy radmat to tmpmat, substituting interpolated values
  // for non-ok cells.
  
  for (long iz = 0; iz < nradz; iz++) {
    for (long iy = 0; iy < nrady; iy++) {
      for (long ix = 0; ix < nradx; ix++) {

        if (isOkDouble( radmat[iz][iy][ix])) {
          tmpmat[iz][iy][ix] = radmat[iz][iy][ix];
        }

        else {
          long numval = 0;
          double tsum = 0;

          for (long inb = 0; inb < nnbr; inb++) {
            // Set jz,jy,jx = indices of the neighbor.
            long jz = iz + izs[inb];
            long jy = iy + iys[inb];
            long jx = ix + ixs[inb];
            if ( jz >= 0 && jz < nradz
		 && jy >= 0 && jy < nrady
		 && jx >= 0 && jx < nradx
		 && isOkDouble( radmat[jz][jy][jx]))
	      {
		numval++;
		tsum += radmat[jz][jy][jx];
	      }
            if (numval > 0)
              tmpmat[iz][iy][ix] = tsum / numval;
            else
              tmpmat[iz][iy][ix] = numeric_limits<double>::quiet_NaN();
          }

        } // else not ok
      } // for ix
    } // for iy
  } // for iz

  // Copy tmpmat back to radmat
  for (long iz = 0; iz < nradz; iz++) {
    for (long iy = 0; iy < nrady; iy++) {
      for (long ix = 0; ix < nradx; ix++) {
        radmat[iz][iy][ix] = tmpmat[iz][iy][ix];
      } // for ix
    } // for iy
  } // for iz

} // end interpMissing

#endif
