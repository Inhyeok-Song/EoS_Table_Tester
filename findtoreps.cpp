#include "NuclearEoS.h"


void findtoreps_bdry( const real x, const real y, const real z,
                      real *found_ltoreps, const real *alltables_mode,
                      const int nx, const int ny, const int nz, const int ntoreps,
                      const real *xt, const real *yt, const real *zt,
                      const real *logtoreps, const int keymode, int *keyerr );



//-------------------------------------------------------------------------------------
// Function    :  findtoreps
// Description :  Finding temperature from different modes
//                -->                 energy mode   (0)
//                                    entropy mode  (2)
//                                    pressure mode (3)
//
// Note        :  1. Use 3D Catmull-Rom cubic interpolation formula
//                   to search the corresponding temperature/energy given 
//                   (rho, (eps, T, entropy, P), Y_e)
//                2. Invoked by nuc_eos_C_short()
//
// Parameter   :  x              : Input vector of first  variable (rho)
//                y              : Input vector of second variable (eps, T, entropy, P)
//                z              : Input vector of third  variable (Y_e)
//                found_ltoreps  : Output log(temp)/log(energy) of interpolated function values
//                alltables_mode : 3d array of tabulated logtemp
//                nx             : X-dimension of table
//                ny             : Y-dimension of table
//                nz             : Z-dimension of table
//                ntoreps        : Size of (temperature/energy) array in the Nuclear Eos table
//                xt             : Vector of x-coordinates of table
//                yt             : Vector of y-coordinates of table
//                zt             : Vector of z-coordinates of table
//                logtoreps      : log(T)/log(energy) array in the table
//                interpol_TL    : interpolation schemes for table look-ups (linear/cubic)
//                keymode        : Which mode we will use
//                                 --> 1: Energy mode   (coming in with internal energy)
//                                     2: Entropy mode  (coming in with entropy)
//                                     3: Pressure mode (coming in with P)
//                keyerr         : Output error
//
// Return      :  found_ltoreps
//-------------------------------------------------------------------------------------
void findtoreps( const real x, const real y, const real z,
                 real *found_ltoreps, const real *alltables_mode,
                 const int nx, const int ny, const int nz, const int ntoreps,
                 const real *xt, const real *yt, const real *zt, const real *logtoreps,
                 const int interpol_TL, const int keymode, int *keyerr )
{


   if      ( interpol_TL == NUC_INTERPOL_LINEAR )
   {
      findtoreps_bdry( x, y, z, found_ltoreps, alltables_mode, nx, ny, nz, ntoreps,
                       xt, yt, zt, logtoreps, keymode, keyerr );
      return;
   }
   else if ( interpol_TL == NUC_INTERPOL_CUBIC  )
   {
      const real *pv = NULL;

      real   dx, dy, dz, dxi, dyi, dzi;
      real   delx, dely, delz;
      real   u[4], v[4], w[4];
      real   r[4], q[4];
      real   vox;
      int    ix, iy, iz;
      int    nxy;

      nxy = nx*ny;

// determine spacing parameters of equidistant (!!!) table
#     if 1
      dx = ( xt[nx-1] - xt[0] ) / (real)(nx-1);
      dy = ( yt[ny-1] - yt[0] ) / (real)(ny-1);
      dz = ( zt[nz-1] - zt[0] ) / (real)(nz-1);

      dxi = (real)1.0 / dx;
      dyi = (real)1.0 / dy;
      dzi = (real)1.0 / dz;
#     endif

#     if 0
      dx = drho;
      dy = dtemp;
      dz = dye;

      dxi = drhoi;
      dyi = dtempi;
      dzi = dyei;
#     endif


      // determine location in table
      ix = (int)( (x - xt[0] )*dxi );
      iy = (int)( (y - yt[0] )*dyi );
      iz = (int)( (z - zt[0] )*dzi );


      // linear interpolation at boundaries
      if ( ix == 0    || iy == 0    || iz == 0 ||
           ix >= nx-2 || iy >= ny-2 || iz >= nz-2 )
      {
         findtoreps_bdry( x, y, z, found_ltoreps, alltables_mode, nx, ny, nz, ntoreps,
                          xt, yt, zt, logtoreps, keymode, keyerr );
         return;
      }


// differences
      delx = ( x - xt[ix] )*dxi;
      dely = ( y - yt[iy] )*dyi;
      delz = ( z - zt[iz] )*dzi;


// factors for Catmull-Rom interpolation
      const real delx2 =  SQR( delx );
      const real dely2 =  SQR( dely );
      const real delz2 =  SQR( delz );
      const real delx3 = CUBE( delx );
      const real dely3 = CUBE( dely );
      const real delz3 = CUBE( delz );

      u[0] = (real)-0.5*delx3 +           delx2 - (real)0.5*delx;
      u[1] = (real) 1.5*delx3 - (real)2.5*delx2 + (real)1.0;
      u[2] = (real)-1.5*delx3 + (real)2.0*delx2 + (real)0.5*delx;
      u[3] = (real) 0.5*delx3 - (real)0.5*delx2;

      v[0] = (real)-0.5*dely3 +           dely2 - (real)0.5*dely;
      v[1] = (real) 1.5*dely3 - (real)2.5*dely2 + (real)1.0;
      v[2] = (real)-1.5*dely3 + (real)2.0*dely2 + (real)0.5*dely;
      v[3] = (real) 0.5*dely3 - (real)0.5*dely2;

      w[0] = (real)-0.5*delz3 +           delz2 - (real)0.5*delz;
      w[1] = (real) 1.5*delz3 - (real)2.5*delz2 + (real)1.0;
      w[2] = (real)-1.5*delz3 + (real)2.0*delz2 + (real)0.5*delz;
      w[3] = (real) 0.5*delz3 - (real)0.5*delz2;


      int iv;
#     if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
      if      ( keymode == NUC_MODE_ENGY ) iv = 0; // temperature table for the energy mode
#     elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
      if      ( keymode == NUC_MODE_TEMP ) iv = 0; // energy table for the temperautre mode
#     endif // #elif NUC_TABLE_MODE ... else ...
      else if ( keymode == NUC_MODE_ENTR ) iv = 1; // temperature/energy table for the entropy mode
      else if ( keymode == NUC_MODE_PRES ) iv = 2; // temperature/energy table for the pressure mode

      vox = (real)0.0;

      pv = alltables_mode + iv + 3*( (ix-1) + (iy-1)*nx + (iz-1)*nxy );

      for (int k=0; k<4; k++)
      {
         q[k] = (real)0.0;

         for (int j=0; j<4; j++)
         {
            r[j] = (real)0.0;

            for (int i=0; i<4; i++)
            {
               r[j] += u[i]* *pv;
               pv   += 3;
            }

            q[k] += v[j]*r[j];
            pv   += 3*nx - 4*3;
         }

         vox += w[k]*q[k];
         pv  += nxy*3 - 4*3*nx;
      } // for (int k=0; k<4; k++)

      *found_ltoreps = vox;


// linear interpolation when cubic interpolations failed
      if ( vox != vox )
      {
//        *keyerr = 20020;
//        return;
//      }
//      else
//      {
//         if ( interpol_TL == NUC_INTERPOL_LINEAR )
//         {
            findtoreps_bdry( x, y, z, found_ltoreps, alltables_mode, nx, ny, nz,
                             ntoreps, xt, yt, zt, logtoreps, keymode, keyerr );
//         }
         return;
      }



   return;


   }


} // FUNCTION : findtoreps



//-------------------------------------------------------------------------------------
// Function    :  findtoreps_bdry
// Description :  Finding temperature from different modes
//                --->                energy   mode    (0)
//                                    temperature mode (1)
//                                    entropy  mode    (2)
//                                    pressure mode    (3)
//
// Note        :  1. Use linear interpolation at boundaries of table to search the
//                   corresponding temperature/energy given (rho, (eps, T, entropy, P), Y_e)
//                2. Invoked by findtoreps()
//
// Parameter   :  x              : Input vector of first  variable (rho)
//                y              : Input vector of second variable (eps, T, entropy, P)
//                z              : Input vector of third  variable (Y_e)
//                found_ltoreps  : Output log(T)/log(energy) of interpolated function values
//                alltables_mode : 3d array of tabulated logenergy
//                nx             : X-dimension of table
//                ny             : Y-dimension of table
//                nz             : Z-dimension of table
//                ntoreps        : Size of temperature/energy array in the Nuclear EoS table
//                xt             : Vector of x-coordinates of table
//                yt             : Vector of y-coordinates of table
//                zt             : Vector of z-coordinates of table
//                logtoreps      : log(temp)/log(energy) array in the table
//                keymode        : Which mode we will use
//                                 0: energy mode      (coming in with eps)
//                                 1: temperature mode (coming in with T)
//                                 2: entropy mode     (coming in with entropy)
//                                 3: pressure mode    (coming in with P)
//                keyerr         : Output error
//
// Return      :  found_ltoreps
//-------------------------------------------------------------------------------------
void findtoreps_bdry( const real x, const real y, const real z,
                      real *found_ltoreps, const real *alltables_mode,
                      const int nx, const int ny, const int nz, const int ntoreps,
                      const real *xt, const real *yt, const real *zt,
                      const real *logtoreps, const int keymode, int *keyerr )
{


// helper variables
  real fh[8], delx, dely, delz, a[8];
  real dx, dy, dz, dxi, dyi, dzi, dxyi, dxzi, dyzi, dxyzi;
  int  ix, iy, iz;


// determine spacing parameters of equidistant (!!!) table
#  if 1
   dx  = ( xt[nx-1] - xt[0] ) / (real)(nx-1);
   dy  = ( yt[ny-1] - yt[0] ) / (real)(ny-1);
   dz  = ( zt[nz-1] - zt[0] ) / (real)(nz-1);

   dxi = (real)1.0 / dx;
   dyi = (real)1.0 / dy;
   dzi = (real)1.0 / dz;
#  endif

#  if 0
   dx = drho;
   dy = dtemp;
   dz = dye;

   dxi = drhoi;
   dyi = dtempi;
   dzi = dyei;
#  endif

   dxyi  = dxi*dyi;
   dxzi  = dxi*dzi;
   dyzi  = dyi*dzi;
   dxyzi = dxi*dyi*dzi;


// determine location in table
   ix = 1 + (int)( ( x - xt[0] )*dxi );
   iy = 1 + (int)( ( y - yt[0] )*dyi );
   iz = 1 + (int)( ( z - zt[0] )*dzi );

   ix = MAX( 1, MIN( ix, nx-1 ) );
   iy = MAX( 1, MIN( iy, ny-1 ) );
   iz = MAX( 1, MIN( iz, nz-1 ) );


// set up aux vars for interpolation
   delx = xt[ix] - x;
   dely = yt[iy] - y;
   delz = zt[iz] - z;

   int idx[8];
   idx[0] = 3*(  (ix  ) + nx*( (iy  ) + ny*(iz  ) )  );
   idx[1] = 3*(  (ix-1) + nx*( (iy  ) + ny*(iz  ) )  );
   idx[2] = 3*(  (ix  ) + nx*( (iy-1) + ny*(iz  ) )  );
   idx[3] = 3*(  (ix  ) + nx*( (iy  ) + ny*(iz-1) )  );
   idx[4] = 3*(  (ix-1) + nx*( (iy-1) + ny*(iz  ) )  );
   idx[5] = 3*(  (ix-1) + nx*( (iy  ) + ny*(iz-1) )  );
   idx[6] = 3*(  (ix  ) + nx*( (iy-1) + ny*(iz-1) )  );
   idx[7] = 3*(  (ix-1) + nx*( (iy-1) + ny*(iz-1) )  );


   int iv;
#  if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
   if      ( keymode == NUC_MODE_ENGY ) iv = 0; // temperature table for the energy mode
#  elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   if      ( keymode == NUC_MODE_TEMP ) iv = 0; // energy table for the temperautre mode
#  endif
   else if ( keymode == NUC_MODE_ENTR ) iv = 1; // temperature/energy table for the entropy mode
   else if ( keymode == NUC_MODE_PRES ) iv = 2; // temperature/energy table for the pressure mode


// set up aux vars for interpolation assuming array ordering (iv, ix, iy, iz)
   fh[0] = alltables_mode[ iv + idx[0] ];
   fh[1] = alltables_mode[ iv + idx[1] ];
   fh[2] = alltables_mode[ iv + idx[2] ];
   fh[3] = alltables_mode[ iv + idx[3] ];
   fh[4] = alltables_mode[ iv + idx[4] ];
   fh[5] = alltables_mode[ iv + idx[5] ];
   fh[6] = alltables_mode[ iv + idx[6] ];
   fh[7] = alltables_mode[ iv + idx[7] ];


// set up coeffs of interpolation polynomical and evaluate function values
   a[0] = fh[0];
   a[1] = dxi  *( fh[1] - fh[0] );
   a[2] = dyi  *( fh[2] - fh[0] );
   a[3] = dzi  *( fh[3] - fh[0] );
   a[4] = dxyi *( fh[4] - fh[1] - fh[2] + fh[0] );
   a[5] = dxzi *( fh[5] - fh[1] - fh[3] + fh[0] );
   a[6] = dyzi *( fh[6] - fh[2] - fh[3] + fh[0] );
   a[7] = dxyzi*( fh[7] - fh[0] + fh[1] + fh[2] +
                  fh[3] - fh[4] - fh[5] - fh[6] );


   *found_ltoreps = a[0]
                  + a[1]*delx
                  + a[2]*dely
                  + a[3]*delz
                  + a[4]*delx*dely
                  + a[5]*delx*delz
                  + a[6]*dely*delz
                  + a[7]*delx*dely*delz;

   if ( *found_ltoreps == NAN ||
        !( *found_ltoreps>=logtoreps[0] && *found_ltoreps<=logtoreps[ntoreps-1] )  ) {
      *keyerr = 665;
   }
   return;

} // FUNCTION : findtoreps_bdry
