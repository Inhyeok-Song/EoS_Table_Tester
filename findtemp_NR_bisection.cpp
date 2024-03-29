#include "NuclearEoS.h"


void nuc_eos_C_linterp_for_temp( const real x, const real y, const real z,
                                 real* f, const real* ft,
                                 const int nx, const int ny, const int nz,
                                 const real* xt, const real* yt, const real* zt,
                                 real* dvardlt, const int keymode );
void bisection( const real lr, const real lt0, const real ye, const real var0, real* ltout,
                const int nrho, const int ntemp, const int nye, const real *alltables,
                const real *logrho, const real *logtemp, const real *yes,
                const int keymode, int* keyerr, const real prec );
real linterp2D( const real* xs, const real* ys, const real* fs, const real x, const real y );



//-------------------------------------------------------------------------------------
// Function    :  findtemp_NR_bisection
// Description :  Finding temperature from different modes
//                ---> energy   mode (0)
//                     entropy  mode (2)
//                     pressure mode (3)
//                by Newton-Raphson and bisection methods
// Note        :  1. Invoked by nuc_eos_C_short()
//
// Parameter   :  lr          : log(rho)
//                lt0         : log(T0) (initial guess of T)
//                ye          : Ye
//                varin       : Input variable (eps, e or P)
//                ltout       : Ouput temperature
//                logrho      : logrho array in the table
//                logtemp     : logtemp array in the table
//                yes         : Ye array in the table
//                keymode     : Which mode we will use
//                              0: energy mode      (coming in with eps)
//                              2: entropy mode     (coming in with entropy)
//                              3: pressure mode    (coming in with P)
//                keyerr      : Output error
//                prec        : tolerance
//
// Return      :  ltout
//-------------------------------------------------------------------------------------
void findtemp_NR_bisection( const real lr, const real lt0, const real ye, const real varin, real *ltout,
                            const int nrho, const int ntemp, const int nye, const real *alltables,
                            const real *logrho, const real *logtemp, const real *yes,
                            const int keymode, int *keyerr, const real prec )
{


// local vars
   int  itmax = 20;               // use at most 20 iterations, then go to bisection
   real dvardlt;                  // derivative dlogprs/dlogT
   real ldt;
   real var, var0, var1;          // temp vars for finding value
   real ltn, lt, lt1;             // temp vars for temperature
   real ltmax = logtemp[ntemp-1]; // max temp
   real ltmin = logtemp[0];       // min temp

   real dtemp, drho, dye;
   real dtempi, drhoi, dyei;


// determine spacing parameters of equidistant (!!!) table
#  if 1
   dtemp = ( logtemp[ntemp-1] - logtemp[0] ) / ( (real)1.0*( ntemp-1 ) );
   drho  = (  logrho[nrho -1] -  logrho[0] ) / ( (real)1.0*( nrho -1 ) );
   dye   = (     yes[nye  -1] -     yes[0] ) / ( (real)1.0*( nye  -1 ) );

   dtempi = (real)1.0 / dtemp;
   drhoi  = (real)1.0 / drho;
   dyei   = (real)1.0 / dye;
#  endif

// setting up some vars
   *keyerr = 0;
   var0    = varin;
   var1    = var0;
   lt      = lt0;
   lt1     = lt;

// step 1: do we already have the right temperature
   nuc_eos_C_linterp_for_temp( lr, lt, ye, &var, alltables, nrho, ntemp, nye,
                               logrho, logtemp, yes, &dvardlt, keymode );

   if( FABS( var-var0 ) < prec*FABS( var0 ) ) {
      *ltout = lt0;
      return;
   }
   lt1  = lt;
   var1 = var;

   int it = 0;
   while( it < itmax ) {

// step 2: check if the two bounding values of the temperature
//         give values that enclose the new values.
      int itemp = MIN( MAX( 1 + (int)( ( lt - logtemp[0] )*dtempi ), 1), ntemp-1 );
      int irho  = MIN( MAX( 1 + (int)( ( lr -  logrho[0] )*drhoi  ), 1), nrho-1  );
      int iye   = MIN( MAX( 1 + (int)( ( ye -     yes[0] )*dyei   ), 1), nye-1   );

      int iv;
      if      ( keymode == NUC_MODE_ENGY ) iv = 1; // energy mode
      else if ( keymode == NUC_MODE_ENTR ) iv = 2; // entropy mode
      else if ( keymode == NUC_MODE_PRES ) iv = 0; // pressure mode

      real vart1, vart2;
// lower vars
      {
// get data at 4 points
         real fs[4];
         // point 0
         int ifs;
         ifs = iv + NUC_TABLE_NVAR*(  (irho-1) + nrho*( (itemp-1) + ntemp*(iye-1) )  );
         fs[0] = alltables[ifs];
         // point 1
         ifs = iv + NUC_TABLE_NVAR*(  (irho  ) + nrho*( (itemp-1) + ntemp*(iye-1) )  );
         fs[1] = alltables[ifs];
         // point 2
         ifs = iv + NUC_TABLE_NVAR*(  (irho-1) + nrho*( (itemp-1) + ntemp*(iye  ) )  );
         fs[2] = alltables[ifs];
         // point 3
         ifs = iv + NUC_TABLE_NVAR*(  (irho  ) + nrho*( (itemp-1) + ntemp*(iye  ) )  );
         fs[3] = alltables[ifs];

         vart1 = linterp2D( &logrho[irho-1], &yes[iye-1], fs, lr, ye );
      }
// upper vars
      {
// get data at 4 points
         real fs[4];
         // point 0
         int ifs;
         ifs = iv + NUC_TABLE_NVAR*(  (irho-1) + nrho*( (itemp) + ntemp*(iye-1) )  );
         fs[0] = alltables[ifs];
         // point 1
         ifs = iv + NUC_TABLE_NVAR*(  (irho  ) + nrho*( (itemp) + ntemp*(iye-1) )  );
         fs[1] = alltables[ifs];
         // point 2
         ifs = iv + NUC_TABLE_NVAR*(  (irho-1) + nrho*( (itemp) + ntemp*(iye  ) )  );
         fs[2] = alltables[ifs];
         // point 3
         ifs = iv + NUC_TABLE_NVAR*(  (irho  ) + nrho*( (itemp) + ntemp*(iye  ) )  );
         fs[3] = alltables[ifs];

         vart2 = linterp2D( &logrho[irho-1], &yes[iye-1], fs, lr, ye );
      }


// Check if we are already bracketing the input internal
// variable. If so, interpolate for new T.
      if ( var0 >= vart1 && var0 <= vart2 ) {
         *ltout = ( logtemp[itemp] - logtemp[itemp-1] )/( vart2 - vart1 )*
                  ( var0 - vart1 ) + logtemp[itemp-1];

         return;

      }



// well, then do a Newton-Raphson step
   ldt  = -( var - var0 )/dvardlt;
   ltn  = MIN( MAX( lt + ldt, ltmin ), ltmax );
   lt1  = lt;
   lt   = ltn;
   var1 = var;

   nuc_eos_C_linterp_for_temp( lr, lt, ye, &var, alltables, nrho, ntemp, nye,
                               logrho, logtemp, yes, &dvardlt, keymode );

   if( FABS( var - var0 ) < prec*FABS( var0 ) ) {
      *ltout = lt;
      return;
   }

// if we are closer than 10^-3  to the
// root (prs-prs0)=0, we are switching to
// the secant method, since the table is rather coarse and the
// derivatives may be garbage.
      if ( FABS( var - var0 ) < 1.0e-3*FABS( var0 ) ) {
         dvardlt = ( var - var1 )/( lt - lt1 );
      }

      it++;
   }

   if ( it >= itmax-1 ) {
// try bisection
      bisection( lr, lt0, ye, var0, ltout, nrho, ntemp, nye,
                 alltables, logrho, logtemp, yes, keymode, keyerr, prec );
      return;
   }


   return;

} // FUNCTION : findtemp_NR_bisection




//-------------------------------------------------------------------------------------
// Function    :  linterp2D
// Description :  2D interpolation for tabulated variables
//
// Note        :  1. Invoked by findtemp_NR_bisection()
//
// Parameter   :  xs   : X coordinates
//                ys   : Y coordinates
//                fs   : Function values
//                 x   : X position
//                 y   : Y position
// Return      :  Interpolated value
//-------------------------------------------------------------------------------------
real linterp2D( const real* xs, const real* ys, const real* fs, const real x, const real y )
{


//  2     3
//
//  0     1
//
// first interpolate in x between 0 and 1, 2 and 3
// then interpolate in y
// assume rectangular grid

   real t1 = ( fs[1] - fs[0] )/( xs[1] - xs[0] )*( x - xs[0] ) + fs[0];
   real t2 = ( fs[3] - fs[2] )/( xs[1] - xs[0] )*( x - xs[0] ) + fs[2];

   return ( t2 - t1 )/( ys[1] - ys[0] )*( y - ys[0] ) + t1;
} // FUNCTION : linterp2D




//-------------------------------------------------------------------------------------
// Function    :  nuc_eos_C_linterp_for_temp
// Description :  Interpolation of a function of three variables in an
//                equidistant(!!!) table.
//                method:  8-point Lagrange linear interpolation formula
//
// Note        :  1. Invoked by findtemp_NR_bisection()
//
// Parameter   :  x       : Input vector of first  variable
//                y       : Input vector of second variable
//                z       : Input vector of third  variable
//                f       : Output vector of interpolated function values
//                ft      : 3d array of tabulated function values
//                nx      : X-dimension of table
//                ny      : Y-dimension of table
//                nz      : Z-dimension of table
//                xt      : Vector of x-coordinates of table
//                yt      : Vector of y-coordinates of table
//                zt      : Vector of z-coordinates of table
//                dvardlt : dvar/dlog(T)
//                keymode : which mode we will use
// Return      :  f
//-------------------------------------------------------------------------------------
void nuc_eos_C_linterp_for_temp( const real x, const real y, const real z,
                                 real* f, const real* ft,
                                 const int nx, const int ny, const int nz,
                                 const real* xt, const real* yt, const real* zt,
                                 real* dvardlt, const int keymode )
{


// helper variables
   real fh[8], delx, dely, delz, a[8];
   real dx, dy, dz, dxi, dyi, dzi, dxyi, dxzi, dyzi, dxyzi;
   int  ix, iy, iz;

// determine spacing parameters of equidistant (!!!) table
   dx = ( xt[nx-1] - xt[0] )/( (real)1.0*(nx-1) );
   dy = ( yt[ny-1] - yt[0] )/( (real)1.0*(ny-1) );
   dz = ( zt[nz-1] - zt[0] )/( (real)1.0*(nz-1) );

   dxi = 1.0/dx;
   dyi = 1.0/dy;
   dzi = 1.0/dz;

   dxyi = dxi*dyi;
   dxzi = dxi*dzi;
   dyzi = dyi*dzi;

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

   idx[0] = NUC_TABLE_NVAR*(  (ix  ) + nx*( (iy  ) + ny*(iz  ) )  );
   idx[1] = NUC_TABLE_NVAR*(  (ix-1) + nx*( (iy  ) + ny*(iz  ) )  );
   idx[2] = NUC_TABLE_NVAR*(  (ix  ) + nx*( (iy-1) + ny*(iz  ) )  );
   idx[3] = NUC_TABLE_NVAR*(  (ix  ) + nx*( (iy  ) + ny*(iz-1) )  );
   idx[4] = NUC_TABLE_NVAR*(  (ix-1) + nx*( (iy-1) + ny*(iz  ) )  );
   idx[5] = NUC_TABLE_NVAR*(  (ix-1) + nx*( (iy  ) + ny*(iz-1) )  );
   idx[6] = NUC_TABLE_NVAR*(  (ix  ) + nx*( (iy-1) + ny*(iz-1) )  );
   idx[7] = NUC_TABLE_NVAR*(  (ix-1) + nx*( (iy-1) + ny*(iz-1) )  );

   int iv;
   if ( keymode == NUC_MODE_ENGY ) iv = 1; // energy mode
   if ( keymode == NUC_MODE_ENTR ) iv = 2; // entropy mode
   if ( keymode == NUC_MODE_PRES ) iv = 0; // pressure mode

// set up aux vars for interpolation
// assuming array ordering (iv, ix, iy, iz)
   fh[0] = ft[iv+idx[0]];
   fh[1] = ft[iv+idx[1]];
   fh[2] = ft[iv+idx[2]];
   fh[3] = ft[iv+idx[3]];
   fh[4] = ft[iv+idx[4]];
   fh[5] = ft[iv+idx[5]];
   fh[6] = ft[iv+idx[6]];
   fh[7] = ft[iv+idx[7]];

// set up coeffs of interpolation polynomical and
// evaluate function values
   a[0] = fh[0];
   a[1] = dxi  * ( fh[1] - fh[0] );
   a[2] = dyi  * ( fh[2] - fh[0] );
   a[3] = dzi  * ( fh[3] - fh[0] );
   a[4] = dxyi * ( fh[4] - fh[1] - fh[2] + fh[0] );
   a[5] = dxzi * ( fh[5] - fh[1] - fh[3] + fh[0] );
   a[6] = dyzi * ( fh[6] - fh[2] - fh[3] + fh[0] );
   a[7] = dxyzi* ( fh[7] - fh[0] + fh[1] + fh[2] +
                   fh[3] - fh[4] - fh[5] - fh[6] );

   *dvardlt = -a[2];

   *f = a[0]
      + a[1] * delx
      + a[2] * dely
      + a[3] * delz
      + a[4] * delx * dely
      + a[5] * delx * delz
      + a[6] * dely * delz
      + a[7] * delx * dely * delz;


   return;

} // FUNCTION : nuc_eos_C_linterp_for_temp



//-------------------------------------------------------------------------------------
// Function    :  bisection
// Description :  Finding temperature from different modes
//                --->                energy   mode (0)
//                                    entropy  mode (2)
//                                    pressure mode (3)
//                by a bisection method
// Note        :  1. Invoked by findtemp_NR_bisection when the Newton-Raphson method failed.
//
// Parameter   :  lr          : log(rho)
//                lt0         : log(T0) (initial guess of T)
//                ye          : Ye
//                var0        : Input variable (eps, e or P)
//                ltout       : Ouput temperature
//                logrho      : logrho array in the table
//                logtemp     : logtemp array in the table
//                yes         : Ye array in the table
//                keymode     : Which mode we will use
//                              --> 1: Temperature mode (coming in with T)
//                                  2: Entropy mode     (coming in with entropy)
//                                  3: Pressure mode    (coming in with P)
//                keyerr      : Output error
//
// Return      :  ltout
//-------------------------------------------------------------------------------------
void bisection( const real lr, const real lt0, const real ye, const real var0,
                real* ltout, const int nrho, const int ntemp, const int nye,
                const real *alltables, const real *logrho, const real *logtemp, const real *yes,
                const int keymode, int* keyerr, const real prec )
{


// iv is the index of the table variable we do the bisection on
   int bcount    = 0;
   int maxbcount = 30;
   int itmax     = 50;

// temporary local vars
   real lt, lt1, lt2;
   real ltmin = logtemp[0];
   real ltmax = logtemp[ntemp-1];
   real f1, f2, fmid, dlt, ltmid;
   real dvardlt;
   real f1a[3] = {0.0};
   real f2a[3] = {0.0};

   int iv;

   if      ( keymode == NUC_MODE_ENGY ) iv = 1;
   else if ( keymode == NUC_MODE_ENTR ) iv = 2;
   else if ( keymode == NUC_MODE_PRES ) iv = 0;


// prepare
   lt = lt0;
   lt1 = LOG10( MIN ( POW( (real)10.0, ltmax ), ( 1.2 ) * ( POW( (real)10.0, lt0 ) ) ) );
   lt2 = LOG10( MAX ( POW( (real)10.0, ltmin ), ( 0.8 ) * ( POW( (real)10.0, lt0 ) ) ) );

   int nvars = 3;
   nuc_eos_C_linterp_some( lr, lt1, ye, f1a, alltables,
                           nrho, ntemp, nye, nvars, logrho, logtemp, yes );

   nuc_eos_C_linterp_some( lr, lt2, ye, f2a, alltables,
                           nrho, ntemp, nye, nvars, logrho, logtemp, yes );

   f1 = f1a[iv] - var0;
   f2 = f2a[iv] - var0;

// iterate until we bracket the right eps, but enforce
// dE/dt > 0, so eps(lt1) > eps(lt2)
#  if 0
   int ifixdeg = 0;
   int ifixdeg_max = 20;
#  endif

   while( f1*f2 >= 0.0 ) {

      lt1 = LOG10( MIN( POW( (real)10.0, ltmax ), ( 1.2 )*( POW( (real)10.0, lt1 ) ) ) );
      lt2 = LOG10( MAX( POW( (real)10.0, ltmin ), ( 0.8 )*( POW( (real)10.0, lt2 ) ) ) );
      nuc_eos_C_linterp_some( lr, lt1, ye, f1a, alltables,
                              nrho, ntemp, nye, nvars, logrho, logtemp, yes );
      nuc_eos_C_linterp_some( lr, lt2, ye, f2a, alltables,
                              nrho, ntemp, nye, nvars, logrho, logtemp, yes );

#  if 0
   // special enforcement of eps(lt1)>eps(lt2)
   while( f1a < f2a && ifixdeg < ifixdeg_max ) {
      lt1 = LOG10( MIN ( POW ( (real)10.0, ltmax ), 1.2*( POW( (real)10.0, lt1 ) ) ) );
      nuc_eos_C_linterp_some( lr, lt1, ye, f1a, alltables,
                              nrho, ntemp, nye, nvars, logrho, logtemp, yes );
      ifixdeg++;
   }
#  endif

      f1 = f1a[iv] - var0;
      f2 = f2a[iv] - var0;

      bcount++;
      if( bcount >= maxbcount ) {
         *keyerr = 668;
         return;
      }

   }

   if( f1 < 0.0 ) {
      lt = lt1;
      dlt = LOG10( POW( (real)10.0, lt2 ) - POW( (real)10.0, lt1 ) );
   } else {
      lt = lt2;
      dlt = LOG10( POW( (real)10.0, lt1 ) - POW( (real)10.0, lt2 ) );
   }

   int it;
   for ( it=0; it<itmax; it++ ) {
      dlt   = LOG10( POW( (real)10.0, dlt )*(real)0.5 );
      ltmid = LOG10( POW( (real)10.0, lt ) + POW( (real)10.0, dlt ) );
      nuc_eos_C_linterp_some( lr, ltmid, ye, f2a, alltables,
                                nrho, ntemp, nye, nvars, logrho, logtemp, yes );
      fmid = f2a[iv] - var0;
      if ( fmid <= 0.0 ) lt=ltmid;

      if( FABS( (real)1.0 - f2a[iv]/var0 ) <= (real)prec ) {
         *ltout = ltmid;
         return;
      }

   }

   if( it >= itmax-1 ) {
      *keyerr = 669;
      return;
   }


} // FUNCTION : bisection
