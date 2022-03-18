#include "NuclearEoS.h"



//-----------------------------------------------------------------------------------------------
// Function    :  nuc_eos_C_short
// Description :  Function to find thermodynamic varibles by searching
//                a pre-calculated nuclear equation of state table
//
// Note        :  1. It will strictly return values in cgs or MeV
//                2. Four modes are supported
//                3. The defalut mode is temperature (1) mode
//                4. In case, three other modes are available for finding temperature
//                   energy      (0) mode
//                   temperature (1) mode
//                   entropy     (2) mode
//                   pressure    (3) mode
//
// Parameter   :  xrho            : Input density (rho (g/cm^3))
//                xtemp           : Input (temperature mode)
//                                  or ouput temperature in MeV
//                xye             : Electron fraction (Y_e)
//                xenr            : Input specific internal energy (energy mode)
//                                  or output specific internal energy
//                xent            : Input (entropy mode)
//                                  or output specific entropy (e)
//                xprs            : Input (pressure mode)
//                                  or output pressure
//                xcs2            : Output sound speed
//                xmunu           : Output chemcial potential
//                energy_shift    : energy_shift
//                nrho            : Size of density array in the Nuclear EoS table
//                ntoreps         : Size of (temperature/energy) array in the Nuclear EoS table (temp/energy-based table)
//                nye             : Size of Y_e array in the Nuclear EoS table
//                nrho_mode       : Size of density array of look-up tables (logtemp_energy.../logenergy_temp...)
//                nmode           : Size of log(eps/temp) arrays of look-up tables
//                                          entropy
//                                          log(P)        array in the Nuclear EoS table
//                                                        for each mode
//                nye_mode        : Size of Y_e array of lookup tables
//                alltables       : Nuclear EoS table
//                alltables_mode  : Auxiliary log(T/eps) arrays for energy/temperature mode
//                                                                  entropy mode
//                                                                  pressure mode
//                logrho          : log(rho) array in the table
//                logtoreps       : log(T) or log(eps) array for (T/eps) mode (temp/energy-based table)
//                yes             : Y_e      array in the table
//                logepsort_mode  : log(eps) or log(T) array for (eps/T) mode (temp/energy-based table)
//                entropy_mode    : entropy  array for entropy mode
//                logpress_mode   : log(P)   array for pressure mode
//                interpol_TL     : interpolation schemes for table look-ups (linear/cubic)
//                interpol_other  : interpolation schemes for other thermodynamic variables (linear/cubic)
//                keymode         : Which mode we will use
//                                  0 : energy mode      (coming in with eps)
//                                  1 : temperature mode (coming in with T)
//                                  2 : entropy mode     (coming in with entropy)
//                                  3 : pressure mode    (coming in with P)
//                keyerr          : Output error
//                                  669 : bisection failed (only in tempearture-based table)
//                                  668 : bisection failed (only in tempearture-based table)
//                                  665 : fail in finding T/e
//                                  101 : Y_e too high
//                                  102 : Y_e too low
//                                  103 : temp too high
//                                  104 : temp too low
//                                  105 : rho too high
//                                  106 : rho too low
//                                  107 : eps too high     (energy-based table)
//                                  108 : eps too low      (energy-based table)
//                                  109 : temp too high    (energy-based table)
//                                  110 : temp too low     (energy-based table)
//                                  111 : entropy too high (energy-based table)
//                                  112 : entropy too low  (energy-based table)
//                                  113 : log(P) too high  (energy-based table)
//                                  114 : log(P) too low   (energy-based table)
//                                  201 : lr  has NaN value
//                                  202 : xye has NaN value
//                                  203 : xenr has NaN value
//                                  204 : xtemp has NaN value
//                                  205 : xent has NaN value
//                                  206 : xprs has NaN value

//                rfeps           : Tolerence for interpolations
//-----------------------------------------------------------------------------------------------
void nuc_eos_C_short( const real xrho, real *xtemp, const real xye,
                      real *xenr, real *xent, real *xprs,
                      real *xcs2, real *xmunu, const real energy_shift,
                      const int nrho, const int ntoreps, const int nye,
                      const int nrho_mode, const int nmode, const int nye_mode,
                      const real *alltables, const real *alltables_mode,
                      const real *logrho, const real *logtoreps, const real *yes,
                      const real *logrho_mode, const real *logepsort_mode,
                      const real *entr_mode, const real *logprss_mode, const real *yes_mode,
                      const int interpol_TL, const int interpol_other,
                      const int keymode, int *keyerr, const real rfeps )
{

// check whether the input density and Ye are within the table
   const real lr = LOG10( xrho );
   *keyerr = 0;


   if ( lr >  logrho[nrho-1] )  {  *keyerr = 105;  return;  }
   if ( lr <  logrho[     0] )  {  *keyerr = 106;  return;  }
   if ( lr != lr             )  {  *keyerr = 201;  return;  }

   if ( xye >  yes  [nye -1] )  {  *keyerr = 101;  return;  }
   if ( xye <  yes  [     0] )  {  *keyerr = 102;  return;  }
   if ( xye != xye           )  {  *keyerr = 202;  return;  }

// find temperature (temperature-based table)
#  if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
   real ltoreps   = LOG10( *xtemp );
   real lt0       = NULL_REAL;
   if ( *xtemp == *xtemp )
   {
      lt0 = ltoreps;
      lt0 = MAX(  MIN( lt0, logtoreps[ntoreps-1] ), logtoreps[0]  );
   }
   else if ( *xtemp != *xtemp )
   {
      lt0 = 1.0;
   }
#  elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   real ltoreps   = LOG10( *xenr + energy_shift );
#  endif // #elif NUC_TABLE_MODE ... else ...
   real var0            = NULL_REAL;
   const real *mode_arr = NULL;

   switch ( keymode )
   {
      case NUC_MODE_ENGY :
      {
#     if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
         const real leps = LOG10( MAX( *xenr + energy_shift, 1.0 ) );
         var0            = leps;
         mode_arr        = logepsort_mode;
         if ( leps > logepsort_mode[nmode-1]  )     {  *keyerr = 107;           }
         if ( leps < logepsort_mode[      0]  )     {  *keyerr = 108;           }
         if ( *xenr != *xenr                  )     {  *keyerr = 203;           }
#     elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
         real leps = LOG10( MAX( *xenr + energy_shift, 1.0 ) );
         if ( leps > logtoreps[ntoreps-1] )         {  *keyerr = 107;  return;  }
         if ( leps < logtoreps[        0] )         {  *keyerr = 108;  return;  }
         if ( *xenr != *xenr              )         {  *keyerr = 203;           }
#     endif // #elif NUC_TABLE_MODE ... else ...
      }
      break;

      case NUC_MODE_TEMP :
      {
#     if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
         real lt  = LOG10( *xtemp );
         if ( lt > logtoreps[ntoreps-1] )           {  *keyerr = 103;  return;  }
         if ( lt < logtoreps[        0] )           {  *keyerr = 104;  return;  }
         if ( *xtemp != *xtemp          )           {  *keyerr = 204;  return;  }
#     elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
         const real lt  = LOG10( *xtemp );
         mode_arr = logepsort_mode;
         var0     = lt;
         if ( lt > logepsort_mode[nmode-1] )        {  *keyerr = 109;  return;  }
         if ( lt < logepsort_mode[      0] )        {  *keyerr = 110;  return;  }
         if ( *xtemp != *xtemp             )        {  *keyerr = 204;  return;  }
#     endif // #elif NUC_TABLE_MODE ... else ...
      }
      break;

      case NUC_MODE_ENTR :
      {
         const real entr = *xent;
         var0            =  entr;
         mode_arr        =  entr_mode;
         if ( entr > entr_mode[nmode-1] )           {  *keyerr = 111;           }
         if ( entr < entr_mode[      0] )           {  *keyerr = 112;           }
         if ( *xent != *xent            )           {  *keyerr = 205;  return;  }
      }
      break;

      case NUC_MODE_PRES :
      {
         const real lprs = LOG10( *xprs );
         var0            = lprs;
         mode_arr        = logprss_mode;
         if ( lprs > logprss_mode[nmode-1] )        {  *keyerr = 113;           }
         if ( lprs < logprss_mode[      0] )        {  *keyerr = 114;           }
         if ( *xprs != *xprs               )        {  *keyerr = 206;  return;  }
      }
      break;
   } // switch ( keymode )


#  if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
   if ( keymode != NUC_MODE_TEMP ) {
      if ( *keyerr == 0 ) {
//    find temperature from energy, entorpy or pressure
         findtoreps( lr, var0, xye, &ltoreps, alltables_mode, nrho_mode, nmode, nye_mode, ntoreps,
                     logrho_mode, mode_arr, yes_mode, logtoreps, interpol_TL, keymode, keyerr );
      }
      if ( *keyerr != 0 ) {
//       fprintf( stdout, "keyerr: %d\n", *keyerr );
//       find temperature by Netwon-Rapshon or bisection if above failed
         *keyerr = 901; return;
         // findtemp_NR_bisection( lr, lt0, xye, var0, &ltoreps, nrho, ntoreps, nye, alltables,
         //                        logrho, logtoreps, yes, keymode, keyerr, rfeps );
      }
   }


#  elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   if ( *keyerr != 0 ) return;
   if ( keymode != NUC_MODE_ENGY )
   {
// find energy from temperature, entorpy or pressure
      findtoreps( lr, var0, xye, &ltoreps, alltables_mode, nrho_mode, nmode, nye_mode, ntoreps,
                  logrho_mode, mode_arr, yes_mode, logtoreps, interpol_TL, keymode, keyerr );
      if ( *keyerr != 0 ) return;
   }
#  endif // #elif NUC_TABLE_MODE ... else ...


   real res[5]; // result array

// linear interolation for other variables
   if      ( interpol_other == NUC_INTERPOL_LINEAR )
   {
      nuc_eos_C_linterp_some( lr, ltoreps, xye, res, alltables,
                              nrho, ntoreps, nye, 5, logrho, logtoreps, yes );
   }
// cubic interpolation for other variables
   else if ( interpol_other == NUC_INTERPOL_CUBIC  )
   {
      nuc_eos_C_cubinterp_some( lr, ltoreps, xye, res, alltables,
                                nrho, ntoreps, nye, 5, logrho, logtoreps, yes );
   }

// assign results
#  if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
   if ( keymode != NUC_MODE_TEMP ) *xtemp = POW( (real)10.0, ltoreps );
   if ( keymode != NUC_MODE_ENGY ) *xenr  = POW( (real)10.0, res[1]  ) - energy_shift;
#  elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   if ( keymode != NUC_MODE_ENGY ) *xenr  = POW( (real)10.0, ltoreps ) - energy_shift;
   if ( keymode != NUC_MODE_TEMP ) *xtemp = POW( (real)10.0, res[1] );
#  endif // #elif NUC_TABLE_MODE ... else ...
   if ( keymode != NUC_MODE_ENTR ) *xent  = res[2];
   if ( keymode != NUC_MODE_PRES ) *xprs  = POW( (real)10.0, res[0] );
   *xmunu = res[3];
   *xcs2  = res[4];


   return;

} // FUNCTION : nuc_eos_C_short
