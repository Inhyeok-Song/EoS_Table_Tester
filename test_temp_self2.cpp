#define H5_USE_16_API 1
#include "hdf5.h"
#include "NuclearEoS.h"

#define MAXCHAR 1000
#define N_MAX   200
#define DATASET "integer scalars"
#define DIM0 15


static int file_is_readable( const char *filename );
static int file_is_readable( const char *filename )
{
   FILE* fp = NULL;
      fp = fopen( filename, "r" );
      if( fp != NULL )
      {
         fclose( fp );
         return 1;
      }
   return 0;
}



extern real  *h_EoS_Table[EOS_NTABLE_MAX];

extern int    g_nrho;
extern int    g_nye;
extern int    g_nrho_mode;
extern int    g_nmode;
extern int    g_nye_mode;
extern double g_energy_shift;

extern real  *g_alltables;
extern real  *g_alltables_mode;
extern real  *g_logrho;
extern real  *g_yes;
extern real  *g_logrho_mode;
extern real  *g_entr_mode;
extern real  *g_logprss_mode;
extern real  *g_yes_mode;

extern int    g_ntemp;
extern real  *g_logtemp;
extern real  *g_logeps_mode;

extern int interpol_TL;
extern int interpol_other;


#if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
//-------------------------------------------------------
// Function    :  nuc_eos_C_testing
// Description :
//-------------------------------------------------------
void nuc_eos_C_testing_temp_self_consistency_whole_domain()
{


   // output files for writing errors
   FILE *output_error[6];

   FILE *ener_err1 = fopen("TEMP_ener_err1.txt", "w");
   FILE *ener_err2 = fopen("TEMP_ener_err2.txt", "w");

   FILE *entr_err1 = fopen("TEMP_entr_err1.txt", "w");
   FILE *entr_err2 = fopen("TEMP_entr_err2.txt", "w");

   FILE *prss_err1 = fopen("TEMP_prss_err1.txt", "w");
   FILE *prss_err2 = fopen("TEMP_prss_err2.txt", "w");

   output_error[0] = ener_err1;
   output_error[1] = entr_err1;
   output_error[2] = prss_err1;
   output_error[3] = ener_err2;
   output_error[4] = entr_err2;
   output_error[5] = prss_err2;

   for ( int i=0; i<6; i++ )
   {
      if ( i<3 )
      {
         fprintf( output_error[i], "#    Num            logrho           logtemp0         logenergy          logpress           entropy                Ye\n", "w" );
      }
      else if ( i>=3 )
      {
         fprintf( output_error[i], "#    Num            logrho           logtemp0         logtemp1           logenergy          logpress           entropy                Ye         rel_error\n", "w" );
      }
   }

   // variables for writing errors
   int  failcount1[4] = {0};   // count error1: fail in finding temperature
   int  failcount2[4] = {0};   // count error2: relative error > 10^-10
   real max_err[4]    = {0.0};
   real rel_err_array[3]  = {0.0};
   int  rel_err_number[3] = {0};


#ifdef FLOAT8
   const real Tolerance = 1e-10;
#else
   const real Tolerance = 1e-6;
#endif


   const int nrho_diff    = 1*g_nrho_mode;// - 1;
   const int nmode_diff   = 1*g_nmode    ;// - 1;
   const int nye_diff     = g_nye_mode;
   const int total_points = nrho_diff*nmode_diff*nye_diff;

   real logrho_diff[nrho_diff];
   real logenergy_diff[nmode_diff];
   real logpress_diff[nmode_diff];
   real entropy_diff[nmode_diff];
   real ye_diff[nmode_diff];

   for ( int i=0; i<nrho_diff; i++ ) {
      logrho_diff[i] = g_logrho_mode[0] + i*( g_logrho_mode[g_nrho_mode-1] - g_logrho_mode[0] )/( nrho_diff - 1 );
   }
   for ( int j=0; j<nmode_diff; j++ ) {
      logenergy_diff[j] = g_logeps_mode[0]  + j*( g_logeps_mode[g_nmode-1]  - g_logeps_mode[0]  )/( nmode_diff -1 );
      logpress_diff[j]  = g_logprss_mode[0] + j*( g_logprss_mode[g_nmode-1] - g_logprss_mode[0] )/( nmode_diff -1 );
      entropy_diff[j]   = g_entr_mode[0]    + j*( g_entr_mode[g_nmode-1]    - g_entr_mode[0]    )/( nmode_diff -1 );
   }
   for ( int k=0; k<nye_diff; k++ ) {
      ye_diff[k] = g_yes_mode[0] + k*( g_yes_mode[g_nye_mode-1] - g_yes_mode[0] )/( g_nye_mode - 1 );
   }


   for ( int i_mode=0; i_mode<3; i_mode++ ) {
      int keymode;
      if      ( i_mode == 0 ) keymode = NUC_MODE_ENGY; // energy  mode
      else if ( i_mode == 1 ) keymode = NUC_MODE_ENTR; // entropy mode
      else if ( i_mode == 2 ) keymode = NUC_MODE_PRES; // pressure mode

      for ( int i=0; i<nrho_diff; i++ ) {
         for ( int j=0; j<nmode_diff; j++ ) { //integer_scalars[0]
            for ( int k=0; k<nye_diff; k++ ) {
               const real dens  = POW( 10.0, logrho_diff[i] );
                     real eint0 = POW( 10.0, logenergy_diff[j]) - g_energy_shift;
                     real entr0 = entropy_diff[j];
                     real pres0 = POW( 10.0, logpress_diff[j]);
               const real ye    = ye_diff[k];
                     real temp  = NULL_REAL;

                     real cs2   = 0.0;
                     real munu  = 0.0;
                     real eint  = eint0;
                     real pres  = pres0;
                     real entr  = entr0;

               int  keyerr1 = 0;
               int  keyerr2 = 0;
               nuc_eos_C_short( dens, &temp, ye, &eint, &entr, &pres, &cs2, &munu, g_energy_shift,
                                g_nrho, g_ntemp, g_nye, g_nrho_mode, g_nmode, g_nye_mode, g_alltables, g_alltables_mode,
                                g_logrho, g_logtemp, g_yes, g_logrho_mode,
                                g_logeps_mode, g_entr_mode, g_logprss_mode, g_yes_mode,
                                interpol_TL, interpol_other, keymode, &keyerr1, Tolerance );
               if ( keyerr1 == 0 ) {
                  real temp0 = temp;
                  real temp1 = temp;
                  real eint1 = NULL_REAL;
                  real entr1 = NULL_REAL;
                  real pres1 = NULL_REAL;

                  nuc_eos_C_short( dens, &temp1, ye, &eint1, &entr1, &pres1, &cs2, &munu, g_energy_shift,
                                   g_nrho, g_ntemp, g_nye, g_nrho_mode, g_nmode, g_nye_mode, g_alltables, g_alltables_mode,
                                   g_logrho, g_logtemp, g_yes, g_logrho_mode,
                                   g_logeps_mode, g_entr_mode, g_logprss_mode, g_yes_mode,
                                   interpol_TL, interpol_other, NUC_MODE_TEMP, &keyerr2, Tolerance );

                  real rel_err = 0.0;
                  switch ( keymode )
                  {
                     case NUC_MODE_ENGY: rel_err = ( eint1-eint0 )/( eint0+g_energy_shift ); break;
                     case NUC_MODE_ENTR: rel_err = ( entr1-entr0 )/entr0;                    break;
                     case NUC_MODE_PRES: rel_err = ( pres1-pres0 )/pres0;                    break;
                  }
                  if ( keyerr2 != 0 ) {
                     // error 1: fail in finding energy from the [i_mode]
                     failcount1[i_mode]++;
                     //fprintf( stdout, "error code: %d\n", keyerr1);
                     fprintf( output_error[i_mode], "%8d    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e\n",
                              failcount1[i_mode], LOG10( dens ), LOG10( temp0 ), LOG10( eint0 + g_energy_shift ), LOG10( pres0 ), entr0, ye );
                  }
                  else if ( rel_err > 1.e-10 ) {
                     // error 2: relative errors between the FLASH1D Data and [i_mode] > 0.1%
                     failcount2[i_mode]++;
                     fprintf( output_error[i_mode+3], "%8d    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e\n",
                              failcount2[i_mode], LOG10( dens ), LOG10( temp0 ), LOG10( temp1 ), LOG10( eint0 + g_energy_shift ), entr0, LOG10( pres0 ), ye, rel_err );
                     rel_err_array[i_mode] = rel_err_array[i_mode] + rel_err;
                     rel_err_number[i_mode] = rel_err_number[i_mode]+1;
                     max_err[i_mode] = ( rel_err > max_err[i_mode] ) ? rel_err: max_err[i_mode];
                  }
               }
            }
         }
      }
   }

   fprintf( stdout, "Maximum relative error in energy mode       :  %10.8E\n", max_err[0] );
   fprintf( stdout, "Mean    relative error in energy mode       :  %10.8E\n", rel_err_array[0]/rel_err_number[0] );
   fprintf( stdout, "-----------------------------------------------------------------------------\n" );
   fprintf( stdout, "Maximum relative error is entropy mode      :  %10.8E\n", max_err[1] );
   fprintf( stdout, "Mean    relative error is entropy mode      :  %10.8E\n", rel_err_array[1]/rel_err_number[1] );
   fprintf( stdout, "-----------------------------------------------------------------------------\n" );
   fprintf( stdout, "Maximum relative error is pressure mode     :  %10.8E\n", max_err[2] );
   fprintf( stdout, "Mean    relative error is pressure mode     :  %10.8E\n", rel_err_array[2]/rel_err_number[2] );
   fprintf( stdout, "-----------------------------------------------------------------------------\n" );

   fprintf( stdout, "fail in energy   mode %6d points and %8.6E %% in total\n", failcount1[0], 1.0e2*failcount1[0]/total_points );
   fprintf( stdout, "fail in entropy  mode %6d points and %8.6E %% in total\n", failcount1[1], 1.0e2*failcount1[1]/total_points );
   fprintf( stdout, "fail in pressure mode %6d points and %8.6E %% in total\n", failcount1[2], 1.0e2*failcount1[2]/total_points );

   fprintf( stdout, "err>1e-10 in energy   mode %6d points and %8.6E %% in total\n", failcount2[0], 1.0e2*failcount2[0]/total_points );
   fprintf( stdout, "err>1e-10 in entropy  mode %6d points and %8.6E %% in total\n", failcount2[1], 1.0e2*failcount2[1]/total_points );
   fprintf( stdout, "err>1e-10 in pressure mode %6d points and %8.6E %% in total\n", failcount2[2], 1.0e2*failcount2[2]/total_points );


   return;

} // FUNCTION : nuc_eos_C_testing

#endif  // NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP
