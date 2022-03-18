#define H5_USE_16_API 1
#include "hdf5.h"
#include "NuclearEoS.h"

#define MAXCHAR 1000
#define N_MAX   200
#define DATASET "integer scalars"
#define DIM0 15

// Catch HDF5 errors
#define HDF5_ERROR(fn_call)                                \
   do {                                                    \
      int _error_code = fn_call;                           \
      if ( _error_code < 0 )                               \
      {                                                    \
         fprintf( stderr,                                  \
                  "HDF5 call '%s' returned error code %d", \
                   #fn_call, _error_code );                \
         abort();                                          \
      }                                                    \
   } while (0)

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
void nuc_eos_C_testing_temp_self_consistency_2()
{

   // test FLASH data with new eps energy-based driver
   const double KelvinToMeV = 8.61732814974493E-11;


   // input files and output files for writing errors
   FILE *input[3];

   FILE *ener_in   = fopen( "TEMP_ener_err2.txt", "r" );
   FILE *entr_in   = fopen( "TEMP_entr_err2.txt", "r" );
   FILE *prss_in   = fopen( "TEMP_prss_err2.txt", "r" );

   input[0] = ener_in;
   input[1] = entr_in;
   input[2] = prss_in;

   char header[255];
   for ( int i=0; i<3; i++ ) fgets( header, sizeof(header), input[i] );

   FILE *output_error[6];

   FILE *ener_err1 = fopen( "TEMP_ener_err1_2.txt", "w" );
   FILE *ener_err2 = fopen( "TEMP_ener_err2_2.txt", "w" );

   FILE *entr_err1 = fopen( "TEMP_entr_err1_2.txt", "w" );
   FILE *entr_err2 = fopen( "TEMP_entr_err2_2.txt", "w" );

   FILE *pres_err1 = fopen( "TEMP_prss_err1_2.txt", "w" );
   FILE *pres_err2 = fopen( "TEMP_prss_err2_2.txt", "w" );

   output_error[0] = ener_err1;
   output_error[1] = entr_err1;
   output_error[2] = pres_err1;
   output_error[3] = ener_err2;
   output_error[4] = entr_err2;
   output_error[5] = pres_err2;



   for ( int i=0; i<6; i++ )
   {
      if ( i<3 )
      {
         fprintf( output_error[i], "#    Num            logrho           logtemp          logenergy          entropy            logpress               Ye\n", "w" );
      }
      else if ( i>=3 )
      {
         fprintf( output_error[i], "#    Num            logrho           logtemp          logenergy          entropy            logpress               Ye\n", "w" );
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

   int  N_data[3] = {0};
   for ( int i_mode=0; i_mode<3; i_mode++ ) {
      int keymode;
      const int keymethod_NB = 1;
      const int keymethod_TL = 2;

      if      ( i_mode == 0 ) keymode = NUC_MODE_ENGY; // energy  mode
      else if ( i_mode == 1 ) keymode = NUC_MODE_ENTR; // entropy mode
      else if ( i_mode == 2 ) keymode = NUC_MODE_PRES; // pressure mode
      while ( !feof( input[i_mode] ) ) {
         real dens;
         real eps;
         real ye;
         real temp;
         real entr;
         real pres;
         real cs2;
         real munu;
         real useless;

         fscanf( input[i_mode], " %d %lf %lf %lf %lf %lf %lf %lf", &N_data[i_mode], &dens, &temp, &eps, &entr, &pres, &ye, &useless );
         if ( feof( input[i_mode] ) ){
            break;
         }


         dens = POW( 10.0, dens );
         temp = POW( 10.0, temp );
         eps  = POW( 10.0, eps ) - g_energy_shift;
         pres = POW( 10.0, pres );

         real eps0  = eps;
         real pres0 = pres;
         real entr0 = entr;
         real temp0 = temp;

         if ( keymode != NUC_MODE_TEMP ) temp = 1.2*temp; // slightly change the temperature for (ener/entropy/pressure) mode

         int keyerr1 = 0;
         int keyerr2 = 0;
         real rel_err = 0.0;
         nuc_eos_C_short( dens, &temp, ye, &eps, &entr, &pres, &cs2, &munu, g_energy_shift,
                          g_nrho, g_ntemp, g_nye, g_nrho_mode, g_nmode, g_nye_mode, g_alltables, g_alltables_mode,
                          g_logrho, g_logtemp, g_yes, g_logrho_mode,
                          g_logeps_mode, g_entr_mode, g_logprss_mode, g_yes_mode,
                          interpol_TL, interpol_other, keymode, &keyerr1, Tolerance );
         real temp1 = temp;
         real eps1  = NULL_REAL;
         real entr1 = NULL_REAL;
         real pres1 = NULL_REAL;
         if ( keyerr1 == 0 ) {
            nuc_eos_C_short( dens, &temp1, ye, &eps1, &entr1, &pres1, &cs2, &munu, g_energy_shift,
                             g_nrho, g_ntemp, g_nye, g_nrho_mode, g_nmode, g_nye_mode, g_alltables, g_alltables_mode,
                             g_logrho, g_logtemp, g_yes, g_logrho_mode,
                             g_logeps_mode, g_entr_mode, g_logprss_mode, g_yes_mode,
                             interpol_TL, interpol_other, NUC_MODE_TEMP, &keyerr2, Tolerance );
            real rel_err = 0.0;
            switch ( keymode )
            {
               case NUC_MODE_ENGY:
               {
                  rel_err = FABS(  LOG10( eps1/eps0 ) / LOG10( eps0+g_energy_shift )  );
                  // rel_err = FABS(  ( eps1-eps0 )/( eps0+g_energy_shift )  );
               } break;
               case NUC_MODE_ENTR:
               {
                  rel_err = FABS( ( entr1-entr0 )/entr0 );
               } break;
               case NUC_MODE_PRES:
               {
                  rel_err = FABS(  LOG10( pres1/pres0 ) )/ LOG10( pres0 );
                  // rel_err = FABS(  ( pres1-pres0 )/pres0  );
               } break;
            }
            if ( rel_err > 1e-16 ) {
               // error 2: relative errors between the FLASH1D Data and [i_mode] > 0.1%
               failcount2[i_mode]++;
               fprintf( output_error[i_mode+3], "%8d    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e\n",
                        failcount2[i_mode], LOG10( dens ), LOG10( temp0 ), LOG10( eps0 + g_energy_shift ), entr0, LOG10( pres0 ), ye, rel_err );
               rel_err_array[i_mode] = rel_err_array[i_mode] + rel_err;
               rel_err_number[i_mode] = rel_err_number[i_mode]+1;
               max_err[i_mode] = ( rel_err > max_err[i_mode] ) ? rel_err: max_err[i_mode];
            }
         }
         else if ( keyerr1 != 0 ) {
            // error 1: fail in finding energy from the [i_mode]
            failcount1[i_mode]++;
            fprintf( stdout, "N_data[%1d]: %8d, dens: %10.8e, temp: %10.8e, ener: %10.8e, entr: %10.8e, prss: %10.8e, Ye: %10.8e\n", i_mode, N_data[i_mode],
                     dens, temp, eps, entr, pres, ye );
            fprintf( stdout, "error code: %d\n", keyerr1 );
            fprintf( output_error[i_mode], "%8d    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e\n",
                     failcount1[i_mode], LOG10( dens ), LOG10( temp0 ), LOG10( eps0 + g_energy_shift ), LOG10( pres0 ), entr0, ye );
         }
      } // end while
   }


 // end while (N_data<Nmax)

   fprintf( stdout, "Maximum relative error in energy mode       :  %10.8E\n", max_err[0] );
   fprintf( stdout, "Mean    relative error in energy mode       :  %10.8E\n", rel_err_array[0]/rel_err_number[0] );
   fprintf( stdout, "-----------------------------------------------------------------------------\n" );
   fprintf( stdout, "Maximum relative error is entropy mode      :  %10.8E\n", max_err[1] );
   fprintf( stdout, "Mean    relative error is entropy mode      :  %10.8E\n", rel_err_array[1]/rel_err_number[1] );
   fprintf( stdout, "-----------------------------------------------------------------------------\n" );
   fprintf( stdout, "Maximum relative error is pressure mode     :  %10.8E\n", max_err[2] );
   fprintf( stdout, "Mean    relative error is pressure mode     :  %10.8E\n", rel_err_array[2]/rel_err_number[2] );
   fprintf( stdout, "-----------------------------------------------------------------------------\n" );

   for ( int i=0; i<3; i++ ) fprintf( stdout, "fail# %6d\n", failcount1[i] );

   return;

} // FUNCTION : nuc_eos_C_testing

#endif  // NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP
