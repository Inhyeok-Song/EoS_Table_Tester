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

extern int    g_neps;
extern real  *g_logeps;
extern real  *g_logtemp_mode;

extern int interpol_TL;
extern int interpol_other;


#if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
//-------------------------------------------------------
// Function    :  nuc_eos_C_testing
// Description :
//-------------------------------------------------------
void nuc_eos_C_testing_engy_self_consistency()
{

   // test FLASH data with new eps energy-based driver
   int   N_data   = 0;
   const double KelvinToMeV = 8.61732814974493E-11;


   // output files for writing errors
   FILE *output_error[6];

   FILE *ener_err1 = fopen( "ENGY_temp_err1.txt", "w" );
   FILE *ener_err2 = fopen( "ENGY_temp_err2.txt", "w" );

   FILE *entr_err1 = fopen( "ENGY_entr_err1.txt", "w" );
   FILE *entr_err2 = fopen( "ENGY_entr_err2.txt", "w" );

   FILE *prss_err1 = fopen( "ENGY_prss_err1.txt", "w" );
   FILE *prss_err2 = fopen( "ENGY_prss_err2.txt", "w" );

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
         fprintf( output_error[i], "#    Num            logrho           logtemp          logenergy          logpress           entropy                Ye\n", "w" );
      }
      else if ( i>=3 )
      {
         fprintf( output_error[i], "#    Num            logrho           logtemp          logenergy          logpress           entropy                Ye         rel_error\n", "w" );
      }
   }

   // variables for writing errors
   int  failcount1[4] = {0};   // count error1: fail in finding temperature
   int  failcount2[4] = {0};   // count error2: relative error > 10^-3
   real max_err[4] = {0.0};


// Use these two defines to easily read in a lot of variables in the same way
// The first reads in one variable of a given type completely
#define READ_FLASH_HDF5( NAME,VAR,TYPE,MEM )                                 \
  do {                                                                       \
     hid_t dataset;                                                          \
     HDF5_ERROR( dataset = H5Dopen( file, NAME ) );                          \
     HDF5_ERROR( H5Dread( dataset, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR ) ); \
     HDF5_ERROR( H5Dclose( dataset ) );                                      \
  } while (0)
#define READ_FLASHDATA_HDF5( NAME, OFF )                                     \
  do {                                                                       \
     hsize_t offset[2]     = { OFF, 0 };                                     \
     H5Sselect_hyperslab( mem3, H5S_SELECT_SET, offset, NULL, var3, NULL );  \
     READ_FLASH_HDF5( NAME, alldata, H5T_NATIVE_DOUBLE, mem3 );              \
  } while (0)

#define READ_FLASH_NODE_HDF5( NAME, VAR, TYPE, MEM)                          \
  do {                                                                       \
     hid_t dataset;                                                          \
     HDF5_ERROR( dataset = H5Dopen( file, NAME ) );                          \
     HDF5_ERROR( H5Dread( dataset, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR ) ); \
     HDF5_ERROR( H5Dclose( dataset ) );                                      \
  } while (0)


   real eps_min  = 100.0;
   real eps_max  = 0.0;
   real entr_min = 100.0;
   real entr_max = 0.0;
   real pres_min = 100.0;
   real pres_max = 0.0;

#ifdef FLOAT8
   const real Tolerance = 1e-10;
#else
   const real Tolerance = 1e-6;
#endif


   // read data
   while ( N_data < 927 ) {

      // read file
      hid_t file;
      hid_t dset;
      hid_t integer_scalars_tid;
      int integer_scalars[DIM0];

      for ( int i=0; i<DIM0; i++ ) integer_scalars[i] = 0;

      // HDF5 file path
      char filepath[512] = "/data/ceag/gamer/201229_s20_ls220/";
      //char filepath[512] = "/data/ceag/gamer/190708_s20v_sfho/";
      char filename[512];
      sprintf( filename, "ccsn1d_hdf5_chk_%04d", N_data );
      strcat( filepath, filename );

      if ( !file_is_readable(filepath) )
      {
        fprintf(stderr, "Could not read nuceos_table_name %s\n",
                filepath);
        abort();
      }
      HDF5_ERROR( file = H5Fopen( filepath, H5F_ACC_RDONLY, H5P_DEFAULT ) );
      HDF5_ERROR( dset = H5Dopen( file, DATASET ) );
      HDF5_ERROR( integer_scalars_tid = H5Tcreate( H5T_COMPOUND, sizeof(int) ) );
      HDF5_ERROR( H5Tinsert( integer_scalars_tid, "value", 0, H5T_NATIVE_INT ) );
      HDF5_ERROR( H5Dread( dset, integer_scalars_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, integer_scalars ) );

      // Prepare HDF5 to read hyperslabs
      const int ix = integer_scalars[4]; //globalnumblock
      const int iy = integer_scalars[0];

      hsize_t table_dims[2]      = {6, ix*iy};
      hsize_t table_dims_node[2] = {1, ix};
      hsize_t var3[2]            = {1, ix*iy};
      hsize_t var3_node[2]       = {1, ix};
      hid_t   mem3               = H5Screate_simple( 2, table_dims,      NULL );
      hid_t   mem3_node          = H5Screate_simple( 2, table_dims_node, NULL );


      real *alldata;
      if ( !(alldata = (real*)malloc( ix*iy*6*sizeof(real) ) ) ) abort();

      int *node;
      if ( !(node = (int*)malloc( ix * 1 * sizeof(int) ) ) ) abort();


      READ_FLASHDATA_HDF5( "dens", 0 );
      READ_FLASHDATA_HDF5( "eint", 1 );
      READ_FLASHDATA_HDF5( "ye  ", 2 );
      READ_FLASHDATA_HDF5( "temp", 3 );
      READ_FLASHDATA_HDF5( "entr", 4 );
      READ_FLASHDATA_HDF5( "pres", 5 );
      READ_FLASH_NODE_HDF5( "node type", node, H5T_NATIVE_INT, H5S_ALL );

      real Dens[ix][iy];
      real Ener[ix][iy];
      real Ye[ix][iy];
      real Temp[ix][iy];
      real Entr[ix][iy];
      real Pres[ix][iy];

      for ( int j=0; j<iy; j++ ) { //integer_scalars[0]
         for ( int i=0; i<ix; i++ ) {
            Dens[i][j] = alldata[j + iy*(i + ix*0)];
            Ener[i][j] = alldata[j + iy*(i + ix*1)] - g_energy_shift;
              Ye[i][j] = alldata[j + iy*(i + ix*2)];
            Temp[i][j] = alldata[j + iy*(i + ix*3)];
            Entr[i][j] = alldata[j + iy*(i + ix*4)];
            Pres[i][j] = alldata[j + iy*(i + ix*5)];
         }
      }


      for ( int i_mode=0; i_mode<3; i_mode++ ) {
         int keymode;
         if      ( i_mode == 0 ) keymode = NUC_MODE_TEMP; // energy  mode
         else if ( i_mode == 1 ) keymode = NUC_MODE_ENTR; // entropy mode
         else if ( i_mode == 2 ) keymode = NUC_MODE_PRES; // pressure mode

         for ( int j=0; j<iy; j++ ) { //integer_scalars[0]
            for ( int i=0; i<ix; i++ ) {
               real rand1 = 0.04*rand()/(real)RAND_MAX + 0.98;
               real rand2 = 0.04*rand()/(real)RAND_MAX + 0.98;
               real rand3 = 0.04*rand()/(real)RAND_MAX + 0.98;
               real rand4 = 0.04*rand()/(real)RAND_MAX + 0.98;
               real rand5 = 0.04*rand()/(real)RAND_MAX + 0.98;
               real rand6 = 0.04*rand()/(real)RAND_MAX + 0.98;


               real dens = Dens[i][j];
               real eps  = Ener[i][j];
               real ye   = Ye[i][j];
               real temp = Temp[i][j]*KelvinToMeV;
               real entr = Entr[i][j];
               real pres = Pres[i][j];
               real cs2  = 0.0;
               real munu = 0.0;

               real eps0  = eps;
               real pres0 = pres;
               real entr0 = entr;
               real temp0 = temp;

               eps_min  = MIN( eps_min, LOG10( eps + g_energy_shift ) );
               eps_max  = MAX( eps_max, LOG10( eps + g_energy_shift ) );
               entr_min = MIN( entr_min, entr );
               entr_max = MAX( entr_max, entr );
               pres_min = MIN( pres_min, LOG10( pres ) );
               pres_max = MAX( pres_max, LOG10( pres ) );

               int keyerr = 0;
               if ( node[i] < 2 ) {
                  real rel_err = 0.0;
                  nuc_eos_C_short( dens, &temp, ye, &eps, &entr, &pres, &cs2, &munu, g_energy_shift,
                                   g_nrho, g_neps, g_nye, g_nrho_mode, g_nmode, g_nye_mode, g_alltables, g_alltables_mode,
                                   g_logrho, g_logeps, g_yes, g_logrho_mode,
                                   g_logtemp_mode, g_entr_mode, g_logprss_mode, g_yes_mode,
                                   interpol_TL, interpol_other, keymode, &keyerr, Tolerance );
                  rel_err = fabs( ( eps - eps0 )/ ( eps0+g_energy_shift ) );

                  if ( keyerr != 0 ) {
                     // error 1: fail in finding energy from the [i_mode]
                     failcount1[i_mode]++;
                     fprintf( output_error[i_mode], "%8d    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e\n",
                              failcount1[i_mode], LOG10( dens ), LOG10( temp0 ), LOG10( eps0 + g_energy_shift ), LOG10( pres0 ), entr0, ye );
                  }
                  else if ( rel_err > 1.e-10 ) { // ) {
                     // error 2: relative errors between the FLASH1D Data and [i_mode] > 0.1%
                     failcount2[i_mode]++;
                     fprintf( output_error[i_mode+3], "%8d    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e    %10.8e\n",
                              failcount2[i_mode], LOG10( dens ), LOG10( temp0 ), LOG10( eps0 + g_energy_shift ), entr0, LOG10( pres0 ), ye, rel_err );
                     max_err[i_mode] = ( rel_err > max_err[i_mode] ) ? rel_err: max_err[i_mode];
                  }
               }
            }
         }
      }

      HDF5_ERROR( H5Sclose( mem3 ) );
      HDF5_ERROR( H5Tclose( integer_scalars_tid ) );
      HDF5_ERROR( H5Dclose( dset ) );
      HDF5_ERROR( H5Fclose( file ) );

      free( alldata );

      N_data++;

   } // end while (N_data<Nmax)

   printf( "min eps : %9.7f\n", eps_min  );
   printf( "max eps : %9.7f\n", eps_max  );
   printf( "min entr: %9.7f\n", entr_min );
   printf( "max entr: %9.7f\n", entr_max );
   printf( "min pres: %9.7f\n", pres_min );
   printf( "max pres: %9.7f\n", pres_max );

   fprintf( stdout, "Maximum relative error in temperature     :  %f\n", max_err[0] );
   fprintf( stdout, "-----------------------------------------------------------------------------\n" );
   fprintf( stdout, "Maximum relative error is entropy mode    :  %f\n", max_err[1] );
   fprintf( stdout, "-----------------------------------------------------------------------------\n" );
   fprintf( stdout, "Maximum relative error is pressure mode   :  %f\n", max_err[2] );
   fprintf( stdout, "-----------------------------------------------------------------------------\n" );
   for ( int i=0; i<3; i++ ) fprintf( stdout, "fail# %6d\n", failcount1[i] );

   return;

} // FUNCTION : nuc_eos_C_testing

#endif  // NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY
