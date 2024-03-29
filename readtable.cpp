#include <stdlib.h>
#include "NuclearEoS.h"


#  define H5_USE_16_API    1
#  include "hdf5.h"

#  ifdef FLOAT8
#    define H5T_GAMER_REAL H5T_NATIVE_DOUBLE
#  else
#    define H5T_GAMER_REAL H5T_NATIVE_FLOAT
#  endif


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

#if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
extern int    g_ntemp;
extern real  *g_logtemp;
extern real  *g_logeps_mode;
#elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
extern int    g_neps;
extern real  *g_logeps;
extern real  *g_logtemp_mode;
#endif



// catch HDF5 errors
#define HDF5_ERROR( fn_call )                                           \
{                                                                       \
   const int _error_code = fn_call;                                     \
   if ( _error_code < 0 )                                               \
   {                                                                    \
      fprintf( stderr,                                                  \
               "Could not read nuceos_table_name %s \n",                \
               nuceos_table_name );                                     \
      abort();                                                          \
   }                                                                    \
}



//-------------------------------------------------------------------------------------
// Function    :  nuc_eos_C_ReadTable
// Description :  Load the EoS table from the disk
//
// Note        :  1. Invoked by EoS_Init_Nuclear()
//
// Parameter   :  nuceos_table_name : Filename
//
// Return      :  EoS tables
//-------------------------------------------------------------------------------------
void nuc_eos_C_ReadTable( char *nuceos_table_name )
{


   fprintf( stdout, "*******************************\n"   );
   fprintf( stdout, "Reading nuc_eos table file:\n"       );
   fprintf( stdout, "%s\n",nuceos_table_name              );
   fprintf( stdout, "*******************************\n\n" );

// use these two macros to easily read in a lot of variables in the same way
// --> the first reads in one variable of a given type completely
#  define READ_EOS_HDF5( NAME, VAR, TYPE, MEM )                                        \
   {                                                                                   \
      hid_t dataset;                                                                   \
      HDF5_ERROR(  dataset = H5Dopen( file, NAME )  );                                 \
      HDF5_ERROR(  H5Dread( dataset, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR )  );        \
      HDF5_ERROR(  H5Dclose( dataset )  );                                             \
   }

#  define READ_EOSTABLE_HDF5( NAME, OFF )                                              \
   {                                                                                   \
      hsize_t offset[2] = { OFF, 0 };                                                  \
      H5Sselect_hyperslab( mem3, H5S_SELECT_SET, offset, NULL, var3, NULL );           \
      READ_EOS_HDF5( NAME, alltables_tmp, H5T_GAMER_REAL, mem3 );                      \
   }

#  define READ_EOSTABLE_MODE_HDF5( NAME,OFF )                                          \
   {                                                                                   \
      hsize_t offset[2] = { OFF, 0 };                                                  \
      H5Sselect_hyperslab( mem3_mode, H5S_SELECT_SET, offset, NULL, var3_mode, NULL ); \
      READ_EOS_HDF5( NAME, alltables_mode_tmp, H5T_GAMER_REAL, mem3_mode );            \
   }



// open file
   hid_t file;
   HDF5_ERROR(  file = H5Fopen( nuceos_table_name, H5F_ACC_RDONLY, H5P_DEFAULT )  );


// read size of tables
   READ_EOS_HDF5( "pointsrho",      &g_nrho,      H5T_NATIVE_INT, H5S_ALL );
#  if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
   READ_EOS_HDF5( "pointstemp",     &g_ntemp,     H5T_NATIVE_INT, H5S_ALL );
#  elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   READ_EOS_HDF5( "pointsenergy",   &g_neps,      H5T_NATIVE_INT, H5S_ALL );
#  endif
   READ_EOS_HDF5( "pointsye",       &g_nye,       H5T_NATIVE_INT, H5S_ALL );
   READ_EOS_HDF5( "pointsrho_mode", &g_nrho_mode, H5T_NATIVE_INT, H5S_ALL );
   READ_EOS_HDF5( "points_mode",    &g_nmode,     H5T_NATIVE_INT, H5S_ALL );
   READ_EOS_HDF5( "pointsye_mode",  &g_nye_mode,  H5T_NATIVE_INT, H5S_ALL );


#  if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
   int n_def_mode = g_ntemp;
#  elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   int n_def_mode = g_neps;
#  endif


// allocate memory for tables
   real *alltables_tmp      = NULL;
   real *alltables_mode_tmp = NULL;

   if (  ! ( alltables_tmp       = (real*)malloc(g_nrho*n_def_mode*g_nye*NUC_TABLE_NVAR*sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
   if (  ! ( g_alltables         = (real*)malloc(g_nrho*n_def_mode*g_nye*NUC_TABLE_NVAR*sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
   if (  ! ( alltables_mode_tmp  = (real*)malloc(g_nrho_mode*g_nmode*g_nye_mode*3      *sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
   if (  ! ( g_alltables_mode    = (real*)malloc(g_nrho_mode*g_nmode*g_nye_mode*3      *sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
   if (  ! ( g_logrho            = (real*)malloc(g_nrho                                *sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
   if (  ! ( g_yes               = (real*)malloc(g_nye                                 *sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
   if (  ! ( g_logrho_mode       = (real*)malloc(g_nrho_mode                           *sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
   if (  ! ( g_entr_mode         = (real*)malloc(g_nmode                               *sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
   if (  ! ( g_logprss_mode      = (real*)malloc(g_nmode                               *sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
   if (  ! ( g_yes_mode          = (real*)malloc(g_nye_mode                            *sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
#  if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
   if (  ! ( g_logtemp           = (real*)malloc(n_def_mode                            *sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
   if (  ! ( g_logeps_mode       = (real*)malloc(g_nmode                               *sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
#  elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   if (  ! ( g_logeps            = (real*)malloc(n_def_mode                            *sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
   if (  ! ( g_logtemp_mode      = (real*)malloc(g_nmode                               *sizeof(real)) )  )
   {
      fprintf( stderr, "Cannot allocate memory for EOS table\n" );
      abort();
   }
#  endif

// prepare HDF5 to read hyperslabs into alltables_tmp[]

   hsize_t table_dims[2]      = { NUC_TABLE_NVAR, g_nrho*n_def_mode*g_nye };
   hsize_t var3[2]            = { 1, g_nrho*n_def_mode*g_nye };
   hid_t   mem3               = H5Screate_simple( 2, table_dims, NULL );

   hsize_t table_dims_mode[2] = { 3, g_nrho_mode*g_nmode*g_nye_mode };
   hsize_t var3_mode[2]       = { 1, g_nrho_mode*g_nmode*g_nye_mode };
   hid_t   mem3_mode          = H5Screate_simple( 2, table_dims_mode, NULL );


// read alltables_tmp[]
   READ_EOSTABLE_HDF5( "logpress",  0 );
#  if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
   READ_EOSTABLE_HDF5( "logenergy", 1 );
#  elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   READ_EOSTABLE_HDF5( "logtemp",   1 );
#  endif
   READ_EOSTABLE_HDF5( "entropy",   2 );
   READ_EOSTABLE_HDF5( "munu",      3 );
   READ_EOSTABLE_HDF5( "cs2",       4 );

// chemical potentials
   READ_EOSTABLE_HDF5( "muhat",     5 );
   READ_EOSTABLE_HDF5( "mu_e",      6 );
   READ_EOSTABLE_HDF5( "mu_p",      7 );
   READ_EOSTABLE_HDF5( "mu_n",      8 );

// compositions
   READ_EOSTABLE_HDF5( "Xa",        9 );
   READ_EOSTABLE_HDF5( "Xh",       10 );
   READ_EOSTABLE_HDF5( "Xn",       11 );
   READ_EOSTABLE_HDF5( "Xp",       12 );

// average nucleus
   READ_EOSTABLE_HDF5( "Abar",     13 );
   READ_EOSTABLE_HDF5( "Zbar",     14 );

// Gamma
   READ_EOSTABLE_HDF5( "gamma",    15 );

// energy for temp, entr modes
#  if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
   READ_EOSTABLE_MODE_HDF5( "logtemp_ener",   0 );
   READ_EOSTABLE_MODE_HDF5( "logtemp_entr",   1 );
   READ_EOSTABLE_MODE_HDF5( "logtemp_prss",   2 );
#  elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   READ_EOSTABLE_MODE_HDF5( "logenergy_temp", 0 );
   READ_EOSTABLE_MODE_HDF5( "logenergy_entr", 1 );
   READ_EOSTABLE_MODE_HDF5( "logenergy_prss", 2 );
#  endif

// read additional tables and variables
   READ_EOS_HDF5( "logrho",         g_logrho,        H5T_GAMER_REAL,    H5S_ALL );
   READ_EOS_HDF5( "ye",             g_yes,           H5T_GAMER_REAL,    H5S_ALL );
   READ_EOS_HDF5( "logrho_mode",    g_logrho_mode,   H5T_GAMER_REAL,    H5S_ALL );
   READ_EOS_HDF5( "entropy_mode",   g_entr_mode,     H5T_GAMER_REAL,    H5S_ALL );
   READ_EOS_HDF5( "logpress_mode",  g_logprss_mode,  H5T_GAMER_REAL,    H5S_ALL );
   READ_EOS_HDF5( "ye_mode",        g_yes_mode,      H5T_GAMER_REAL,    H5S_ALL );
   READ_EOS_HDF5( "energy_shift",  &g_energy_shift,  H5T_NATIVE_DOUBLE, H5S_ALL );
#  if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
   READ_EOS_HDF5( "logtemp",        g_logtemp,       H5T_GAMER_REAL,    H5S_ALL );
   READ_EOS_HDF5( "logenergy_mode", g_logeps_mode,   H5T_GAMER_REAL,    H5S_ALL );
#  elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   READ_EOS_HDF5( "logenergy",      g_logeps,        H5T_GAMER_REAL,    H5S_ALL );
   READ_EOS_HDF5( "logtemp_mode",   g_logtemp_mode,  H5T_GAMER_REAL,    H5S_ALL );
#  endif

   HDF5_ERROR(  H5Sclose( mem3      )  );
   HDF5_ERROR(  H5Sclose( mem3_mode )  );
   HDF5_ERROR(  H5Fclose( file      )  );


// change ordering of g_alltables[] so that the table kind is the fastest changing index
   for ( int iv=0; iv<NUC_TABLE_NVAR; iv++ )
   for ( int k=0;  k<g_nye;           k++  )
   for ( int j=0;  j<n_def_mode;      j++  )
   for ( int i=0;  i<g_nrho;          i++  )
   {
      const long indold = i + g_nrho*( j + n_def_mode*(k + g_nye*iv) );
      const long indnew = iv + NUC_TABLE_NVAR*( i + g_nrho*(j + n_def_mode*k) );

      g_alltables[indnew] = alltables_tmp[indold];
   }

   for ( int iv=0; iv<3;          iv++ )
   for ( int k=0;  k<g_nye_mode;  k++  )
   for ( int j=0;  j<g_nmode;     j++  )
   for ( int i=0;  i<g_nrho_mode; i++  )
   {
      const long indold = i + g_nrho_mode*( j + g_nmode*(k + g_nye_mode*iv) );
      const long indnew = iv + 3*( i + g_nrho_mode*(j + g_nmode*k) );

      g_alltables_mode[indnew] = alltables_mode_tmp[indold];
   }


// free memory of temporary arrays
   free( alltables_tmp      );
   free( alltables_mode_tmp );


// set the EoS table pointers
   h_EoS_Table[NUC_TAB_ALL      ] = g_alltables;
   h_EoS_Table[NUC_TAB_ALL_MODE ] = g_alltables_mode;
   h_EoS_Table[NUC_TAB_RHO      ] = g_logrho;
#  if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
   h_EoS_Table[NUC_TAB_TORE     ] = g_logtemp;
#  elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   h_EoS_Table[NUC_TAB_TORE     ] = g_logeps;
#  endif
   h_EoS_Table[NUC_TAB_YE       ] = g_yes;
   h_EoS_Table[NUC_TAB_RHO_MODE ] = g_logrho_mode;
#  if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
   h_EoS_Table[NUC_TAB_EORT_MODE] = g_logeps_mode;
#  elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   h_EoS_Table[NUC_TAB_EORT_MODE] = g_logtemp_mode;
#  endif
   h_EoS_Table[NUC_TAB_ENTR_MODE] = g_entr_mode;
   h_EoS_Table[NUC_TAB_PRES_MODE] = g_logprss_mode;
   h_EoS_Table[NUC_TAB_YE_MODE  ] = g_yes_mode;



} // FUNCTION : nuc_eos_C_ReadTable
