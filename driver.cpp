#include "NuclearEoS.h"


real  *h_EoS_Table[EOS_NTABLE_MAX];

int    g_nrho;
int    g_nye;
int    g_nrho_mode;
int    g_nmode;
int    g_nye_mode;
double g_energy_shift;

real  *g_alltables;
real  *g_alltables_mode;
real  *g_logrho;
real  *g_yes;
real  *g_logrho_mode;
real  *g_entr_mode;
real  *g_logprss_mode;
real  *g_yes_mode;

#if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
int    g_ntemp;
real  *g_logtemp;
real  *g_logeps_mode;
#elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
int    g_neps;
real  *g_logeps;
real  *g_logtemp_mode;
#endif

int interpol_other = NUC_INTERPOL_CUBIC;
int interpol_TL    = NUC_INTERPOL_CUBIC;


int main(void)
{


   // read nuclear eos table
# if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
      nuc_eos_C_ReadTable( "/data/ihsong/EoS_Table/SFHo_Tables/SFHo_temp_111_cub.h5" );
# elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
      nuc_eos_C_ReadTable( "../LS220_eps_111_111_linear.h5" );
# endif


   // self consistency testing
   fprintf( stdout, "-----------------------------------------------------------------------------\n" );
   fprintf( stdout, "Self Consistency Test\n" );
   #if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
   nuc_eos_C_testing_temp_self_consistency();
   // nuc_eos_C_testing_temp_self_consistency_2();
   #elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   nuc_eos_C_testing_engy_self_consistency();
   #endif
   fprintf( stdout, "-----------------------------------------------------------------------------\n\n" );


   // temperature error in energy-based table and driver
#  if ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
   fprintf( stdout, "-----------------------------------------------------------------------------\n" );
   fprintf( stdout, "Temperature Error Test in Energy-Based Driver\n" );
   nuc_eos_C_testing_engy_temp_err();
   fprintf( stdout, "-----------------------------------------------------------------------------\n" );
#  endif


   return 0;
}
