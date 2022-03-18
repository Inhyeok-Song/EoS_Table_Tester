#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Macro.h"

#define NUC_TABLE_MODE NUC_TABLE_MODE_TEMP


#define NUC_TABLE_NVAR       16     // number of variables in the EoS table lookup
#define NUC_TABLE_NPTR       10     // number of table pointers to be sent to GPU


// auxiliary array indices
#define NUC_AUX_ESHIFT        0     // AuxArray_Flt: energy_shift
#define NUC_AUX_DENS2CGS      1     // AuxArray_Flt: convert density    to cgs
#define NUC_AUX_PRES2CGS      2     // AuxArray_Flt: convert pressure   to cgs
#define NUC_AUX_VSQR2CGS      3     // AuxArray_Flt: convert velocity^2 to cgs
#define NUC_AUX_PRES2CODE     4     // AuxArray_Flt: convert pressure   to code unit
#define NUC_AUX_VSQR2CODE     5     // AuxArray_Flt: convert velocity^2 to code unit
#define NUC_AUX_KELVIN2MEV    6     // AuxArray_Flt: convert kelvin     to MeV
#define NUC_AUX_MEV2KELVIN    7     // AuxArray_Flt: convert MeV        to kelvin
#define NUC_AUX_M_kB          8     // AuxArray_Flt: mean molecular weight*atomic mass unit/
                                    //               Bolzmann constant*(UNIT_E/UNIT_M)

#define NUC_AUX_NRHO          0     // AuxArray_Int: nrho
#define NUC_AUX_NTORE         1     // AuxArray_Int: ntemp/neps
#define NUC_AUX_NYE           2     // AuxArray_Int: nye
#define NUC_AUX_NRHO_MODE     3     // AuxArray_Int: nrho_mode
#define NUC_AUX_NMODE         4     // AuxArray_Int: nmode
#define NUC_AUX_NYE_MODE      5     // AuxArray_Int: nye_mode

// Interpolation schemes
#define NUC_AUX_INT_TL        6     // AuxArray_Int: interpolation scheme for table look-ups
#define NUC_AUX_INT_OTHER     7     // AuxArray_Int: interpolation scheme for other thermodynamic variables


// table indices
#define NUC_TAB_ALL           0     // alltables
#define NUC_TAB_ALL_MODE      1     // alltables_mode
#define NUC_TAB_RHO           2     // logrho
#define NUC_TAB_TORE          3     // logtemp/logenergy
#define NUC_TAB_YE            4     // yes
#define NUC_TAB_RHO_MODE      5     // logrho_mode
#define NUC_TAB_EORT_MODE     6     // logenergy_mode/logtemp_mode
#define NUC_TAB_ENTR_MODE     7     // entr_mode
#define NUC_TAB_PRES_MODE     8     // logprss_mode
#define NUC_TAB_YE_MODE       9     // yes_mode


// EoS modes
#define NUC_MODE_ENGY         0     // energy mode
#define NUC_MODE_TEMP         1     // temperature mode
#define NUC_MODE_ENTR         2     // entropy mode
#define NUC_MODE_PRES         3     // pressure mode


// Tolerance for Newton-Raphson or bisection method in temperature driver



// functions prototype
void nuc_eos_C_ReadTable( char *nuceos_table_name );

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
                      const int keymode, int *keyerr, const real rfeps );

void nuc_eos_C_cubinterp_some( const real x, const real y, const real z,
                               real *output_vars, const real *alltables,
                               const int nx, const int ny, const int nz, const int nvars,
                               const real *xt, const real *yt, const real *zt );

void nuc_eos_C_linterp_some( const real x, const real y, const real z,
                             real *output_vars, const real *alltables,
                             const int nx, const int ny, const int nz, const int nvars,
                             const real *xt, const real *yt, const real *zt );

void findtoreps( const real x, const real y, const real z,
                 real *found_lt, const real *alltables_mode,
                 const int nx, const int ny, const int nz, const int ntemp,
                 const real *xt, const real *yt, const real *zt, const real *logtoreps,
                 const int interpol_TL, const int keymode, int *keyerr );

void findtemp_NR_bisection( const real lr, const real lt0, const real ye, const real varin, real *ltout,
                            const int nrho, const int ntemp, const int nye, const real *alltables,
                            const real *logrho, const real *logtemp, const real *yes,
                            const int keymode, int *keyerrt, const real prec );


// testing
#if   ( NUC_TABLE_MODE == NUC_TABLE_MODE_TEMP )
void nuc_eos_C_testing_temp_self_consistency();
void nuc_eos_C_testing_temp_self_consistency_2();
void nuc_eos_C_testing_temp_self_consistency_whole_domain();
#elif ( NUC_TABLE_MODE == NUC_TABLE_MODE_ENGY )
void nuc_eos_C_testing_engy_self_consistency();
void nuc_eos_C_testing_engy_temp_err();
#endif
