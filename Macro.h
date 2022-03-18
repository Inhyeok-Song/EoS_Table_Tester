#define FLOAT8


// single/double precision
#ifdef FLOAT8
typedef double real;
#else
typedef float  real;
#endif


// max/min functions
#define MAX( a, b )     (  ( (a) > (b) ) ? (a) : (b)  )
#define MIN( a, b )     (  ( (a) < (b) ) ? (a) : (b)  )


// square/cube function
#define SQR(  a )       ( (a)*(a)     )
#define CUBE( a )       ( (a)*(a)*(a) )


// single/double-precision mathematic functions
#ifdef FLOAT8
#  define  FABS( a )        fabs( a )
#  define  SQRT( a )        sqrt( a )
#  define   SIN( a )         sin( a )
#  define   COS( a )         cos( a )
#  define   LOG( a )         log( a )
#  define LOG10( a )       log10( a )
#  define   EXP( a )         exp( a )
#  define  ATAN( a )        atan( a )
#  define FLOOR( a )       floor( a )
#  define  FMAX( a, b )     fmax( a, b )
#  define  FMIN( a, b )     fmin( a, b )
#  define   POW( a, b )      pow( a, b )
#  define  FMOD( a, b )     fmod( a, b )
#  define ATAN2( a, b )    atan2( a, b )
#else
#  define  FABS( a )        fabsf( a )
#  define  SQRT( a )        sqrtf( a )
#  define   SIN( a )         sinf( a )
#  define   COS( a )         cosf( a )
#  define   LOG( a )         logf( a )
#  define LOG10( a )       log10f( a )
#  define   EXP( a )         expf( a )
#  define  ATAN( a )        atanf( a )
#  define FLOOR( a )       floorf( a )
#  define  FMAX( a, b )     fmaxf( a, b )
#  define  FMIN( a, b )     fminf( a, b )
#  define   POW( a, b )      powf( a, b )
#  define  FMOD( a, b )     fmodf( a, b )
#  define ATAN2( a, b )    atan2f( a, b )
#endif


#ifndef NULL_REAL
#  define NULL_REAL        __FLT_MAX__
#endif


// nuclear eos table default mode
#define NUC_TABLE_MODE_TEMP 1
#define NUC_TABLE_MODE_ENGY 2


//nuclear eos table size
# define EOS_NAUX_MAX           20    // EoS_AuxArray_Flt/Int[]
# define EOS_NTABLE_MAX         20    // *_EoS_Table[]


// nuclear eos interpolation schemes
#define NUC_INTERPOL_LINEAR 1
#define NUC_INTERPOL_CUBIC  2

