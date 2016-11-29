#ifndef COSMO_MACROS
#define COSMO_MACROS


/************************************************/
/* Variables that can be changed at compilation */
/* time using, eg, the -D option for gcc.       */
/************************************************/

// physical box size (eg., in hubble units)
// eg; L = H_LEN_FRAC = N*dx
#ifndef H_LEN_FRAC
  #define H_LEN_FRAC 0.5
#endif

// compile using reference integrator?
#ifndef USE_REFERENCE_FRW
  #define USE_REFERENCE_FRW true
#endif

// Stencil order
#ifndef STENCIL_ORDER
  #define STENCIL_ORDER 4
#endif

#ifndef REFINEMENT_ORDER
  #define REFINEMENT_ORDER 1
#endif

// evolve shift as well? (if not, assumed to be zero)
#ifndef USE_BSSN_SHIFT
  #define USE_BSSN_SHIFT false
#endif

// Gamma-driver gauge settings (must turn on bssn_shift as well)
#ifndef USE_GAMMA_DRIVER
  #define USE_GAMMA_DRIVER false
#endif

// normalize conformal metric and time-derivative
#ifndef NORMALIZE_GAMMAIJ_AIJ
  #define NORMALIZE_GAMMAIJ_AIJ true
#endif

// Numerical Error Damping strength parameters
#ifndef BS_H_DAMPING_AMPLITUDE
  #define BS_H_DAMPING_AMPLITUDE 0.0
#endif
#ifndef JM_K_DAMPING_AMPLITUDE
  #define JM_K_DAMPING_AMPLITUDE 0.0
#endif

// Optionally exclude some second-order terms
#ifndef EXCLUDE_SECOND_ORDER_SMALL
  #define EXCLUDE_SECOND_ORDER_SMALL false
#endif
#ifndef EXCLUDE_SECOND_ORDER_FRW
  #define EXCLUDE_SECOND_ORDER_FRW false
#endif

// Optionally compile without raytracing, multigrid classes
#ifndef USE_MULTIGRID
  #define USE_MULTIGRID true
#endif
#ifndef USE_COSMOTRACE
  #define USE_COSMOTRACE true
#endif

//Potential types
#ifndef USE_COSMO_CONST_POTENTIAL
  #define USE_COSMO_CONST_POTENTIAL true
  #ifndef COSMO_CONST
    #define COSMO_CONST 0.0003
  #endif
#endif

/*****************************************/
/* Additional variable/macro definitions */
/*****************************************/

// not really tested:
#define USE_Z4c_DAMPING false
#if USE_Z4c_DAMPING
  #define Z4c_K1_DAMPING_AMPLITUDE 0.01
  #define Z4c_K2_DAMPING_AMPLITUDE 0.0
#else
  #define Z4c_K1_DAMPING_AMPLITUDE 0.0
  #define Z4c_K2_DAMPING_AMPLITUDE 0.0
#endif

#define STENCIL_CONCATENATOR(function, order) function ## order
#define STENCIL_EVALUATOR(function, order) STENCIL_CONCATENATOR(function, order)
#define STENCIL_ORDER_FUNCTION(function) STENCIL_EVALUATOR(function, STENCIL_ORDER)

#define PI (4.0*atan(1.0))
#define SIGN(x) (((x) < 0.0) ? -1 : ((x) > 0.0))
#define pw2(x) ((x)*(x))
#define C_RE(c) ((c)[0])
#define C_IM(c) ((c)[1])
#define ROUND_2_IDXT(f) ((idx_t)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))


// index with designated grid number
#define H_INDEX(i,j,k,nx,ny,nz) (((i+nx)%(nx))*(ny)*(nz) + ((j+ny)%(ny))*(nz) + (k+(nz))%(nz))




// RK4 method, using 4 "registers".  One for the "_p"revious step data, one
// for the data being "_a"ctively used for calculation, one for the
// Runge-Kutta "_c"oefficient being calculated, and lastly the "_f"inal
// result of the calculation.

#define VAR_CREATE(field)  \                                  
        boost::shared_ptr<pdat::CellVariable<real_t>> field
        
#define VAR_INIT(field)  \
        field = boost::shared_ptr<pdat::CellData<real_t>> (  \
          new pdat::CellVariable<real_t>(dim, "field",hier::IntVector::getOne(dim), 1))
        

#define RK4_PDATA_CREATE(field)  \                                  
        boost::shared_ptr<pdat::CellData<real_t>> field##_a_pdata;  \
        boost::shared_ptr<pdat::CellData<real_t>> field##_p_pdata;  \
        boost::shared_ptr<pdat::CellData<real_t>> field##_c_pdata;  \
        boost::shared_ptr<pdat::CellData<real_t>> field##_f_pdata

#define RK4_MDA_ACCESS_CREATE(field)                          \
  arr_t field##_a;  \
  arr_t field##_p;  \
  arr_t field##_c;  \
  arr_t field##_f

        
#define RK4_IDX_CREATE(name) \
        idx_t name##_a_idx, name##_c_idx, name##_p_idx, name##_s_idx


#define RK4_ARRAY_ALLOC(field) \
        level->allocatePatchData(bssnSim->field##_a_idx);  \
        level->allocatePatchData(bssnSim->field##_p_idx);  \
        level->allocatePatchData(bssnSim->field##_f_idx);  \
        level->allocatePatchData(bssnSim->field##_c_idx)


#define RK4_SET_LOCAL_VALUES(name) \
        bd->name = name##_a(bd->i, bd->j, bd->k);


#define GEN1_IDX_CREATE(name) \
        idx_t name##_a_idx

#define GEN1_PDATA_CREATE(field)  \                                  
        boost::shared_ptr<pdat::CellData<real_t>> field##_a_pdata;  \

#define GEN1_MDA_ACCESS_CREATE(field)                          \
  arr_t field##_a;   \


        
#define GEN1_SET_LOCAL_VALUES(name) \
        bd->name = name##_a(bd->i, bd->j, bd->k);


// macros for summing
#define COSMO_SUMMATION_1(MACRO) \
  ( MACRO(1) + MACRO(2) + MACRO(3) )

#define COSMO_SUMMATION_1_ARGS(MACRO, ...) \
  ( MACRO(1, __VA_ARGS__) + MACRO(2, __VA_ARGS__) + MACRO(3, __VA_ARGS__) )

#define COSMO_SUMMATION_2(MACRO) ( \
  MACRO(1, 1) + MACRO(1, 2) + MACRO(1, 3) \
  + MACRO(2, 1) + MACRO(2, 2) + MACRO(2, 3) \
  + MACRO(3, 1) + MACRO(3, 2) + MACRO(3, 3) \
  )

#define COSMO_SUMMATION_2_ARGS(MACRO, ...) ( \
  MACRO(1, 1, __VA_ARGS__) + MACRO(1, 2, __VA_ARGS__) + MACRO(1, 3, __VA_ARGS__) \
  + MACRO(2, 1, __VA_ARGS__) + MACRO(2, 2, __VA_ARGS__) + MACRO(2, 3, __VA_ARGS__) \
  + MACRO(3, 1, __VA_ARGS__) + MACRO(3, 2, __VA_ARGS__) + MACRO(3, 3, __VA_ARGS__) \
  )

#define COSMO_SUMMATION_3(MACRO) ( \
  MACRO(1, 1, 1) + MACRO(1, 1, 2) + MACRO(1, 1, 3) \
  + MACRO(1, 2, 1) + MACRO(1, 2, 2) + MACRO(1, 2, 3) \
  + MACRO(1, 3, 1) + MACRO(1, 3, 2) + MACRO(1, 3, 3) \
  + MACRO(2, 1, 1) + MACRO(2, 1, 2) + MACRO(2, 1, 3) \
  + MACRO(2, 2, 1) + MACRO(2, 2, 2) + MACRO(2, 2, 3) \
  + MACRO(2, 3, 1) + MACRO(2, 3, 2) + MACRO(2, 3, 3) \
  + MACRO(3, 1, 1) + MACRO(3, 1, 2) + MACRO(3, 1, 3) \
  + MACRO(3, 2, 1) + MACRO(3, 2, 2) + MACRO(3, 2, 3) \
  + MACRO(3, 3, 1) + MACRO(3, 3, 2) + MACRO(3, 3, 3) \
  )

#define COSMO_SUMMATION_3_ARGS(MACRO, ...) ( \
  MACRO(1, 1, 1, __VA_ARGS__) + MACRO(1, 1, 2, __VA_ARGS__) + MACRO(1, 1, 3, __VA_ARGS__) \
  + MACRO(1, 2, 1, __VA_ARGS__) + MACRO(1, 2, 2, __VA_ARGS__) + MACRO(1, 2, 3, __VA_ARGS__) \
  + MACRO(1, 3, 1, __VA_ARGS__) + MACRO(1, 3, 2, __VA_ARGS__) + MACRO(1, 3, 3, __VA_ARGS__) \
  + MACRO(2, 1, 1, __VA_ARGS__) + MACRO(2, 1, 2, __VA_ARGS__) + MACRO(2, 1, 3, __VA_ARGS__) \
  + MACRO(2, 2, 1, __VA_ARGS__) + MACRO(2, 2, 2, __VA_ARGS__) + MACRO(2, 2, 3, __VA_ARGS__) \
  + MACRO(2, 3, 1, __VA_ARGS__) + MACRO(2, 3, 2, __VA_ARGS__) + MACRO(2, 3, 3, __VA_ARGS__) \
  + MACRO(3, 1, 1, __VA_ARGS__) + MACRO(3, 1, 2, __VA_ARGS__) + MACRO(3, 1, 3, __VA_ARGS__) \
  + MACRO(3, 2, 1, __VA_ARGS__) + MACRO(3, 2, 2, __VA_ARGS__) + MACRO(3, 2, 3, __VA_ARGS__) \
  + MACRO(3, 3, 1, __VA_ARGS__) + MACRO(3, 3, 2, __VA_ARGS__) + MACRO(3, 3, 3, __VA_ARGS__) \
  )

#endif
