#ifndef COSMO_MACROS
#define COSMO_MACROS

#include "cosmo_types.h"
/************************************************/
/* Variables that can be changed at compilation */
/* time using, eg, the -D option for gcc.       */
/************************************************/


#ifndef DIM
  #define DIM 3
#endif

#ifndef STENCIL_ORDER
  #define STENCIL_ORDER 4
#endif

#ifndef HALF_STENCIL_ORDER
  #define HALF_STENCIL_ORDER 2
#endif

#ifndef STENCIL_ORDER_WIDTH
  #define STENCIL_ORDER_WIDTH STENCIL_ORDER
#endif


#ifndef REFINEMENT_ORDER
  #define REFINEMENT_ORDER 1
#endif

#ifndef USE_SOMMERFIELD_BOUNDARY
  #define USE_SOMMERFIELD_BOUNDARY 1
#endif

#ifndef USE_EXPANSION
  #define USE_EXPANSION 0
#endif

// evolve shift as well? (if not, assumed to be zero)
#ifndef USE_BSSN_SHIFT
  #define USE_BSSN_SHIFT true
#endif

// Gamma-driver gauge settings (must turn on bssn_shift as well)
#ifndef USE_GAMMA_DRIVER
  #define USE_GAMMA_DRIVER true
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
#ifndef USE_Z4C
  #define USE_Z4C false
#endif


#define STENCIL_CONCATENATOR(function, order) function ## order
#define STENCIL_EVALUATOR(function, order) STENCIL_CONCATENATOR(function, order)
#define STENCIL_ORDER_FUNCTION(function) STENCIL_EVALUATOR(function, STENCIL_ORDER)
#define HALF_STENCIL_ORDER_FUNCTION(function) STENCIL_EVALUATOR(function, HALF_STENCIL_ORDER)

#define PI (4.0*atan(1.0))
#define SIGN(x) (((x) < 0.0) ? -1 : ((x) > 0.0))
#define pw2(x) ((x)*(x))
#define pw3(x) ((x)*(x)*(x))
#define C_RE(c) ((c)[0])
#define C_IM(c) ((c)[1])
#define ROUND_2_IDXT(f) ((idx_t)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))


#define DECLARE_REAL_T(name) \
  real_t name


#define VAR_CREATE(field)                               \
  boost::shared_ptr<pdat::CellVariable<real_t>> field

#define VAR_INIT(field)                                         \
  field = boost::shared_ptr<pdat::CellVariable<real_t>> (       \
    new pdat::CellVariable<real_t>(dim, #field, 1))

#define RK4_PDATA_CREATE(field, type)                                   \
  boost::shared_ptr<pdat::CellData<real_t>> field##_##type##_pdata

#define RK4_MDA_ACCESS_CREATE(field,type)       \
  arr_t field##_##type

#define RK4_PDATA_ALL_CREATE(field)                             \
  boost::shared_ptr<pdat::CellData<real_t>> field##_a_pdata;    \
  boost::shared_ptr<pdat::CellData<real_t>> field##_s_pdata;    \
  boost::shared_ptr<pdat::CellData<real_t>> field##_p_pdata;    \
  boost::shared_ptr<pdat::CellData<real_t>> field##_k1_pdata;   \
  boost::shared_ptr<pdat::CellData<real_t>> field##_k2_pdata;   \
  boost::shared_ptr<pdat::CellData<real_t>> field##_k3_pdata;   \
  boost::shared_ptr<pdat::CellData<real_t>> field##_k4_pdata

#define RK4_MDA_ACCESS_ALL_CREATE(field)        \
  arr_t field##_a;                              \
  arr_t field##_s;                              \
  arr_t field##_p;                              \
  arr_t field##_k1;                             \
  arr_t field##_k2;                             \
  arr_t field##_k3;                             \
  arr_t field##_k4

// RK4 method, using 4 "registers" plus _p represents "previous,
// _a represents "active", _s represents "scratch".
// for the data being "_a"ctively used for calculation
// "_p" stores results on previous step 
// "_s" is used when you need a temperary register to store somthing
#define RK4_IDX_ALL_CREATE(name)                                        \
  idx_t name##_a_idx, name##_s_idx, name##_p_idx, name##_k1_idx,        \
    name##_k2_idx, name##_k3_idx, name##_k4_idx

#define RK4_IDX_CREATE(name, type)              \
  idx_t name##_##type##_idx

#define REG_TO_CONTEXT(field, context, name, width)  \
  field##_##name##_idx =  \
    variable_db->registerVariableAndContext(  \
      field,  \
      context,  \
      hier::IntVector(dim, width))



#define RK4_ARRAY_ZERO(field, hcellmath)                \
  hcellmath.setToScalar(field##_a_idx, 0, 0);  \
  hcellmath.setToScalar(field##_s_idx, 0, 0);  \
  hcellmath.setToScalar(field##_p_idx, 0, 0);  \
  hcellmath.setToScalar(field##_k1_idx, 0, 0); \
  hcellmath.setToScalar(field##_k2_idx, 0, 0); \
  hcellmath.setToScalar(field##_k3_idx, 0, 0); \
  hcellmath.setToScalar(field##_k4_idx, 0, 0)

#define RK4_ARRAY_ALLOC(field)                 \
  level->allocatePatchData(field##_a_idx);     \
  level->allocatePatchData(field##_s_idx);     \
  level->allocatePatchData(field##_p_idx);     \
  level->allocatePatchData(field##_k1_idx);    \
  level->allocatePatchData(field##_k2_idx);    \
  level->allocatePatchData(field##_k3_idx);    \
  level->allocatePatchData(field##_k4_idx)

#define EXTRA_ARRAY_ZERO(field, hcellmath)                \
  hcellmath.setToScalar(field##_a_idx, 0, 0);  


#define EXTRA_ARRAY_ALLOC(field)                        \
  level->allocatePatchData(field##_a_idx)



#define ADD_TO_LIST(field, list)        \
  list.push_back(field##_a_idx)

#define RK4_SET_LOCAL_VALUES(name)              \
  bd->name = name##_a(bd->i, bd->j, bd->k);

//only add active field idx
#define ADD_VAR_TO_LIST(field)                          \
  variable_id_list.push_back(bssnSim->field##_a_idx)

#define GEN1_IDX_CREATE(name)                   \
  idx_t name##_a_idx

#define GEN1_PDATA_CREATE(field)                                \
  boost::shared_ptr<pdat::CellData<real_t>> field##_a_pdata;    \

#define GEN1_MDA_ACCESS_CREATE(field)           \
  arr_t field##_a;                              \



#define GEN1_SET_LOCAL_VALUES(name)             \
  bd->name = name##_a(bd->i, bd->j, bd->k);

#define RK4_INIT_L_K1(field)  \
  field##_s(i,j,k) =  field##_k1(i,j,k) / 2.0

#define RK4_INIT_L_K2(field)  \
  field##_s(i,j,k) =         \
    (3.0*field##_k1(i,j,k) + 4.0*field##_k2(i,j,k)  \
     + 2.0*field##_k3(i,j,k) - field##_k4(i,j,k))/16.0

#define RK4_INIT_L_K3(field)  \
    field##_s(i,j,k) =         \
      (3.0*field##_k1(i,j,k) + 2.0*field##_k2(i,j,k) \
       + 4.0*field##_k3(i,j,k) - field##_k4(i,j,k))/16.0

#define RK4_INIT_L_K4(field)  \
    field##_s(i,j,k) =         \
      (field##_k2(i,j,k) + field##_k3(i,j,k))/4.0


#define RK4_INIT_R_K1(field)  \
  field##_s(i,j,k) =  (field##_k2(i,j,k) + field##_k3(i,j,k)) / 4.0

#define RK4_INIT_R_K2(field)  \
  field##_s(i,j,k) =         \
    (-field##_k1(i,j,k) + 4.0*field##_k2(i,j,k)  \
     + 2.0*field##_k3(i,j,k) + 3.0*field##_k4(i,j,k))/16.0

#define RK4_INIT_R_K3(field)  \
    field##_s(i,j,k) =         \
      (-field##_k1(i,j,k) + 2.0*field##_k2(i,j,k) \
       + 4.0*field##_k3(i,j,k) + 3.0*field##_k4(i,j,k))/16.0

#define RK4_INIT_R_K4(field)  \
    field##_s(i,j,k) = field##_k4(i,j,k)/2.0

#define REGISTER_SPACE_REFINE_A(field, refiner,refine_op)       \
  refiner.registerRefine(field##_a_idx,                         \
                         field##_a_idx,                         \
                         field##_a_idx,                         \
                         refine_op)

#define REGISTER_SPACE_REFINE_S(field, refiner,refine_op)       \
  refiner.registerRefine(field##_s_idx,                         \
                         field##_s_idx,                         \
                         field##_s_idx,                         \
                         refine_op)

#define REGISTER_COARSEN_A(field,coarsener,coarsen_op)       \
  coarsener.registerCoarsen(field##_a_idx,                   \
                            field##_a_idx,                   \
                            coarsen_op,                      \
                            NULL)


#define RK4_FINALIZE_FIELD_1(field)      \
  field##_k1(i,j,k) =  field##_s(i,j,k);  \
  field##_a(i,j,k) = field##_p(i,j,k) + field##_k1(i,j,k)/2.0  

#define RK4_FINALIZE_FIELD_2(field)              \
  field##_k2(i,j,k) = field##_s(i,j,k);  \
  field##_a(i,j,k) = field##_p(i,j,k) + field##_k2(i,j,k)/2.0  

#define RK4_FINALIZE_FIELD_3(field) \
  field##_k3(i,j,k) = field##_s(i,j,k);  \
  field##_a(i,j,k) = field##_p(i,j,k) + field##_k3(i,j,k)  

#define RK4_FINALIZE_FIELD_4(field) \
  field##_k4(i,j,k) = field##_s(i,j,k);    \
  field##_a(i,j,k) =  field##_p(i,j,k) +                                  \
    (field##_k4(i,j,k) + 2.0*field##_k3(i,j,k) + 2.0*field##_k2(i,j,k) + field##_k1(i,j,k))/6.0  

#define COPY_A_TO_P(field)  \
  hcellmath.copyData(field##_p_idx, field##_a_idx, 0)

#define PDATA_INIT(field, type)   \
  field##_##type##_pdata =  \
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  \
      patch->getPatchData(field##_##type##_idx))

#define MDA_ACCESS_INIT(field,type)   \
  field##_##type = pdat::ArrayDataAccess::access<DIM, double>(  \
    field##_##type##_pdata->getArrayData())


#define PDATA_ALL_INIT(field)  \
  field##_p_pdata =  \
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  \
      patch->getPatchData(field##_p##_idx));  \
  field##_s_pdata =  \
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  \
      patch->getPatchData(field##_s##_idx));  \
  field##_a_pdata =  \
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  \
      patch->getPatchData(field##_a##_idx));  \
  field##_k1_pdata =  \
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  \
      patch->getPatchData(field##_k1##_idx));             \
  field##_k2_pdata =                                      \
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(    \
      patch->getPatchData(field##_k2##_idx));               \
  field##_k3_pdata =  \
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  \
      patch->getPatchData(field##_k3##_idx));             \
  field##_k4_pdata =  \
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(  \
      patch->getPatchData(field##_k4##_idx))




#define MDA_ACCESS_ALL_INIT(field)  \
  field##_p = pdat::ArrayDataAccess::access<DIM, double>(  \
    field##_p_pdata->getArrayData());                    \
  field##_a = pdat::ArrayDataAccess::access<DIM, double>(  \
    field##_a_pdata->getArrayData());                    \
  field##_s = pdat::ArrayDataAccess::access<DIM, double>(  \
    field##_s_pdata->getArrayData());                    \
  field##_k1 = pdat::ArrayDataAccess::access<DIM, double>(  \
    field##_k1_pdata->getArrayData());                    \
  field##_k2 = pdat::ArrayDataAccess::access<DIM, double>(  \
    field##_k2_pdata->getArrayData());                    \
  field##_k3 = pdat::ArrayDataAccess::access<DIM, double>(  \
    field##_k3_pdata->getArrayData());                    \
  field##_k4 = pdat::ArrayDataAccess::access<DIM, double>(  \
    field##_k4_pdata->getArrayData())                  

#define SET_PATCH_TIME(field, from_t, to_t)        \
  patch->getPatchData(field##_a_idx)->setTime(to_t); \
  patch->getPatchData(field##_p_idx)->setTime(from_t)  \

#define SET_LEVEL_TIME(field, from_t, to_t) \
  level->setTime(to_t, field##_a_idx);		    \
  level->setTime(from_t, field##_p_idx)			    \


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
