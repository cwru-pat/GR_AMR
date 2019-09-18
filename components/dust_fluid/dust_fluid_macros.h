#ifndef DUST_FLUID_MACROS
#define DUST_FLUID_MACROS


#define DUST_FLUID_APPLY_TO_FIELDS_ARGS(function, ...)        \
  function(DF_D, __VA_ARGS__);                           \
  function(DF_S1, __VA_ARGS__);                            \
  function(DF_S2, __VA_ARGS__);                          \
  function(DF_S3, __VA_ARGS__);                          \
  function(DF_E, __VA_ARGS__);              

#define DUST_FLUID_APPLY_TO_FIELDS(function)          \
  function(DF_D);                                \
  function(DF_S1);                                 \
  function(DF_S2);                               \
  function(DF_S3);                               \
  function(DF_E);              


#define DUST_FLUID_APPLY_TO_DERIVED_FIELDS_ARGS(function, ...)  \
  function(F01, __VA_ARGS__);                           \
  function(F02, __VA_ARGS__);                           \
  function(F03, __VA_ARGS__);                           \
  function(F11, __VA_ARGS__);                           \
  function(F12, __VA_ARGS__);                           \
  function(F13, __VA_ARGS__);                           \
  function(F21, __VA_ARGS__);                           \
  function(F22, __VA_ARGS__);                           \
  function(F23, __VA_ARGS__);                           \
  function(F31, __VA_ARGS__);                           \
  function(F32, __VA_ARGS__);                           \
  function(F33, __VA_ARGS__);                           \
  function(F41, __VA_ARGS__);                           \
  function(F42, __VA_ARGS__);                           \
  function(F43, __VA_ARGS__);                           


#define DUST_FLUID_APPLY_TO_DERIVED_FIELDS(function)          \
  function(F01);                           \
  function(F02);                           \
  function(F03);                           \
  function(F11);                           \
  function(F12);                           \
  function(F13);                           \
  function(F21);                           \
  function(F22);                           \
  function(F23);                           \
  function(F31);                           \
  function(F32);                           \
  function(F33);                           \
  function(F41);                           \
  function(F42);                           \
  function(F43);                           




#define DUST_FLUID_RK_EVOLVE_PT_FIELD(field)               \
  field##_s(i,j,k) = ev_##field(&bd, &dd, dx) * dt;

#define DUST_FLUID_RK_EVOLVE_BD_FIELD(field)               \
  field##_s(i,j,k) = ev_##field##_bd(&bd, &dd, dx, l_idx, codim) * dt;


#define DUST_FLUID_RK_EVOLVE_PT \
  DUST_FLUID_APPLY_TO_FIELDS(DUST_FLUID_RK_EVOLVE_PT_FIELD)

#define DUST_FLUID_RK_EVOLVE_BD \
  DUST_FLUID_APPLY_TO_FIELDS(DUST_FLUID_RK_EVOLVE_BD_FIELD)



#endif
