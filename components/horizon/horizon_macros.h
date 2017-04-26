#ifndef HORIZON_MACROS
#define HORIZON_MACROS


#define SCALAR_APPLY_TO_FIELDS_ARGS(function, ...)        \
  function(F, __VA_ARGS__);                           

#define SCALAR_APPLY_TO_FIELDS(function)          \
  function(F);                                

#define SCALAR_APPLY_TO_FIELDS_ARGS(function, ...)        \
  function(F, __VA_ARGS__);                           

#define SCALAR_APPLY_TO_GEN1_EXTRAS(function)          \
  function(F);                                

#define SCALAR_RK_EVOLVE_PT_FIELD(field)               \
  field##_s(i,j,k) = ev_##field(&bd, &sd, dx) * dt;

#define SCALAR_RK_EVOLVE_PT \
  SCALAR_APPLY_TO_FIELDS(SCALAR_RK_EVOLVE_PT_FIELD)



#endif
