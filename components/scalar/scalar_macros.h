#ifndef SCALAR_MACROS
#define SCALAR_MACROS


#define SCALAR_APPLY_TO_FIELDS_ARGS(function, ...)        \
  function(phi, __VA_ARGS__);                           \
  function(Pi, __VA_ARGS__);                            

#define SCALAR_APPLY_TO_FIELDS(function)          \
  function(phi);                                \
  function(Pi);                                 

#define SCALAR_RK_EVOLVE_PT_FIELD(field)               \
  field##_s(i,j,k) = ev_##field(&bd, &sd, dx) * dt;

#define SCALAR_RK_EVOLVE_BD_FIELD(field)               \
  field##_s(i,j,k) = ev_##field##_bd(&bd, &sd, dx, l_idx, codim) * dt;


#define SCALAR_RK_EVOLVE_PT \
  SCALAR_APPLY_TO_FIELDS(SCALAR_RK_EVOLVE_PT_FIELD)

#define SCALAR_RK_EVOLVE_BD \
  SCALAR_APPLY_TO_FIELDS(SCALAR_RK_EVOLVE_BD_FIELD)



#endif
