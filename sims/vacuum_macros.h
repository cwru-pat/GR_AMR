#ifndef VACUUM_MACROS
#define VACUUM_MACROS

#define VAC_REGISTER_SPACE_REFINE_A(field, refiner,refine_op)   \
  refiner.registerRefine(bssnSim->field##_a_idx,                \
                         bssnSim->field##_a_idx,                \
                         bssnSim->field##_a_idx,                \
                         refine_op)
#define VAC_REGISTER_SPACE_REFINE_P(field, refiner,refine_op)   \
  refiner.registerRefine(bssnSim->field##_p_idx,                \
                         bssnSim->field##_p_idx,                \
                         bssnSim->field##_p_idx,                \
                         refine_op)

 

#define VAC_REGISTER_COARSEN_A(field,coarsener,coarsen_op)      \
  coarsener.registerCoarsen(bssnSim->field##_a_idx,             \
                            bssnSim->field##_f_idx,             \
                            coarsen_op,                         \
                            NULL)



#define VAC_REGISTER_COARSEN_P(field,coarsener,coarsen_op)      \
  coarsener.registerCoarsen(bssnSim->field##_p_idx,             \
                            bssnSim->field##_f_idx,             \
                            coarsen_op,                         \
                            NULL)
#endif
