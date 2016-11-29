#define VAC_REGISTER_SPACE_REFINE_A(field, refiner,refine_op)  \
  refiner##.registerRefine(bssnSim->field##_a_idx,  \
                           bssnSim->field##_a_idx,  \
                           bssnSim->field##_a_idx,  \
                           refine_op)
#define VAC_REGISTER_SPACE_REFINE_P(field, refiner,refine_op)  \
  refiner##.registerRefine(bssnSim->field##_p_idx,  \
                           bssnSim->field##_p_idx,  \
                           bssnSim->field##_p_idx,  \
                           refine_op)

#define VAC_REGISTER_SPACE_REFINE_C(field, refiner,refine_op)  \
  refiner##.registerRefine(bssnSim->field##_c_idx,  \
                           bssnSim->field##_c_idx,  \
                           bssnSim->field##_c_idx,  \
                           refine_op)

#define VAC_REGISTER_SPACE_REFINE_F(field, refiner,refine_op)  \
  refiner##.registerRefine(bssnSim->field##_f_idx,  \
                           bssnSim->field##_f_idx,  \
                           bssnSim->field##_f_idx,  \
                           refine_op)
 

#define VAC_REGISTER_INTER_LEVEL_REFINE_A(field, refiner,space_refine_op, time_refine_op) \
  refiner##.registerRefine(bssnSim->field##_a_idx,  \
                           bssnSim->field##_a_idx,  \
                           bssnSim->field##_p_idx,  \
                           bssnSim->field##_f_idx,  \
                           bssnSim->field##_a_idx,  \
                           space_refine_op,  \
                           time_refine_op)
#define VAC_REGISTER_INTER_LEVEL_REFINE_P(field, refiner,space_refine_op, time_refine_op) \
  refiner##.registerRefine(bssnSim->field##_p_idx,  \
                           bssnSim->field##_p_idx,  \
                           bssnSim->field##_p_idx,  \
                           bssnSim->field##_f_idx,  \
                           bssnSim->field##_p_idx,  \
                           space_refine_op,  \
                           time_refine_op)
#define VAC_REGISTER_INTER_LEVEL_REFINE_C(field, refiner,space_refine_op, time_refine_op) \
  refiner##.registerRefine(bssnSim->field##_c_idx,  \
                           bssnSim->field##_c_idx,  \
                           bssnSim->field##_p_idx,  \
                           bssnSim->field##_f_idx,  \
                           bssnSim->field##_c_idx,  \
                           space_refine_op,  \
                           time_refine_op)
#define VAC_REGISTER_INTER_LEVEL_REFINE_F(field, refiner,space_refine_op, time_refine_op) \
  refiner##.registerRefine(bssnSim->field##_f_idx,  \
                           bssnSim->field##_f_idx,  \
                           bssnSim->field##_p_idx,  \
                           bssnSim->field##_f_idx,  \
                           bssnSim->field##_f_idx,  \
                           space_refine_op,  \
                           time_refine_op)

#define VAC_REGISTER_COARSEN_A(field,coarsener,coarsen_op)  \
  coarsener##.registerCoarsen(bssnSim->field##_a_idx,  \
                              bssnSim->field##_f_idx,  \
                              coarsen_op,     \
                              NULL)

#define VAC_REGISTER_COARSEN_C(field,coarsener,coarsen_op)     \
  coarsener##.registerCoarsen(bssnSim->field##_c_idx,  \
                              bssnSim->field##_f_idx,  \
                              coarsen_op,     \
                              NULL)

#define VAC_REGISTER_COARSEN_F(field,coarsener,coarsen_op)     \
  coarsener##.registerCoarsen(bssnSim->field##_f_idx,  \
                              bssnSim->field##_f_idx,  \
                              coarsen_op,     \
                              NULL)

#define VAC_REGISTER_COARSEN_P(field,coarsener,coarsen_op)     \
  coarsener##.registerCoarsen(bssnSim->field##_p_idx,  \
                              bssnSim->field##_f_idx,  \
                              coarsen_op,     \
                              NULL)

