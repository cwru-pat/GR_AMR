/* evaluate.h -- evaluate a molecule (as a linear combination of the data) */
/* $Header$ */

/*
 * prerequisite headers:
 *	"cctk.h"
 *	"../InterpLocalUniform.h"
 *	"structs.h"
 */

fp AEILocalInterp_eval_1dcube2(const struct coeffs_struct_1d_cube_size2 *coeffs,
			       const struct data_struct_1d_cube_size2 *data);
fp AEILocalInterp_eval_1dcube3(const struct coeffs_struct_1d_cube_size3 *coeffs,
			       const struct data_struct_1d_cube_size3 *data);
fp AEILocalInterp_eval_1dcube4(const struct coeffs_struct_1d_cube_size4 *coeffs,
			       const struct data_struct_1d_cube_size4 *data);
fp AEILocalInterp_eval_1dcube5(const struct coeffs_struct_1d_cube_size5 *coeffs,
			       const struct data_struct_1d_cube_size5 *data);
fp AEILocalInterp_eval_1dcube6(const struct coeffs_struct_1d_cube_size6 *coeffs,
			       const struct data_struct_1d_cube_size6 *data);
fp AEILocalInterp_eval_1dcube7(const struct coeffs_struct_1d_cube_size7 *coeffs,
			       const struct data_struct_1d_cube_size7 *data);

fp AEILocalInterp_eval_2dcube2(const struct coeffs_struct_2d_cube_size2 *coeffs,
			       const struct data_struct_2d_cube_size2 *data);
fp AEILocalInterp_eval_2dcube3(const struct coeffs_struct_2d_cube_size3 *coeffs,
			       const struct data_struct_2d_cube_size3 *data);
fp AEILocalInterp_eval_2dcube4(const struct coeffs_struct_2d_cube_size4 *coeffs,
			       const struct data_struct_2d_cube_size4 *data);
fp AEILocalInterp_eval_2dcube5(const struct coeffs_struct_2d_cube_size5 *coeffs,
			       const struct data_struct_2d_cube_size5 *data);
fp AEILocalInterp_eval_2dcube6(const struct coeffs_struct_2d_cube_size6 *coeffs,
			       const struct data_struct_2d_cube_size6 *data);

fp AEILocalInterp_eval_3dcube2(const struct coeffs_struct_3d_cube_size2 *coeffs,
			       const struct data_struct_3d_cube_size2 *data);
fp AEILocalInterp_eval_3dcube3(const struct coeffs_struct_3d_cube_size3 *coeffs,
			       const struct data_struct_3d_cube_size3 *data);
fp AEILocalInterp_eval_3dcube4(const struct coeffs_struct_3d_cube_size4 *coeffs,
			       const struct data_struct_3d_cube_size4 *data);
fp AEILocalInterp_eval_3dcube5(const struct coeffs_struct_3d_cube_size5 *coeffs,
			       const struct data_struct_3d_cube_size5 *data);
fp AEILocalInterp_eval_3dcube6(const struct coeffs_struct_3d_cube_size6 *coeffs,
			       const struct data_struct_3d_cube_size6 *data);
