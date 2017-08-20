/* store.h -- store molecule coefficients in Jacobian */
/* $Header$ */

/*
 * prerequisite headers:
 *	"cctk.h"
 *	"../InterpLocalUniform.h"
 *	"structs.h"
 */

void AEILocalInterp_store_1dcube2
      (fp factor, const struct coeffs_struct_1d_cube_size2 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride);
void AEILocalInterp_store_1dcube3
      (fp factor, const struct coeffs_struct_1d_cube_size3 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride);
void AEILocalInterp_store_1dcube4
      (fp factor, const struct coeffs_struct_1d_cube_size4 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride);
void AEILocalInterp_store_1dcube5
      (fp factor, const struct coeffs_struct_1d_cube_size5 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride);
void AEILocalInterp_store_1dcube6
      (fp factor, const struct coeffs_struct_1d_cube_size6 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride);
void AEILocalInterp_store_1dcube7
      (fp factor, const struct coeffs_struct_1d_cube_size7 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride);

void AEILocalInterp_store_2dcube2
      (fp factor, const struct coeffs_struct_2d_cube_size2 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride);
void AEILocalInterp_store_2dcube3
      (fp factor, const struct coeffs_struct_2d_cube_size3 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride);
void AEILocalInterp_store_2dcube4
      (fp factor, const struct coeffs_struct_2d_cube_size4 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride);
void AEILocalInterp_store_2dcube5
      (fp factor, const struct coeffs_struct_2d_cube_size5 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride);
void AEILocalInterp_store_2dcube6
      (fp factor, const struct coeffs_struct_2d_cube_size6 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride);

void AEILocalInterp_store_3dcube2
      (fp factor, const struct coeffs_struct_3d_cube_size2 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride, int Jacobian_mk_stride);
void AEILocalInterp_store_3dcube3
      (fp factor, const struct coeffs_struct_3d_cube_size3 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride, int Jacobian_mk_stride);
void AEILocalInterp_store_3dcube4
      (fp factor, const struct coeffs_struct_3d_cube_size4 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride, int Jacobian_mk_stride);
void AEILocalInterp_store_3dcube5
      (fp factor, const struct coeffs_struct_3d_cube_size5 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride, int Jacobian_mk_stride);
void AEILocalInterp_store_3dcube6
      (fp factor, const struct coeffs_struct_3d_cube_size6 *coeffs,
       fp Jacobian_ptr[],
       int Jacobian_mi_stride, int Jacobian_mj_stride, int Jacobian_mk_stride);
