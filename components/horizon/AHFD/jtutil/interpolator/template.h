/*@@
  @file     template.h
  @date     22 Jan 2002
  @author   Jonathan Thornburg
  @desc     prototype for template function in "template.c"
            (also used for function pointers in "InterpLocalUniform.c")
  @enddesc
  @version  $Header$
  @@*/
int FUNCTION_NAME(/***** coordinate system *****/
		  const CCTK_REAL coord_origin[],
		  const CCTK_REAL coord_delta[],
		  /***** interpolation points *****/
		  int N_interp_points,
		  int interp_coords_type_code,
		  const void* const interp_coords[],
		  const CCTK_INT N_boundary_points_to_omit[],
		  const CCTK_REAL boundary_off_centering_tolerance[],
		  const CCTK_REAL boundary_extrapolation_tolerance[],
		  /***** input arrays *****/
		  int N_input_arrays,
		  const CCTK_INT input_array_offsets[],
		  const CCTK_INT input_array_strides[],
		  const CCTK_INT input_array_min_subscripts[],
		  const CCTK_INT input_array_max_subscripts[],
		  const CCTK_INT input_array_type_codes[],
		  const void* const input_arrays[],
		  /***** output arrays *****/
		  int N_output_arrays,
		  const CCTK_INT output_array_type_codes[],
		  void* const output_arrays[],
		  /***** operation info *****/
		  const CCTK_INT operand_indices[],
		  const CCTK_INT operation_codes[],
		  /***** debugging *****/
		  int debug, FILE* log_fp,
		  /***** other return results *****/
		  struct error_info* error_info,
		  struct molecule_structure_flags* molecule_structure_flags,
		  struct molecule_min_max_m_info* molecule_min_max_m_info,
		  CCTK_INT* const molecule_positions[],
		  struct Jacobian_info* Jacobian_info);
