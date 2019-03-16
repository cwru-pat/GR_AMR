#include "AHFD.h"

using namespace SAMRAI;
using namespace cosmo;

namespace AHFinderDirect
{

Horizon::Horizon(const std::shared_ptr<hier::PatchHierarchy>& hierarchy_in,
                 BSSN * bssn_in,
                 const tbox::Dimension& dim_in,
                 std::shared_ptr<tbox::Database> database_in,
                 const char * visit_d_name,
                 int w_idx_in):
  hierarchy(hierarchy_in),
  bssn(bssn_in),
  AHFD_db(database_in),
  dim(dim_in),
  w_idx(w_idx_in),
  ILUCG__error_tolerance(1e-10),
  ILUCG__limit_CG_iterations(true),
  UMFPACK__N_II_iterations(0),
  Jacobian_perturbation_amplitude(1e-6),
  test_all_Jacobian_compute_methods(true),
  coordinate_system_name("cart3d"),
  geometry_interpolator_name("Hermite polynomial interpolation"),
  geometry_interpolator_pars("order=3 boundary_off_centering_tolerance={1.0e-10 1.0e-10 1.0e-10 1.0e-10 1.0e-10 1.0e-10} boundary_extrapolation_tolerance={0.0 0.0 0.0 0.0 0.0 0.0}"),
  hardwire_Schwarzschild_EF_geometry(false),
  geometry__Schwarzschild_EF__mass(1.0),
  geometry__Schwarzschild_EF__x_posn(0.0),
  geometry__Schwarzschild_EF__y_posn(0.0),
  geometry__Schwarzschild_EF__z_posn(0.0),
  geometry__Schwarzschild_EF__epsilon(1e-9),
  geometry__Schwarzschild_EF__Delta_xyz(1e-6),
  check_that_h_is_finite(true),
  check_that_geometry_is_finite(true),
  interpatch_interpolator_name("Lagrange polynomial interpolation"),
  interpatch_interpolator_pars("order=5"),
  surface_interpolator_name("Lagrange polynomial interpolation"),
  surface_interpolator_pars("order=2 boundary_off_centering_tolerance={1.0e-10 1.0e-10 1.0e-10 1.0e-10} boundary_extrapolation_tolerance={0.0 0.0 0.0 0.0}"),
  integral_method("automatic choice"),
  find_every(AHFD_db->getIntegerWithDefault("find_every",1)),
  find_every_individual{0,-1},
  find_after_individual{0,0},
  dont_find_after_individual{0,-1},
  find_after_individual_time{0,0},
  dont_find_after_individual_time{0,0.0},
  disable_horizon{false,false},
  method("find horizons"),
  N_horizons(AHFD_db->getIntegerWithDefault("N_horizons", 1)),
  which_horizon_to_announce_centroid(0),
  which_surface_to_store_info{0,-1},
  which_surface_to_store_info_by_name{"",""},
  verbose_level("algorithm highlights"),
  print_timing_stats(false),
  want_expansion_gradients(false),
  surface_definition{"","expansion"},
  surface_modification{"","none"},
  surface_selection{"","definition"},
  desired_value{0,0.0},
  use_pretracking{0,0},
  pretracking_value{0,1.0},
  pretracking_minimum_value{0,0.0},
  pretracking_maximum_value{0,10.0},
  pretracking_delta{0,1.0},
  pretracking_minimum_delta{0,1e-4},
  pretracking_maximum_delta{0,1.0},
  depends_on{0,0},
  desired_value_factor{0,1.0},
  desired_value_offset{0,0.0},
  shiftout_factor{0,1.0},
  smoothing_factor{0,0.0},
  initial_guess_method{"","coordinate sphere"},
  reset_horizon_after_not_finding{0,1},
  initial_guess__read_from_named_file__file_name{"","h.gp"},
  initial_guess__Kerr_Kerr__x_posn{0,0.0},
  initial_guess__Kerr_Kerr__y_posn{0,0.0},
  initial_guess__Kerr_Kerr__z_posn{0,0.0},
  initial_guess__Kerr_Kerr__mass{0,0},
  initial_guess__Kerr_Kerr__spin{0,0.6},
  initial_guess__Kerr_KerrSchild__x_posn{0,0.0},
  initial_guess__Kerr_KerrSchild__y_posn{0,0.0},
  initial_guess__Kerr_KerrSchild__z_posn{0,0.0},
  initial_guess__Kerr_KerrSchild__mass{0,1.0},
  initial_guess__Kerr_KerrSchild__spin{0,0.6},
  initial_guess__coord_sphere__x_center{0,0.0},
  initial_guess__coord_sphere__y_center{0,0.0},
  initial_guess__coord_sphere__z_center{0,0.0},
  initial_guess__coord_sphere__radius{0,2.0},
  initial_guess__coord_ellipsoid__x_center{0,0.0},
  initial_guess__coord_ellipsoid__y_center{0,0.0},
  initial_guess__coord_ellipsoid__z_center{0,0.0},
  initial_guess__coord_ellipsoid__x_radius{0,2.0},
  initial_guess__coord_ellipsoid__y_radius{0,2.0},
  initial_guess__coord_ellipsoid__z_radius{0,2.0},
  output_BH_diagnostics(true),
  BH_diagnostics_directory(""),
  BH_diagnostics_base_file_name("BH_diagnostics"),
  BH_diagnostics_file_name_extension("gp"),
  output_h_every(1),
  output_Theta_every(0),
  output_mean_curvature_every(0),
  output_ASCII_files(1),
  output_HDF5_files(0),
  output_ghost_zones_for_h(0),
  ASCII_gnuplot_file_name_extension("gp"),
  HDF5_file_name_extension("h5"),
  h_directory(""),
  h_base_file_name("h"),
  Theta_base_file_name("Theta"),
  mean_curvature_base_file_name("mean_curvature"),
  Delta_h_base_file_name("Delta_h"),
  h_min_digits(0),
  output_OpenDX_control_files(1),
  OpenDX_control_file_name_extension("dx"),
  output_initial_guess(0),
  debugging_output_at_each_Newton_iteration(0),
  Jacobian_base_file_name("Jacobian.dat"),
  set_mask_for_all_horizons(0),
  set_mask_for_individual_horizon{0,0},
  mask_radius_multiplier(0.8),
  mask_radius_offset(-5.0),
  mask_buffer_thickness(5.0),
  mask_is_noshrink(1),
  min_horizon_radius_points_for_mask(-1e10),
  set_old_style_mask(1),
  set_new_style_mask(0),
  old_style_mask_gridfn_name("SpaceMask::emask"),
  old_style_mask_inside_value(0.0),
  old_style_mask_buffer_value(0.5),
  old_style_mask_outside_value(1.0),
  new_style_mask_gridfn_name("SpaceMask::space_mask"),
  new_style_mask_bitfield_name("mask"),
  new_style_mask_inside_value("inside"),
  new_style_mask_buffer_value("buffer"),
  new_style_mask_outside_value("outside"),
  warn_level__point_outside__initial(1),
  warn_level__point_outside__subsequent(2),
  warn_level__skipping_finite_check(3),
  warn_level__nonfinite_geometry(1),
  warn_level__gij_not_positive_definite__initial(2),
  warn_level__gij_not_positive_definite__subsequent(2),
  max_Newton_iterations__initial(AHFD_db->getIntegerWithDefault("max_Newton_iterations_initial",20)),
  max_Newton_iterations__subsequent(AHFD_db->getIntegerWithDefault("max_Newton_iterations_subsequent",10)),
  max_allowable_Delta_h_over_h(0.1),
  max_allowable_horizon_radius{0.0,1e10},
  Theta_norm_for_convergence(AHFD_db->getDoubleWithDefault("Theta_norm_for_convergence",1e-7)),
  max_allowable_Theta(1e10),
  max_allowable_Theta_growth_iterations(0),
  max_allowable_Theta_nonshrink_iterations(0),
  origin_x{0,0},
  origin_y{0,0},
  origin_z{0,0},
  move_origins(0),
  reshape_while_moving(0),
  predict_origin_movement(0),
  track_origin_from_grid_scalar{0,0},
  track_origin_source_x{"",""},
  track_origin_source_y{"",""},
  track_origin_source_z{"",""},
  patch_system_type{"","full sphere"},
  N_zones_per_right_angle{0,18},
  max_N_zones_per_right_angle(18),
  ghost_zone_width(2),
  patch_overlap_width(1),
  Jacobian_compute_method("symbolic differentiation with finite diff d/dr"),
  Jacobian_store_solve_method("row-oriented sparse matrix/ILUCG")
{
  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));

  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  const double * lower = &grid_geometry.getXLower()[0];
  const double * upper = &grid_geometry.getXUpper()[0];

  
  for(int i = 0 ; i < DIM; i++)
  {
    domain_lower[i] = lower[i];
    domain_upper[i] = upper[i];
  }

  AHFD_VEC_INIT(find_every_individual, -1);
  AHFD_VEC_INIT(find_after_individual, 0);
  AHFD_VEC_INIT(dont_find_after_individual,-1);
  AHFD_VEC_INIT(find_after_individual_time,0);
  AHFD_VEC_INIT(dont_find_after_individual_time,0);
  AHFD_VEC_INIT(disable_horizon,false);
  AHFD_VEC_INIT(which_surface_to_store_info,-1);
  AHFD_VEC_INIT(which_surface_to_store_info_by_name,"");
  AHFD_VEC_INIT(surface_definition,"expansion");
  AHFD_VEC_INIT(surface_modification,"none");
  AHFD_VEC_INIT(surface_selection,"definition");
  AHFD_VEC_INIT(desired_value,0);
  AHFD_VEC_INIT(use_pretracking,0);
  AHFD_VEC_INIT(pretracking_value,1.0);
  AHFD_VEC_INIT(pretracking_minimum_value,0);
  AHFD_VEC_INIT(pretracking_maximum_value,10.0);
  AHFD_VEC_INIT(pretracking_delta,1.0);
  AHFD_VEC_INIT(pretracking_minimum_delta,1e-4);
  AHFD_VEC_INIT(pretracking_maximum_delta,1.0);
  AHFD_VEC_INIT(depends_on,0);
  AHFD_VEC_INIT(desired_value_factor,1.0);
  AHFD_VEC_INIT(desired_value_offset,0);
  AHFD_VEC_INIT(shiftout_factor,1.0);
  AHFD_VEC_INIT(smoothing_factor,0.0);
  AHFD_VEC_INIT(reset_horizon_after_not_finding,1);
  AHFD_VEC_INIT(initial_guess__read_from_named_file__file_name,"h.gp");
  AHFD_VEC_INIT(initial_guess__coord_sphere__x_center,0);
  AHFD_VEC_INIT(initial_guess__coord_sphere__y_center,0);
  AHFD_VEC_INIT(initial_guess__coord_sphere__z_center,0);
  AHFD_VEC_INIT(initial_guess__coord_sphere__radius,2.0);
  AHFD_VEC_INIT(set_mask_for_individual_horizon,0);
  AHFD_VEC_INIT(max_allowable_horizon_radius, (domain_upper[0]-domain_lower[0]) / 2.0);
  AHFD_VEC_INIT(origin_x,0);
  AHFD_VEC_INIT(origin_y,0);
  AHFD_VEC_INIT(origin_z,0);
  AHFD_VEC_INIT(track_origin_from_grid_scalar,0);
  AHFD_VEC_INIT(track_origin_source_x,"");
  AHFD_VEC_INIT(track_origin_source_y,"");
  AHFD_VEC_INIT(track_origin_source_z,"");
  AHFD_VEC_INIT(patch_system_type,"full sphere");
  AHFD_VEC_INIT(N_zones_per_right_angle,18);

  
  
  getcwd(cur_directory, FILENAME_MAX);

  strcat(cur_directory,"/");
  strcat(cur_directory,visit_d_name);
  strcat(cur_directory,".visit");
  
  
  /*********initializing initial guess******************************/

  static std::string initial_guess_type;
  initial_guess_type = AHFD_db->getStringWithDefault("initial_guess_type", "coordinate sphere");

  AHFD_VEC_INIT(initial_guess_method,initial_guess_type.c_str());

  if(initial_guess_type == "coordinate sphere")
  {
    AHFD_db->getDoubleArray("origin_x", origin_x+1, N_horizons);
    AHFD_db->getDoubleArray("origin_y", origin_y+1, N_horizons);
    AHFD_db->getDoubleArray("origin_z", origin_z+1, N_horizons);

    AHFD_db->getDoubleArray(
      "sphere_x_center", initial_guess__coord_sphere__x_center+1, N_horizons);
    AHFD_db->getDoubleArray(
      "sphere_y_center", initial_guess__coord_sphere__y_center+1, N_horizons);
    AHFD_db->getDoubleArray(
      "sphere_z_center", initial_guess__coord_sphere__z_center+1, N_horizons);

    AHFD_db->getDoubleArray(
      "sphere_radius",initial_guess__coord_sphere__radius+1, N_horizons);
  }
  else if(initial_guess_type == "coordinate ellipsoid")
  {
    AHFD_db->getDoubleArray("origin_x", origin_x+1, N_horizons);
    AHFD_db->getDoubleArray("origin_y", origin_y+1, N_horizons);
    AHFD_db->getDoubleArray("origin_z", origin_z+1, N_horizons);

    AHFD_db->getDoubleArray(
      "ellipsoid_x_center", initial_guess__coord_ellipsoid__x_center+1, N_horizons);
    AHFD_db->getDoubleArray(
      "ellipsoid_y_center", initial_guess__coord_ellipsoid__y_center+1, N_horizons);
    AHFD_db->getDoubleArray(
      "ellipsoid_z_center", initial_guess__coord_ellipsoid__z_center+1, N_horizons);

    AHFD_db->getDoubleArray(
      "ellipsoid_x_radius",initial_guess__coord_ellipsoid__x_radius+1, N_horizons);
    AHFD_db->getDoubleArray(
      "ellipsoid_y_radius",initial_guess__coord_ellipsoid__y_radius+1, N_horizons);
    AHFD_db->getDoubleArray(
      "ellipsoid_z_radius",initial_guess__coord_ellipsoid__z_radius+1, N_horizons);

  }
  else
  {
    TBOX_ERROR("Unrecgnized horizon initial guess type!\n");
  }
  /*********initializing other settings******************************/
  if(AHFD_db->keyExists("find_after_individual"))
    AHFD_db->getIntegerArray(
      "find_after_individual",find_after_individual+1, N_horizons);
  
  if(AHFD_db->keyExists("max_allowable_horizon_radius"))
    AHFD_db->getDoubleArray(
      "max_allowable_horizon_radius",max_allowable_horizon_radius+1, N_horizons);

  if(AHFD_db->keyExists("N_zones_per_right_angle"))
    AHFD_db->getIntegerArray(
      "N_zones_per_right_angle",N_zones_per_right_angle+1, N_horizons);
  
  
}

int Horizon::allocate_horizons_to_processor(int N_procs, int my_proc,
                                   int N_horizons, bool multiproc_flag,
                                   const CCTK_INT depends_on[],
                                   horizon_sequence& my_hs,
                                   const struct verbose_info& verbose_info)
{
  const int N_active_procs = multiproc_flag ? jtutil::min(N_procs, N_horizons)
    : 1;

  //
  // Implementation note:
  // We allocate the horizons to active processors in round-robin order.
  //
  std::vector<int> proc_of_horizon (N_horizons+1);
  for (int hn = 1 ; hn <= N_horizons ; ++hn)
  {
    proc_of_horizon.at(hn) = -1;
  }

  int proc = 0;
  for (int hn = 1 ; hn <= N_horizons ; ++hn)
  {
    if (depends_on[hn] < 0 || depends_on[hn] > N_horizons) {
      CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "horizon %d depends on a horizon with the illegal index %d",
                 hn, int(depends_on[hn]));
    } else if (depends_on[hn] == hn) {
      CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "horizon %d depends on itself",
                 hn);
    } else if (depends_on[hn] > hn) {
      CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "horizon %d depends on a horizon with a larger index %d",
                 hn, int(depends_on[hn]));
    }
    const int this_horizons_proc
      = depends_on[hn] == 0 ? proc : proc_of_horizon.at(depends_on[hn]);
    assert (this_horizons_proc >= 0 && this_horizons_proc < N_procs);
    proc_of_horizon.at(hn) = this_horizons_proc;
    if (verbose_info.print_algorithm_highlights)
      then CCTK_VInfo(CCTK_THORNSTRING,
                      "   allocating horizon %d to processor #%d",
                      hn, this_horizons_proc);
    if (this_horizons_proc == my_proc)
      then my_hs.append_hn(hn);
    if (++proc >= N_active_procs)
      then proc = 0;
  }

  return N_active_procs;
}

enum verbose_level
Horizon::decode_verbose_level(const char verbose_level_string[])
{
if	(STRING_EQUAL(verbose_level_string, "physics highlights"))
   then return verbose_level__physics_highlights;
else if (STRING_EQUAL(verbose_level_string, "physics details"))
   then return verbose_level__physics_details;
else if (STRING_EQUAL(verbose_level_string, "algorithm highlights"))
   then return verbose_level__algorithm_highlights;
else if (STRING_EQUAL(verbose_level_string, "algorithm details"))
   then return verbose_level__algorithm_details;
else if (STRING_EQUAL(verbose_level_string, "algorithm debug"))
   then return verbose_level__algorithm_debug;
else	CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
"decode_verbose_level(): unknown verbose_level_string=\"%s\"!",
		   verbose_level_string);			/*NOTREACHED*/
}

  
enum method
Horizon::decode_method(const char method_string[])
{
if	(STRING_EQUAL(method_string, "evaluate expansions"))
   then return method__evaluate_expansions;
else if (STRING_EQUAL(method_string, "test expansion Jacobians"))
   then return method__test_expansion_Jacobians;
else if (STRING_EQUAL(method_string, "find horizons"))
   then return method__find_horizons;
else	CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
		    "decode_method(): unknown method_string=\"%s\"!",
		    method_string);				/*NOTREACHED*/
}

enum initial_guess_method
Horizon::decode_initial_guess_method(const char initial_guess_method_string[])
{
if	(STRING_EQUAL(initial_guess_method_string, "read from named file"))
   then return initial_guess__read_from_named_file;
else if (STRING_EQUAL(initial_guess_method_string, "read from h file"))
   then return initial_guess__read_from_h_file;
else if (STRING_EQUAL(initial_guess_method_string, "Kerr/Kerr"))
   then return initial_guess__Kerr_Kerr;
else if (STRING_EQUAL(initial_guess_method_string, "Kerr/Kerr-Schild"))
   then return initial_guess__Kerr_KerrSchild;
else if (STRING_EQUAL(initial_guess_method_string, "coordinate sphere"))
   then return initial_guess__coord_sphere;
else if (STRING_EQUAL(initial_guess_method_string, "coordinate ellipsoid"))
   then return initial_guess__coord_ellipsoid;
else	CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   decode_initial_guess_method():\n"
"        unknown initial_guess_method_string=\"%s\"!",
		   initial_guess_method_string);		/*NOTREACHED*/
}

  
void Horizon::set_initial_guess_parameters(struct AH_data& AH_data, const int hn,
                                  const fp ini_origin_x, const fp ini_origin_y, const fp ini_origin_z)
{  
  AH_data.initial_guess_info.method
    = decode_initial_guess_method(initial_guess_method[hn]);
  AH_data.initial_guess_info.reset_horizon_after_not_finding
    = reset_horizon_after_not_finding[hn];
  // ... read from named file
  AH_data.initial_guess_info.read_from_named_file_info.file_name
    = initial_guess__read_from_named_file__file_name[hn];

  if (!track_origin_from_grid_scalar[hn]) {

    // ... Kerr/Kerr
    AH_data.initial_guess_info.Kerr_Kerr_info.x_posn
      = initial_guess__Kerr_Kerr__x_posn[hn];
    AH_data.initial_guess_info.Kerr_Kerr_info.y_posn
      = initial_guess__Kerr_Kerr__y_posn[hn];
    AH_data.initial_guess_info.Kerr_Kerr_info.z_posn
      = initial_guess__Kerr_Kerr__z_posn[hn];
    AH_data.initial_guess_info.Kerr_Kerr_info.mass
      = initial_guess__Kerr_Kerr__mass[hn];
    AH_data.initial_guess_info.Kerr_Kerr_info.spin
      = initial_guess__Kerr_Kerr__spin[hn];
    // ... Kerr/Kerr-Schild
    AH_data.initial_guess_info.Kerr_KerrSchild_info.x_posn
      = initial_guess__Kerr_KerrSchild__x_posn[hn];
    AH_data.initial_guess_info.Kerr_KerrSchild_info.y_posn
      = initial_guess__Kerr_KerrSchild__y_posn[hn];
    AH_data.initial_guess_info.Kerr_KerrSchild_info.z_posn
      = initial_guess__Kerr_KerrSchild__z_posn[hn];
    AH_data.initial_guess_info.Kerr_KerrSchild_info.mass
      = initial_guess__Kerr_KerrSchild__mass[hn];
    AH_data.initial_guess_info.Kerr_KerrSchild_info.spin
      = initial_guess__Kerr_KerrSchild__spin[hn];
    // ... coordinate sphere
    AH_data.initial_guess_info.coord_sphere_info.x_center
      = initial_guess__coord_sphere__x_center[hn];
    AH_data.initial_guess_info.coord_sphere_info.y_center
      = initial_guess__coord_sphere__y_center[hn];
    AH_data.initial_guess_info.coord_sphere_info.z_center
      = initial_guess__coord_sphere__z_center[hn];
    AH_data.initial_guess_info.coord_sphere_info.radius
      = initial_guess__coord_sphere__radius[hn];
    // ... coordinate ellipsoid
    AH_data.initial_guess_info.coord_ellipsoid_info.x_center
      = initial_guess__coord_ellipsoid__x_center[hn];
    AH_data.initial_guess_info.coord_ellipsoid_info.y_center
      = initial_guess__coord_ellipsoid__y_center[hn];
    AH_data.initial_guess_info.coord_ellipsoid_info.z_center
      = initial_guess__coord_ellipsoid__z_center[hn];
    AH_data.initial_guess_info.coord_ellipsoid_info.x_radius
      = initial_guess__coord_ellipsoid__x_radius[hn];
    AH_data.initial_guess_info.coord_ellipsoid_info.y_radius
      = initial_guess__coord_ellipsoid__y_radius[hn];
    AH_data.initial_guess_info.coord_ellipsoid_info.z_radius
      = initial_guess__coord_ellipsoid__z_radius[hn];

  } else {

    // ... Kerr/Kerr
    AH_data.initial_guess_info.Kerr_Kerr_info.x_posn
      = ini_origin_x;
    AH_data.initial_guess_info.Kerr_Kerr_info.y_posn
      = ini_origin_y;
    AH_data.initial_guess_info.Kerr_Kerr_info.z_posn
      = ini_origin_z;
    AH_data.initial_guess_info.Kerr_Kerr_info.mass
      = initial_guess__Kerr_Kerr__mass[hn];
    AH_data.initial_guess_info.Kerr_Kerr_info.spin
      = initial_guess__Kerr_Kerr__spin[hn];
    // ... Kerr/Kerr-Schild
    AH_data.initial_guess_info.Kerr_KerrSchild_info.x_posn
      = ini_origin_x;
    AH_data.initial_guess_info.Kerr_KerrSchild_info.y_posn
      = ini_origin_y;
    AH_data.initial_guess_info.Kerr_KerrSchild_info.z_posn
      = ini_origin_z;
    AH_data.initial_guess_info.Kerr_KerrSchild_info.mass
      = initial_guess__Kerr_KerrSchild__mass[hn];
    AH_data.initial_guess_info.Kerr_KerrSchild_info.spin
      = initial_guess__Kerr_KerrSchild__spin[hn];
    // ... coordinate sphere
    AH_data.initial_guess_info.coord_sphere_info.x_center
      = ini_origin_x;
    AH_data.initial_guess_info.coord_sphere_info.y_center
      = ini_origin_y;
    AH_data.initial_guess_info.coord_sphere_info.z_center
      = ini_origin_z;
    AH_data.initial_guess_info.coord_sphere_info.radius
      = initial_guess__coord_sphere__radius[hn];
    // ... coordinate ellipsoid
    AH_data.initial_guess_info.coord_ellipsoid_info.x_center
      = ini_origin_x;
    AH_data.initial_guess_info.coord_ellipsoid_info.y_center
      = ini_origin_y;
    AH_data.initial_guess_info.coord_ellipsoid_info.z_center
      = ini_origin_z;
    AH_data.initial_guess_info.coord_ellipsoid_info.x_radius
      = initial_guess__coord_ellipsoid__x_radius[hn];
    AH_data.initial_guess_info.coord_ellipsoid_info.y_radius
      = initial_guess__coord_ellipsoid__y_radius[hn];
    AH_data.initial_guess_info.coord_ellipsoid_info.z_radius
      = initial_guess__coord_ellipsoid__z_radius[hn];

  }

}
  
void Horizon::AHFinderDirect_setup()
{
  int need_zones = 0;
  for (int n=1; n<N_horizons; ++n) {
    need_zones = jtutil::max(need_zones, int(N_zones_per_right_angle[n]));
  }
  if (need_zones > max_N_zones_per_right_angle) {
    CCTK_VWarn (FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                "AHFinderDirect_setup(): "
                "The parameter max_N_zones_per_right_angle must be at least the maximum of all N_zones_per_right_angle[] parameters.  "
                "Set max_N_zones_per_right_angle to %d or higher to continue.",
                need_zones);      /*NOTREACHED*/
  }
  for (int n=1; n<N_horizons; ++n) {
    if (depends_on[n] != 0) {
      assert (depends_on[n] >= 1 && depends_on[n] < n);
      if (N_zones_per_right_angle[n] != N_zones_per_right_angle[depends_on[n]]) {
        CCTK_VWarn (FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "AHFinderDirect_setup(): "
                    "The parameter N_zones_per_right_angle must be the same for a horizon as for the horizon on which it depends.  "
                    "Horizon %d depends on horizon %d, but they have different resolutions.",
                    n, int(depends_on[n])); /*NOTREACHED*/
      }
    }
  }

  bool find_individual_is_set = false;
  for (int n = 1 ; n < N_horizons ; ++n)
  {
    if (find_every_individual[n] > 0)
      then find_individual_is_set = true;
  }
  if (move_origins && find_individual_is_set)
    then {
      CCTK_VWarn (CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Find_every_individual is currently not compatible with moving origins via the move_origins parameter.  "
                  "Make sure simultaneous moving horizons are not being calculated at different frequencies." );
    }

  //
  // basic setup
  //
  state.method = decode_method(method);

  state.error_info.warn_level__point_outside__initial
    = warn_level__point_outside__initial;
  state.error_info.warn_level__point_outside__subsequent
    = warn_level__point_outside__subsequent;
  state.error_info.warn_level__skipping_finite_check
    = warn_level__skipping_finite_check;
  state.error_info.warn_level__nonfinite_geometry
    = warn_level__nonfinite_geometry;
  state.error_info.warn_level__gij_not_positive_definite__initial
    = warn_level__gij_not_positive_definite__initial;
  state.error_info.warn_level__gij_not_positive_definite__subsequent
    = warn_level__gij_not_positive_definite__subsequent;

  struct verbose_info& verbose_info = state.verbose_info;
  verbose_info.verbose_level = decode_verbose_level(verbose_level);
  verbose_info.print_physics_highlights
    = (state.verbose_info.verbose_level >= verbose_level__physics_highlights);
  verbose_info.print_physics_details
    = (state.verbose_info.verbose_level >= verbose_level__physics_details);
  verbose_info.print_algorithm_highlights
    = (state.verbose_info.verbose_level >= verbose_level__algorithm_highlights);
  verbose_info.print_algorithm_details
    = (state.verbose_info.verbose_level >= verbose_level__algorithm_details);
  verbose_info.print_algorithm_debug
    = (state.verbose_info.verbose_level >= verbose_level__algorithm_debug);

  // TODOMARKS
  // set timer
  // state.timer_handle = (print_timing_stats != 0) ? CCTK_TimerCreate("finding apparent horizons") : -1;

  state.N_procs = hierarchy->getMPI().getSize();
  state.my_proc = hierarchy->getMPI().getRank();
  
  state.N_horizons = N_horizons;
  state.N_active_procs = 0;	// dummy value, will be set properly later
  CCTK_VInfo(CCTK_THORNSTRING,
             "           to search for %d horizon%s on %d processor%s",
             state.N_horizons, ((state.N_horizons == 1) ? "" : "s"),
             state.N_procs, ((state.N_procs == 1) ? "" : "s"));


//
// Cactus grid info
//
 if (verbose_info.print_algorithm_highlights)
    then CCTK_VInfo(CCTK_THORNSTRING, "   setting up Cactus grid info");
 //struct cactus_grid_info& cgi = state.cgi;
  //cgi.GH = cctkGH;
  //cgi.coord_system_handle = CCTK_CoordSystemHandle(coordinate_system_name);
  //cgi.use_Cactus_conformal_metric = false;	// dummy value, may change later

  // TODOMARKS
  // May not need this
  // cgi.mask_varindex    = Cactus_gridfn_varindex("AHFinderDirect::ahmask");
  // cgi.g_dd_11_varindex = Cactus_gridfn_varindex("ADMBase::gxx");
  // cgi.g_dd_12_varindex = Cactus_gridfn_varindex("ADMBase::gxy");
  // cgi.g_dd_13_varindex = Cactus_gridfn_varindex("ADMBase::gxz");
  // cgi.g_dd_22_varindex = Cactus_gridfn_varindex("ADMBase::gyy");
  // cgi.g_dd_23_varindex = Cactus_gridfn_varindex("ADMBase::gyz");
  // cgi.g_dd_33_varindex = Cactus_gridfn_varindex("ADMBase::gzz");
  // cgi.K_dd_11_varindex = Cactus_gridfn_varindex("ADMBase::kxx");
  // cgi.K_dd_12_varindex = Cactus_gridfn_varindex("ADMBase::kxy");
  // cgi.K_dd_13_varindex = Cactus_gridfn_varindex("ADMBase::kxz");
  // cgi.K_dd_22_varindex = Cactus_gridfn_varindex("ADMBase::kyy");
  // cgi.K_dd_23_varindex = Cactus_gridfn_varindex("ADMBase::kyz");
  // cgi.K_dd_33_varindex = Cactus_gridfn_varindex("ADMBase::kzz");
  // cgi.psi_varindex     = Cactus_gridfn_varindex("StaticConformal::psi");


  //
  // geometry info
  //
 // TODOMARKS
  if (verbose_info.print_algorithm_highlights)
    then CCTK_VInfo(CCTK_THORNSTRING, "   setting up geometry interpolator");
  struct geometry_info& gi = state.gi;
  // gi.hardwire_Schwarzschild_EF_geometry
  //   = (hardwire_Schwarzschild_EF_geometry != 0);

  // gi.operator_handle = CCTK_InterpHandle(geometry_interpolator_name);
  // if (gi.operator_handle < 0)
  //   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
  //                   "AHFinderDirect_setup(): couldn't find interpolator \"%s\"!",
  //                   geometry_interpolator_name);		/*NOTREACHED*/
  gi.param_table_handle = Util_TableCreateFromString(geometry_interpolator_pars);
  if (gi.param_table_handle < 0)
    then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "AHFinderDirect_setup(): bad geometry-interpolator parameter(s) \"%s\"!",
                    geometry_interpolator_pars);		/*NOTREACHED*/

  // gi.geometry__Schwarzschild_EF__mass     = geometry__Schwarzschild_EF__mass;
  // gi.geometry__Schwarzschild_EF__x_posn   = geometry__Schwarzschild_EF__x_posn;
  // gi.geometry__Schwarzschild_EF__y_posn   = geometry__Schwarzschild_EF__y_posn;
  // gi.geometry__Schwarzschild_EF__z_posn   = geometry__Schwarzschild_EF__z_posn;
  // gi.geometry__Schwarzschild_EF__epsilon  = geometry__Schwarzschild_EF__epsilon;
  // gi.geometry__Schwarzschild_EF__Delta_xyz= geometry__Schwarzschild_EF__Delta_xyz;

  // gi.check_that_h_is_finite        = (check_that_h_is_finite        != 0);
  // gi.check_that_geometry_is_finite = (check_that_geometry_is_finite != 0);
  // gi.mask_is_noshrink		 = (mask_is_noshrink != 0);

  //
  // Jacobian info
  //
  struct Jacobian_info& Jac_info = state.Jac_info;
  Jac_info.Jacobian_compute_method
    = decode_Jacobian_compute_method(Jacobian_compute_method);
  Jac_info.Jacobian_store_solve_method
    = decode_Jacobian_store_solve_method(Jacobian_store_solve_method);
  Jac_info.perturbation_amplitude = Jacobian_perturbation_amplitude;


  //
  // solver info
  //
  struct solver_info& solver_info = state.solver_info;
  solver_info.debugging_output_at_each_Newton_iteration
    = (debugging_output_at_each_Newton_iteration != 0);
  solver_info.linear_solver_pars.ILUCG_pars.error_tolerance
    = ILUCG__error_tolerance;
  solver_info.linear_solver_pars.ILUCG_pars.limit_CG_iterations
    = (ILUCG__limit_CG_iterations != 0);
  solver_info.linear_solver_pars.UMFPACK_pars.N_II_iterations
    = UMFPACK__N_II_iterations;
  solver_info.max_Newton_iterations__initial
    = max_Newton_iterations__initial;
  solver_info.max_Newton_iterations__subsequent
    = max_Newton_iterations__subsequent;
  solver_info.max_allowable_Delta_h_over_h = max_allowable_Delta_h_over_h;
  solver_info.Theta_norm_for_convergence   = Theta_norm_for_convergence;
  solver_info.max_allowable_Theta          = max_allowable_Theta;
  solver_info.max_allowable_Theta_growth_iterations
    = max_allowable_Theta_growth_iterations;
  solver_info.max_allowable_Theta_nonshrink_iterations
    = max_allowable_Theta_nonshrink_iterations;
  // ... horizon numbers run from 1 to N_horizons inclusive
  //     so the array size is N_horizons+1
  solver_info.max_allowable_horizon_radius = new fp[state.N_horizons+1];
 {
   for (int hn = 0 ; hn <= N_horizons ; ++hn)
   {
     solver_info.max_allowable_horizon_radius[hn]
       = max_allowable_horizon_radius[hn];
   }
 }
 solver_info.want_expansion_gradients = want_expansion_gradients;

 //
 // I/O info
 //
 struct IO_info& IO_info = state.IO_info;
 IO_info.output_ASCII_files = (output_ASCII_files != 0);
 IO_info.output_HDF5_files = (output_HDF5_files != 0);
 IO_info.output_initial_guess = (output_initial_guess != 0);
 IO_info.output_h_every     = output_h_every;
 IO_info.output_Theta_every = output_Theta_every;
 IO_info.output_mean_curvature_every = output_mean_curvature_every;
 IO_info.output_h     = false;	// dummy value
 IO_info.output_Theta = false;	// dummy value
 IO_info.output_mean_curvature = false;	// dummy value

 IO_info.output_BH_diagnostics              = (output_BH_diagnostics != 0);
 IO_info.BH_diagnostics_directory
   = (strlen(BH_diagnostics_directory) == 0)
   ? /* IO:: */ /*out_dir*/ cur_directory
   : BH_diagnostics_directory;
 IO_info.BH_diagnostics_base_file_name      = BH_diagnostics_base_file_name;
 IO_info.BH_diagnostics_file_name_extension = BH_diagnostics_file_name_extension;

 IO_info.output_ghost_zones_for_h  = (output_ghost_zones_for_h != 0);
 IO_info.ASCII_gnuplot_file_name_extension = ASCII_gnuplot_file_name_extension;
 IO_info.HDF5_file_name_extension          = HDF5_file_name_extension;
 IO_info.h_directory
   = (strlen(h_directory) == 0)
   ? /* IO:: */ /*out_dir*/cur_directory
   : h_directory;
 IO_info.h_base_file_name         = h_base_file_name;
 IO_info.Theta_base_file_name     = Theta_base_file_name;
 IO_info.mean_curvature_base_file_name     = mean_curvature_base_file_name;
 IO_info.Delta_h_base_file_name   = Delta_h_base_file_name;
 IO_info.h_min_digits             = h_min_digits;
 IO_info.Jacobian_base_file_name  = Jacobian_base_file_name;
 IO_info.output_OpenDX_control_files  = (output_OpenDX_control_files != 0);
 IO_info.OpenDX_control_file_name_extension = OpenDX_control_file_name_extension;
 IO_info.time_iteration = 0;
 IO_info.time           = 0.0;

 //
 // other misc setup
 //
 state.BH_diagnostics_info.integral_method
   = patch::decode_integration_method(integral_method);


 //
 // mask parameters
 //
 // struct mask_info& mask_info = state.mask_info;
 // mask_info.set_mask_for_any_horizon = false;
 // // ... horizon numbers run from 1 to N_horizons inclusive
 // //     so the array size is N_horizons+1
 // mask_info.set_mask_for_this_horizon = new bool[N_horizons+1];
 // {
 //   for (int hn = 1 ; hn <= N_horizons ; ++hn)
 //   {
 //     mask_info.set_mask_for_this_horizon[hn]
 //       = (set_mask_for_all_horizons != 0)
 //       || (set_mask_for_individual_horizon[hn] != 0);
 //     mask_info.set_mask_for_any_horizon
 //       |= mask_info.set_mask_for_this_horizon[hn];
 //   }
 // }
 // if (mask_info.set_mask_for_any_horizon)
 //   then {
 //     mask_info.radius_multiplier  = mask_radius_multiplier;
 //     mask_info.radius_offset      = mask_radius_offset;
 //     mask_info.buffer_thickness   = mask_buffer_thickness;
 //     mask_info.mask_is_noshrink   = mask_is_noshrink;
 //     mask_info.min_horizon_radius_points_for_mask
 //       = min_horizon_radius_points_for_mask;
 //     mask_info.set_old_style_mask = (set_old_style_mask != 0);
 //     mask_info.set_new_style_mask = (set_new_style_mask != 0);
 //     if (mask_info.set_old_style_mask)
 //       then {
 //         struct mask_info::old_style_mask_info& osmi
 //           = mask_info.old_style_mask_info;
 //         osmi.gridfn_name     = old_style_mask_gridfn_name;
 //         osmi.gridfn_varindex = Cactus_gridfn_varindex(osmi.gridfn_name);
 //         osmi.gridfn_dataptr  = NULL;	// dummy value; fixup later
 //         osmi.inside_value  = old_style_mask_inside_value;
 //         osmi.buffer_value  = old_style_mask_buffer_value;
 //         osmi.outside_value = old_style_mask_outside_value;
 //       }
 //     if (mask_info.set_new_style_mask)
 //       then {
 //         struct mask_info::new_style_mask_info& nsmi
 //           = mask_info.new_style_mask_info;
 //         nsmi.gridfn_name     = new_style_mask_gridfn_name;
 //         nsmi.gridfn_varindex = Cactus_gridfn_varindex(nsmi.gridfn_name);
 //         nsmi.gridfn_dataptr  = NULL;	// dummy value; fixup later
 //         nsmi.bitfield_name   = new_style_mask_bitfield_name;
 //         nsmi.bitfield_bitmask = 0;	// dummy value; fixup later
 //         nsmi.inside_value    = new_style_mask_inside_value;
 //         nsmi.buffer_value    = new_style_mask_buffer_value;
 //         nsmi.outside_value   = new_style_mask_outside_value;
 //         nsmi.inside_bitvalue = 0;	// dummy value; fixup later
 //         nsmi.buffer_bitvalue = 0;	// dummy value; fixup later
 //         nsmi.outside_bitvalue = 0;	// dummy value; fixup later
 //       }
 //   }


 //
 // (genuine) horizon sequence for this processor
 //
 state.my_hs = new horizon_sequence(state.N_horizons);
 horizon_sequence& hs = *state.my_hs;

 //
 // if we're going to actually find horizons
 //    we spread the horizons over multiple processors for maximum efficiency,
 // otherwise (we're just doing testing/debugging computations, so)
 //    we allocate all the horizons to processor #0 for simplicity
 //
 const bool multiproc_flag = (state.method == method__find_horizons);
 state.N_active_procs
   = allocate_horizons_to_processor(state.N_procs, state.my_proc,
				    state.N_horizons, multiproc_flag,
                                    depends_on,
				    hs,
				    verbose_info);

 // ... horizon numbers run from 1 to N_horizons inclusive
 //     so the array size is N_horizons+1
 state.AH_data_array = new AH_data*[N_horizons+1];
 {
   for (int hn = 0 ; hn <= N_horizons ; ++hn)
   {
     state.AH_data_array[hn] = NULL;
   }
 }


 //
 // horizon-specific info for each horizon
 //

 // set up the interpatch interpolator
 if (verbose_info.print_algorithm_highlights)
   then CCTK_VInfo(CCTK_THORNSTRING, "   setting up interpatch interpolator");
 const int ip_interp_handle = 0;
 //const int ip_interp_handle = CCTK_InterpHandle(interpatch_interpolator_name);
 if (ip_interp_handle < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "AHFinderDirect_setup(): couldn't find interpatch interpolator \"%s\"!",
        	   interpatch_interpolator_name);		/*NOTREACHED*/
 //const int ip_interp_param_table_handle = 0;
 const int ip_interp_param_table_handle
   = Util_TableCreateFromString(interpatch_interpolator_pars);
 if (ip_interp_param_table_handle < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "AHFinderDirect_setup(): bad interpatch-interpolator parameter(s) \"%s\"!",
        	   interpatch_interpolator_pars);		/*NOTREACHED*/

 // set up the surface interpolator if it's going to be used
 int surface_interp_handle = -1;
 int surface_interp_param_table_handle = -1;

 if (strlen(surface_interpolator_name) > 0)
   then {
     if (verbose_info.print_algorithm_highlights)
       then CCTK_VInfo(CCTK_THORNSTRING,
                       "   setting up surface interpolator");
     //surface_interp_handle = CCTK_InterpHandle(surface_interpolator_name);
     surface_interp_handle = 0;
     if (surface_interp_handle < 0)
       then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "AHFinderDirect_setup(): couldn't find surface interpolator \"%s\"!",
                       surface_interpolator_name);		/*NOTREACHED*/
     surface_interp_param_table_handle
       = Util_TableCreateFromString(surface_interpolator_pars);
     if (surface_interp_param_table_handle < 0)
       then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "AHFinderDirect_setup(): bad surface-interpolator parameter(s) \"%s\"!",
                       surface_interpolator_pars);		/*NOTREACHED*/
   }
 // setup all horizons on this processor,
 // with full-fledged patch systems for genuine horizons
 // and skeletal patch systems for others
 {
   for (int hn = 1 ; hn <= hs.N_horizons() ; ++hn)
   {
     const bool genuine_flag = hs.is_hn_genuine(hn);

     state.AH_data_array[hn] = new AH_data;
     struct AH_data& AH_data = *state.AH_data_array[hn];

     if (verbose_info.print_algorithm_highlights)
       then CCTK_VInfo(CCTK_THORNSTRING,
                       "   setting up %s data structures for horizon %d",
                       (genuine_flag ? "full-fledged" : "skeletal"),
                       hn);

     // decide what type of patch system this one should be

     const enum patch_system::patch_system_type ps_type =
       patch_system::type_of_name(patch_system_type[hn]); 

     // const enum patch_system::patch_system_type ps_type
     //   = STRING_EQUAL(patch_system_type[hn], "match Cactus grid symmetry")
     //   ? // choose a patch system type based on grid:: parameters
     //   // ... we inherit from grid, and we ask for some of its
     //   //     parameters in our param.ccl file; these appear as
     //   //     ordinary Cactus parameters here, so (eg)
     //   //     grid::domain is just "domain" here
     //   choose_patch_system_type(/* grid:: */ domain,
     //                            /* grid:: */ bitant_plane,
     //                            /* grid:: */ quadrant_direction,
     //                            /* grid:: */ rotation_axis,
     //                            origin_x[hn],
     //                            origin_y[hn],
     //                            origin_z[hn])
     //   : patch_system::type_of_name(patch_system_type[hn]);
     // create the patch system
     AH_data.ps_ptr
       = new patch_system(origin_x[hn], origin_y[hn], origin_z[hn],
                          ps_type,
                          ghost_zone_width, patch_overlap_width,
                          N_zones_per_right_angle[hn],
                          gfns::nominal_min_gfn,
                          (genuine_flag ? gfns::nominal_max_gfn
                           : gfns::skeletal_nominal_max_gfn),
                          gfns::ghosted_min_gfn, gfns::ghosted_max_gfn,
                          ip_interp_handle, ip_interp_param_table_handle,
                          surface_interp_handle,
                          surface_interp_param_table_handle,
                          true, verbose_info.print_algorithm_details);

     patch_system& ps = *AH_data.ps_ptr;
     if (genuine_flag)
       then ps.set_gridfn_to_constant(0.0, gfns::gfn__zero);
     if (genuine_flag)
       then ps.set_gridfn_to_constant(1.0, gfns::gfn__one);

     AH_data.Jac_ptr = genuine_flag
       ? new_Jacobian(Jac_info.Jacobian_store_solve_method,
                      ps,
                      verbose_info.print_algorithm_details)
       : NULL;

     AH_data.compute_info.surface_definition =
       STRING_EQUAL(surface_definition[hn], "expansion")
       ? definition_expansion
       : STRING_EQUAL(surface_definition[hn], "inner expansion")
       ? definition_inner_expansion
       : STRING_EQUAL(surface_definition[hn], "mean curvature")
       ? definition_mean_curvature
       : STRING_EQUAL(surface_definition[hn], "expansion product")
       ? definition_expansion_product
       : (CCTK_WARN (0, "internal error"), definition_error);
     AH_data.compute_info.surface_modification =
       STRING_EQUAL(surface_modification[hn], "none")
       ? modification_none
       : STRING_EQUAL(surface_modification[hn], "radius")
       ? modification_radius
       : STRING_EQUAL(surface_modification[hn], "radius^2")
       ? modification_radius2
#if 0
       : STRING_EQUAL(surface_modification[hn], "mean radius")
       ? modification_mean_radius
       : STRING_EQUAL(surface_modification[hn], "areal radius")
       ? modification_areal_radius
#endif
       : (CCTK_WARN (0, "internal error"), modification_error);
     AH_data.compute_info.surface_selection = 
       STRING_EQUAL(surface_selection[hn], "definition")
       ? selection_definition
       : STRING_EQUAL(surface_selection[hn], "mean coordinate radius")
       ? selection_mean_coordinate_radius
       : STRING_EQUAL(surface_selection[hn], "areal radius")
       ? selection_areal_radius
       : STRING_EQUAL(surface_selection[hn], "expansion times mean coordinate radius")
       ? selection_expansion_mean_coordinate_radius
       : STRING_EQUAL(surface_selection[hn], "expansion times areal radius")
       ? selection_expansion_areal_radius
       : (CCTK_WARN (0, "internal error"), selection_error);
     AH_data.compute_info.desired_value = desired_value[hn];

     AH_data.move_origins               = move_origins;

     AH_data.use_pretracking            = use_pretracking[hn];
     AH_data.pretracking_max_iterations = pretracking_max_iterations[hn];

     AH_data.pretracking_value          = pretracking_value[hn];
     AH_data.pretracking_minimum_value  = pretracking_minimum_value[hn];
     AH_data.pretracking_maximum_value  = pretracking_maximum_value[hn];
     AH_data.pretracking_delta          = pretracking_delta[hn];
     AH_data.pretracking_minimum_delta  = pretracking_minimum_delta[hn];
     AH_data.pretracking_maximum_delta  = pretracking_maximum_delta[hn];

     AH_data.depends_on           = depends_on[hn];
     AH_data.desired_value_factor = desired_value_factor[hn];
     AH_data.desired_value_offset = desired_value_offset[hn];

     AH_data.shiftout_factor = shiftout_factor[hn];
     AH_data.smoothing_factor = smoothing_factor[hn];

     // AH_data.initial_find_flag = genuine_flag;
     // AH_data.really_initial_find_flag = AH_data.initial_find_flag;
     AH_data.initial_find_flag = true;
     AH_data.really_initial_find_flag = AH_data.initial_find_flag;

     if (genuine_flag)
       then {
         if (verbose_info.print_algorithm_details)
           then CCTK_VInfo(CCTK_THORNSTRING,
			   "      setting initial guess parameters etc");
         set_initial_guess_parameters(AH_data, hn, /* irrelevant here; leave at zero */0, 0, 0);
       }

     AH_data.search_flag = false;
     AH_data.found_flag = false;
     AH_data.h_files_written = false;
     AH_data.BH_diagnostics_fileptr = NULL;
   }
 }

}


bool Horizon::broadcast_status(
		      int N_procs, int N_active_procs,
		      int my_proc, bool my_active_flag,
		      int hn, int iteration,
		      enum expansion_status effective_expansion_status,
		      fp mean_horizon_radius, fp infinity_norm,
		      bool found_this_horizon, bool I_need_more_iterations,
		      struct iteration_status_buffers& isb)
{
  assert( my_proc >= 0 );
  assert( my_proc < N_procs );

  //
  // We do the broadcast via a Cactus reduction operation (this is a KLUDGE,
  // but until Cactus gets a generic interprocessor communications mechanism
  // it's the best we can do without mucking with MPI ourself).  To do the
  // gather via a reduction, we set up a 2-D buffer whose entries are all
  // 0s, except that on each processor the [my_proc] row has the values we
  // want to gather.  Then we do a sum-reduction of the buffers across
  // processors.  For the actual reduction we treat the buffers as 1-D
  // arrays on each processor (this slightly simplifies the code).
  //
  // To reduce overheads, we do the entire operation with a single (CCTK_REAL)
  // Cactus reduce, casting values to CCTK_REAL as necessary (and encoding
  // Boolean flags into signs of known-to-be-positive values.
  //
  // Alas MPI (and thus Cactus) requires the input and output reduce buffers
  // to be distinct, so we need two copies of the buffers on each processor.
  //

  // the buffers are actually 2-D arrays; these are the column numbers
  // ... if we wanted to, we could pack hn, iteration, and
  //     effective_expansion_status all into a single (64-bit)
  //     floating-point value, but it's not worth the trouble...
  enum	{
    // CCTK_INT buffer
    buffer_var__hn = 0,	// also encodes found_this_horizon flag
    // in sign: +=true, -=false
    buffer_var__iteration,	// also encodes I_need_more_iterations flag
				// in sign: +=true, -=false
    buffer_var__expansion_status,
    buffer_var__mean_horizon_radius,
    buffer_var__Theta_infinity_norm,
    N_buffer_vars // no comma
  };

  //
  // allocate buffers if this is the first use
  //
  if (isb.hn_buffer == NULL)
    then {
      isb.hn_buffer                  = new int [N_active_procs];
      isb.iteration_buffer           = new int [N_active_procs];
      isb.expansion_status_buffer = new enum expansion_status[N_active_procs];
      isb.mean_horizon_radius_buffer = new fp  [N_active_procs];
      isb.Theta_infinity_norm_buffer = new fp  [N_active_procs];
      isb.found_horizon_buffer       = new bool[N_active_procs];

      isb.send_buffer_ptr    = new jtutil::array2d<CCTK_REAL>
        (0, N_active_procs-1,
         0, N_buffer_vars-1);
      isb.receive_buffer_ptr = new jtutil::array2d<CCTK_REAL>
        (0, N_active_procs-1,
         0, N_buffer_vars-1);
    }
  jtutil::array2d<CCTK_REAL>& send_buffer    = *isb.send_buffer_ptr;
  jtutil::array2d<CCTK_REAL>& receive_buffer = *isb.receive_buffer_ptr;
  //
  // pack this processor's values into the reduction buffer
  //
  jtutil::zero_C_array(send_buffer.N_array(), send_buffer.data_array());
  if (my_active_flag)
    then {
      assert( send_buffer.is_valid_i(my_proc) );
      assert( hn >= 0 );		// encoding scheme assumes this
      assert( iteration > 0 );	// encoding scheme assumes this
      send_buffer(my_proc, buffer_var__hn)
        = found_this_horizon ? +hn : -hn;
      send_buffer(my_proc, buffer_var__iteration)
        = I_need_more_iterations ? +iteration : -iteration;
      send_buffer(my_proc, buffer_var__expansion_status)
        = int(effective_expansion_status);
      send_buffer(my_proc, buffer_var__mean_horizon_radius)
        = mean_horizon_radius;
      send_buffer(my_proc, buffer_var__Theta_infinity_norm)
        = infinity_norm;
    }

  //
  // do the reduction
  //


  hierarchy->getMPI().Allreduce(
    static_cast<void*>(send_buffer.data_array()),
    static_cast<void*>(receive_buffer.data_array()),
    send_buffer.N_array(),
    MPI_DOUBLE,
    MPI_SUM);

  //
  // unpack the reduction buffer back to the high-level result buffers and
  // compute the inclusive-or of the broadcast I_need_more_iterations flags
  //
  bool any_proc_needs_more_iterations = false;
  for (int proc = 0 ; proc < N_active_procs ; ++proc)
  {
    const int hn_temp = static_cast<int>(
      receive_buffer(proc, buffer_var__hn)
    );

    isb.hn_buffer[proc] = jtutil::abs(hn_temp);
    isb.found_horizon_buffer[proc] = (hn_temp > 0);

    const int iteration_temp = static_cast<int>(
      receive_buffer(proc, buffer_var__iteration)
    );
    isb.iteration_buffer[proc] = jtutil::abs(iteration_temp);
    const bool proc_needs_more_iterations = (iteration_temp > 0);
    any_proc_needs_more_iterations |= proc_needs_more_iterations;

    isb.expansion_status_buffer[proc]
      = static_cast<enum expansion_status>(
        static_cast<int>(
          receive_buffer(proc, buffer_var__expansion_status)
        )
      );

    isb.mean_horizon_radius_buffer[proc]
      = receive_buffer(proc, buffer_var__mean_horizon_radius);
    isb.Theta_infinity_norm_buffer[proc]
      = receive_buffer(proc, buffer_var__Theta_infinity_norm);
  }

  return any_proc_needs_more_iterations;
  
}

//******************************************************************************

//
// This function (which must be called on *every* processor) broadcasts
// the BH diagnostics and (ghosted) horizon shape from a specified processor
// to all processors.
//
// The present implementation of this function uses the Cactus reduction
// API.  If AHFinderDirect is ported to some other software environment,
// it's probably best to re-implement this function on top of whatever
// interprocessor-broadcast facility that environment provides.
//
// Arguments:
// GH --> The Cactus grid hierarchy.
// broadcast_flag = true on the broadcasting processor,
//		    false on all other processors
// BH_diagnostics = On the broadcasting processor, this is the BH diagnostics
//		    to broadcast; on all other processors, it's set to the
//		    broadcast BH diagnostics.
// ps = On the broadcasting processor,  gfn__h  is broadcast;
//      on all other processors,  gfn__h  is set to the broadcast values.
// horizon_buffers = Internal buffers for use in the broadcast;
//		     if  N_buffer == 0  then we set N_buffer and allocate
//		     the buffers.
//

void Horizon::broadcast_horizon_data(
			    bool broadcast_flag, bool broadcast_horizon_shape,
                            struct AH_data& AH_data,
			    struct BH_diagnostics& BH_diagnostics,
			    patch_system& ps,
			    struct horizon_buffers& horizon_buffers)
{
//
// Implementation notes:
//
// We do the send via a Cactus sum-reduce where the data are 0 on
// all processors except the sending one.
//
// To reduce the interprocessor-communications cost, we actually only
// broadcast the nominal-grid horizon shape; we then do a synchronize
// operation on the patch system to recover the full ghosted-grid shape.
//

if (horizon_buffers.N_buffer == 0)
   then {
	// allocate the buffers
	horizon_buffers.N_buffer
		= BH_diagnostics::N_buffer
		  + (broadcast_horizon_shape ? ps.N_grid_points() : 0)
                  + 4;
	horizon_buffers.send_buffer
		= new CCTK_REAL[horizon_buffers.N_buffer];
	horizon_buffers.receive_buffer
		= new CCTK_REAL[horizon_buffers.N_buffer];
	}

if (broadcast_flag)
   then {
	// pack the data to be broadcast into the send buffer
	BH_diagnostics.copy_to_buffer(horizon_buffers.send_buffer);
        int posn = BH_diagnostics::N_buffer;
	if (broadcast_horizon_shape)
	   then {
			for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
			{
			const patch& p = ps.ith_patch(pn);
				for (int irho = p.min_irho() ;
				     irho <= p.max_irho() ;
				     ++irho)
				{
				for (int isigma = p.min_isigma() ;
				     isigma <= p.max_isigma() ;
				     ++isigma)
				{
				assert( posn < horizon_buffers.N_buffer );
				horizon_buffers.send_buffer[posn++]
				  = p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
				}
				}
			}
		}
        horizon_buffers.send_buffer[posn++] = AH_data.initial_find_flag;
        horizon_buffers.send_buffer[posn++] = AH_data.really_initial_find_flag;
        horizon_buffers.send_buffer[posn++] = AH_data.search_flag;
        horizon_buffers.send_buffer[posn++] = AH_data.found_flag;
        assert( posn == horizon_buffers.N_buffer );
	}
   else jtutil::zero_C_array(horizon_buffers.N_buffer,
			     horizon_buffers.send_buffer);


 hierarchy->getMPI().Allreduce(
   static_cast<void*>(horizon_buffers.send_buffer),
   static_cast<void*>(horizon_buffers.receive_buffer),
   horizon_buffers.N_buffer,
   MPI_DOUBLE,
   MPI_SUM);

if (!broadcast_flag)
   then {
	// unpack the data from the receive buffer
	BH_diagnostics.copy_from_buffer(horizon_buffers.receive_buffer);
        ps.origin_x(BH_diagnostics.origin_x);
        ps.origin_y(BH_diagnostics.origin_y);
        ps.origin_z(BH_diagnostics.origin_z);
	int posn = BH_diagnostics::N_buffer;
	if (broadcast_horizon_shape)
	   then {
			for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
			{
			patch& p = ps.ith_patch(pn);
				for (int irho = p.min_irho() ;
				     irho <= p.max_irho() ;
				     ++irho)
				{
				for (int isigma = p.min_isigma() ;
				     isigma <= p.max_isigma() ;
				     ++isigma)
				{
				assert( posn < horizon_buffers.N_buffer );
				p.ghosted_gridfn(gfns::gfn__h, irho,isigma)
				   = horizon_buffers.receive_buffer[posn++];
				}
				}
			}

		// recover the full ghosted-grid horizon shape
		// (we only broadcast the nominal-grid shape)
		ps.synchronize();
		}
        AH_data.initial_find_flag        = horizon_buffers.receive_buffer[posn++];
        AH_data.really_initial_find_flag = horizon_buffers.receive_buffer[posn++];
        AH_data.search_flag              = horizon_buffers.receive_buffer[posn++];
        AH_data.found_flag               = horizon_buffers.receive_buffer[posn++];
	assert( posn == horizon_buffers.N_buffer );
	}
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function is called on processor #0 to print the status of the
// Newton iteration on each active processor.
//
// Arguments:
// N_active_procs = The number of active processors.
// isb = The high-level buffers here give the information to be printed.
//

void Horizon::print_status(int N_active_procs,
		  const struct iteration_status_buffers& isb)
{
  for (int proc = 0 ; proc < N_active_procs ; ++proc)
  {
    // don't print anything for processors doing dummy evaluations
    if (isb.hn_buffer[proc] == 0)
      then continue;

    if (isb.expansion_status_buffer[proc] == expansion_success)
      then CCTK_VInfo(CCTK_THORNSTRING,
                      "   proc %d/horizon %d:it %d r_grid=%#.3g ||Theta||=%.1e",
                      proc, isb.hn_buffer[proc],
                      isb.iteration_buffer[proc],
                      double(isb.mean_horizon_radius_buffer[proc]),
                      double(isb.Theta_infinity_norm_buffer[proc]));
    else CCTK_VInfo(CCTK_THORNSTRING,
                    "   proc %d/horizon %d: %s",
                    proc, isb.hn_buffer[proc],
                    expansion_status_string(
                      isb.expansion_status_buffer[proc]
                    ));
  }
}

//******************************************************************************

//
// This function takes the Newton step, scaling it down if it's too large.
//
// Arguments:
// ps = The patch system containing the gridfns h and Delta_h.
// mean_horizon_radius = ||h||_mean
// max_allowable_Delta_h_over_h = The maximum allowable
//				     ||Delta_h||_infinity / ||h||_mean
//				  Any step over this is internally clamped
//				  (scaled down) to this size.
//

void Horizon::Newton_step(patch_system& ps,
		 fp mean_horizon_radius, fp max_allowable_Delta_h_over_h,
		 const struct verbose_info& verbose_info)
{
//
// compute scale factor (1 for small steps, <1 for large steps)
//

const fp max_allowable_Delta_h
	= max_allowable_Delta_h_over_h * mean_horizon_radius;

jtutil::norm<fp> Delta_h_norms;
ps.gridfn_norms(gfns::gfn__Delta_h, Delta_h_norms);
const fp max_Delta_h = Delta_h_norms.infinity_norm();

const fp scale = (max_Delta_h <= max_allowable_Delta_h)
		 ? 1.0
		 : max_allowable_Delta_h / max_Delta_h;

if (verbose_info.print_algorithm_details)
   then {

	if (scale == 1.0)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "h += Delta_h (rms-norm=%.1e, infinity-norm=%.1e)",
			   Delta_h_norms.rms_norm(),
			   Delta_h_norms.infinity_norm());
	   else CCTK_VInfo(CCTK_THORNSTRING,
			   "h += %g * Delta_h (infinity-norm clamped to %.2g)",
			   scale,
			   scale * Delta_h_norms.infinity_norm());
	}


//
// take the Newton step (scaled if necessary)
//
	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		p.ghosted_gridfn(gfns::gfn__h, irho,isigma)
			-= scale * p.gridfn(gfns::gfn__Delta_h, irho,isigma);
		}
		}
	}

if (ps.N_additional_points())
	{
	const int np = ps.N_grid_points();
	const int gnp = ps.ghosted_N_grid_points();
	ps.ghosted_gridfn_data(gfns::gfn__h)[gnp]
		-= scale * ps.gridfn_data(gfns::gfn__Delta_h)[np];
	}
}

//******************************************************************************
//
// This function reads a coordinate origin from a grid scalar, 
// and sets the patch system's origin to that new value
//
//
void Horizon::track_origin(patch_system& ps, 
                  struct AH_data* const AH_data_ptr, 
                  const int hn, const bool print_algorithm_details)
{
  // TODOMARKS
  // figuring out new way of doing thsi
if (AH_data_ptr->depends_on == 0) {
        // move the origin as specified in the grid scalars
        // fp const * const ox =
        //   static_cast<CCTK_REAL const *>
        //   (CCTK_VarDataPtr (cctkGH, 0, track_origin_source_x[hn]));
        // assert (ox);
        // fp const * const oy =
        //   static_cast<CCTK_REAL const *>
        //   (CCTK_VarDataPtr (cctkGH, 0, track_origin_source_y[hn]));
        // assert (oy);
        // fp const * const oz =
        //   static_cast<CCTK_REAL const *>
        //   (CCTK_VarDataPtr (cctkGH, 0, track_origin_source_z[hn]));
        // assert (oz);
        // if (print_algorithm_details) {
        //   std::cout << "AHF tracked position ox " << *ox << " oy " << *oy << " oz " << *oz << std::endl;
  tbox::pout<<"Entering the function track_origin, not sure whether do it right!\n";
        ps.origin_x (0);
        ps.origin_y (0);
        ps.origin_z (0);
        };
}

//******************************************************************************

//
// This function tries to finds each horizon assigned to this processor,
// by solving Theta(h) = 0 via Newton's method.  For each horizon found,
// it computes the BH diagnostics, optionally prints them (via CCTK_VInfo()),
// and optionally writes them to a BH diagnostics output file.  It also
// optionally writes the horizon shape itself to a data file.
//
// This function must be called synchronously across all processors,
// and it operates synchronously, returning only when every processor
// is done with all its horizons.  See ./README.parallel for a discussion
// of the parallel/multiprocessor issues and algorithm.
//
 void Horizon::Newton(
	    int N_procs, int N_active_procs, int my_proc,
	    horizon_sequence& hs, struct AH_data* const AH_data_array[],
	    const struct cactus_grid_info& cgi,
	    const struct geometry_info& gi,
	    const struct Jacobian_info& Jacobian_info,
	    const struct solver_info& solver_info,
	    const struct IO_info& IO_info,
	    const struct BH_diagnostics_info& BH_diagnostics_info,
	    bool broadcast_horizon_shape,
	    const struct error_info& error_info,
	    const struct verbose_info& verbose_info,
	    struct iteration_status_buffers& isb)
{

const bool my_active_flag = hs.has_genuine_horizons();

// print out which horizons we're finding on this processor
if (hs.has_genuine_horizons())
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "proc %d: searching for horizon%s %s/%d",
		   my_proc,
		   (hs.my_N_horizons() > 1 ? "s" : ""),
		   hs.sequence_string(","), int(N_horizons));
   else CCTK_VInfo(CCTK_THORNSTRING,
		   "proc %d: dummy horizon(s) only",
		   my_proc);

    //
    // each pass through this loop finds a single horizon,
    // or does corresponding dummy-horizon calls
    //
    // note there is no explicit exit test, rather we exit from the middle
    // of the loop (only) when all processors are done with all their genuine
    // horizons
    //
    for (int hn = hs.init_hn() ; ; hn = hs.next_hn())
    {
    if (verbose_info.print_algorithm_details)
       then CCTK_VInfo(CCTK_THORNSTRING,
		       "Newton_solve(): processor %d working on horizon %d",
		       my_proc, hn);

    // only try to find horizons every  find_every  time steps
    const bool horizon_is_genuine =
      hs.is_genuine() && AH_data_array[hn]->search_flag;
    // this is only a pessimistic approximation
    const bool there_is_another_genuine_horizon = hs.is_next_genuine();
    if (verbose_info.print_algorithm_details)
       then {
	    CCTK_VInfo(CCTK_THORNSTRING,
		       "                horizon_is_genuine=%d",
		       int(horizon_is_genuine));
	    CCTK_VInfo(CCTK_THORNSTRING,
		       "                there_is_another_genuine_horizon=%d",
		       int(there_is_another_genuine_horizon));
	    }

    struct AH_data* AH_data_ptr
	                  = horizon_is_genuine ? AH_data_array[hn]   : NULL;

    patch_system* const  ps_ptr = horizon_is_genuine
				  ? AH_data_ptr->ps_ptr : NULL;
    Jacobian*     const Jac_ptr = horizon_is_genuine
				  ? AH_data_ptr->Jac_ptr: NULL;

    if (horizon_is_genuine) {
      // deal with dependent horizons
      if (AH_data_ptr->depends_on > 0) {
        assert (AH_data_ptr->depends_on != hn);
        assert (AH_data_ptr->depends_on < hn);
        // check that the other horizon is on the same processor!
        AH_data *AH_other_ptr = AH_data_array[AH_data_ptr->depends_on];
        assert (AH_other_ptr);
        AH_data_ptr->compute_info.desired_value
          = AH_other_ptr->compute_info.desired_value
          * AH_data_ptr->desired_value_factor
          + AH_data_ptr->desired_value_offset;
        const int gnp = ps_ptr->ghosted_N_grid_points();
        assert (AH_other_ptr->ps_ptr->ghosted_N_grid_points() == gnp);
        memcpy (ps_ptr->ghosted_gridfn_data(gfns::gfn__h),
                AH_other_ptr->ps_ptr->ghosted_gridfn_data(gfns::gfn__h),
                gnp * sizeof(fp));
        ps_ptr->origin_x (AH_other_ptr->ps_ptr->origin_x());
        ps_ptr->origin_y (AH_other_ptr->ps_ptr->origin_y());
        ps_ptr->origin_z (AH_other_ptr->ps_ptr->origin_z());
      }
    }

    if (horizon_is_genuine) {
      if (AH_data_ptr->move_origins
          && AH_data_ptr->depends_on == 0
          && AH_data_ptr->found_flag)
      {
        // move the origin to the centre
        patch_system& ps = *ps_ptr;
        fp cx=0, cy=0, cz=0;    // change of origin
        fp sx=0, sy=0, sz=0;    // change of shape
        if (! predict_origin_movement) {
          size_t N=0;
          for (int pn=0; pn<ps.N_patches(); ++pn) {
            patch& p = ps.ith_patch(pn);
            for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
              for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
                const fp rho = p.rho_of_irho (irho);
                const fp sigma = p.sigma_of_isigma (isigma);
                const fp r = p.ghosted_gridfn (gfns::gfn__h, irho, isigma);
                fp x, y, z;
                p.xyz_of_r_rho_sigma (r, rho, sigma, x, y, z);
                cx += x;
                cy += y;
                cz += z;
                ++N;
              }
            }
          }
          assert (N > 0);
          cx /= N;
          cy /= N;
          cz /= N;
          sx=cx; sy=cy; sz=cz;
        } else {TBOX_ERROR("Cannot predict centroid movement yet!\n");                // if predict_origin_movement
          // if (ah_centroid_valid[hn-1]) {
          //   if (verbose_info.print_algorithm_details) {
          //     std::cout << "AHF centroid x " << (ah_centroid_x[hn-1]) << " y " << (ah_centroid_y[hn-1]) << " z " << (ah_centroid_z[hn-1]) << " t " << (ah_centroid_t[hn-1]) << std::endl;
          //   }
          //   if (ah_centroid_valid_p[hn-1]) {
          //     // have two previous origins: linear extrapolation
          //     if (verbose_info.print_algorithm_details) {
          //       std::cout << "AHF centroid_p x " << (ah_centroid_x_p[hn-1]) << " y " << (ah_centroid_y_p[hn-1]) << " z " << (ah_centroid_z_p[hn-1]) << " t " << (ah_centroid_t_p[hn-1]) << std::endl;
          //     }
          //     fp const dt   = ah_centroid_t[hn-1] - ah_centroid_t_p[hn-1];
          //     fp const dt_n = cctk_time - ah_centroid_t  [hn-1];
          //     fp const dt_p = cctk_time - ah_centroid_t_p[hn-1];
          //     fp const timescale =
          //       fabs (dt) + fabs (dt_p) + fabs (dt_n) +
          //       0.001 * fabs (cctk_delta_time);
          //     if (10 * fabs (dt  ) < timescale ||
          //         10 * fabs (dt_p) < timescale ||
          //         10 * fabs (dt_n) < timescale)
          //     {
          //       // the times are too similar
          //       if (verbose_info.print_algorithm_details) {
          //         std::cout << "AHF toosim" << std::endl;
          //       }
          //       cx = ah_centroid_x[hn-1];
          //       cy = ah_centroid_y[hn-1];
          //       cz = ah_centroid_z[hn-1];
          //     } else {
          //       fp const fact_n = + dt_p / dt;
          //       fp const fact_p = - dt_n / dt;
          //       if (fabs (fact_n) > 5 || fabs (fact_p) > 5) {
          //         // don't trust a large extrapolation
          //         if (verbose_info.print_algorithm_details) {
          //           std::cout << "AHF notrust" << std::endl;
          //         }
          //         cx = ah_centroid_x[hn-1];
          //         cy = ah_centroid_y[hn-1];
          //         cz = ah_centroid_z[hn-1];
          //       } else {
          //         if (verbose_info.print_algorithm_details) {
          //           std::cout << "AHF extrap fn " << fact_n << " fp " << fact_p << std::endl;
          //         }
          //         cx = fact_n * ah_centroid_x[hn-1] + fact_p * ah_centroid_x_p[hn-1];
          //         cy = fact_n * ah_centroid_y[hn-1] + fact_p * ah_centroid_y_p[hn-1];
          //         cz = fact_n * ah_centroid_z[hn-1] + fact_p * ah_centroid_z_p[hn-1];
          //         if (verbose_info.print_algorithm_details) {
          //           std::cout << "AHF xp " << (ah_centroid_x_p[hn-1]) << " x " << (ah_centroid_x[hn-1]) << " cx " << cx << std::endl;
          //         }
          //       }
          //     }
          //   } else {
          //     // have one previous origin: constant extrapolation
          //     if (verbose_info.print_algorithm_details) {
          //       std::cout << "AHF const" << std::endl;
          //     }
          //     cx = ah_centroid_x[hn-1];
          //     cy = ah_centroid_y[hn-1];
          //     cz = ah_centroid_z[hn-1];
          //   }
          //   if (verbose_info.print_algorithm_details) {
          //     std::cout << "AHF predicted position cx " << cx << " cy " << cy << " cz " << cz << std::endl;
          //   }
          //   cx -= ps.origin_x();
          //   cy -= ps.origin_y();
          //   cz -= ps.origin_z();
          // }
          // sx = ah_centroid_x[hn-1] - ps.origin_x();
          // sy = ah_centroid_y[hn-1] - ps.origin_y();
          // sz = ah_centroid_z[hn-1] - ps.origin_z();
        } // if predict_origin_movement
        switch (ps.type()) {
        case patch_system::patch_system__full_sphere:
          break;                // do nothing
        case patch_system::patch_system__plus_z_hemisphere:
          cz = 0; sz = 0; break;
        case patch_system::patch_system__plus_xy_quadrant_mirrored:
          cx = cy = 0; sx = sy = 0; break;
        case patch_system::patch_system__plus_xy_quadrant_rotating:
          cx = cy = 0; sx = sy = 0; break;
        case patch_system::patch_system__plus_xz_quadrant_mirrored:
          cx = cz = 0; sx = sz = 0; break;
        case patch_system::patch_system__plus_xz_quadrant_rotating:
          cx = cz = 0; sx = sz = 0; break;
        case patch_system::patch_system__plus_xyz_octant_mirrored:
          cx = cy = cz = 0; sx = sy = sz = 0; break;
        case patch_system::patch_system__plus_xyz_octant_rotating:
          cx = cy = cz = 0; sx = sy = sz = 0; break;
        default: assert(0);
        }
        ps.origin_x (ps.origin_x() + cx);
        ps.origin_y (ps.origin_y() + cy);
        ps.origin_z (ps.origin_z() + cz);
        if (reshape_while_moving) {
          for (int pn=0; pn<ps.N_patches(); ++pn) {
            patch& p = ps.ith_patch(pn);
            for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
              for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
                const fp rho = p.rho_of_irho (irho);
                const fp sigma = p.sigma_of_isigma (isigma);
                fp & r = p.ghosted_gridfn (gfns::gfn__h, irho, isigma);
                fp x, y, z;
                p.xyz_of_r_rho_sigma (r, rho, sigma, x, y, z);
                fp const proj = (sx*x + sy*y + sz*z) / r;
                r -= proj;
              }
            }
          }
          ps.synchronize();
        }
      }
      if (track_origin_from_grid_scalar[hn]) {
         track_origin(*ps_ptr, AH_data_ptr, hn, verbose_info.print_algorithm_details);
      }
      
      // modify the initial guess
      jtutil::norm<fp> norms;
      ps_ptr->ghosted_gridfn_norms (gfns::gfn__h, norms);
      // smoothing:
      ps_ptr->scale_ghosted_gridfn
        (1.0 - AH_data_ptr->smoothing_factor, gfns::gfn__h);
      ps_ptr->add_to_ghosted_gridfn
        (AH_data_ptr->smoothing_factor * norms.mean(), gfns::gfn__h);
      // enlarging:
      ps_ptr->scale_ghosted_gridfn
        (AH_data_ptr->shiftout_factor, gfns::gfn__h);
    }

    bool I_am_pretracking = horizon_is_genuine && AH_data_ptr->use_pretracking;
    bool I_was_pretracking = false;
    bool pretracking_have_upper_bound = false;
    bool pretracking_have_lower_bound = false;
    bool pretracking_was_successful = false;
    fp const old_pretracking_value = I_am_pretracking ? AH_data_ptr->pretracking_value : 0.0;
    fp pretracking_upper_bound;
    fp pretracking_lower_bound;
    bool pretracking_have_horizon_info;
    fp pretracking_mean_expansion;
    for (int pretracking_iter = 0;
         I_am_pretracking
           ? pretracking_iter < AH_data_ptr->pretracking_max_iterations
           : true;
         ++pretracking_iter)
      {
        bool found_this_horizon;
        if (I_am_pretracking || I_was_pretracking) {
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Pretracking: iteration %d", pretracking_iter);
          if (pretracking_have_lower_bound) {
            if (pretracking_have_upper_bound) {
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Pretracking: looking for value %g in [%g,%g]",
                         double(AH_data_ptr->pretracking_value),
                           double(pretracking_lower_bound),
                         double(pretracking_upper_bound));
            } else {
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Pretracking: looking for value %g in [%g,*]",
                         double(AH_data_ptr->pretracking_value),
                         double(pretracking_lower_bound));
            }
          } else {
            if (pretracking_have_upper_bound) {
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Pretracking: looking for value %g in [*,%g]",
                         double(AH_data_ptr->pretracking_value),
                           double(pretracking_upper_bound));
            } else {
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Pretracking: looking for value %g in [*,*]",
                         double(AH_data_ptr->pretracking_value));
            }
          }
          AH_data_ptr->compute_info.desired_value = AH_data_ptr->pretracking_value;
          ps_ptr->ghosted_gridfn_copy (gfns::gfn__h, gfns::gfn__save_h);
          pretracking_have_horizon_info = false;
        }

        if (! I_am_pretracking) {
          if (horizon_is_genuine) {
            if (! AH_data_ptr->initial_guess_info.reset_horizon_after_not_finding) {
              // save the surface for possible backtracking later
              ps_ptr->ghosted_gridfn_copy (gfns::gfn__h, gfns::gfn__save_h);
            }
          }
        }

	struct what_to_compute dummy_compute_info;
	struct what_to_compute & compute_info =
		horizon_is_genuine
		? AH_data_ptr->compute_info
		: dummy_compute_info;

    if (horizon_is_genuine) {
      if (ps_ptr->N_additional_points()) {
        const int gnp = ps_ptr->ghosted_N_grid_points();
        ps_ptr->ghosted_gridfn_data(gfns::gfn__h)[gnp] = 0;
      }
    }

    const int max_iterations
       = horizon_is_genuine
	 ? (AH_data_ptr->initial_find_flag
	    ? solver_info.max_Newton_iterations__initial
	    : solver_info.max_Newton_iterations__subsequent)
	 : INT_MAX;

    int num_Theta_growth_iterations = 0;
    fp previous_Theta_norm = 1.0e30;
    int num_Theta_nonshrink_iterations = 0;
    fp best_Theta_norm = 1.0e30;

	//
	// each pass through this loop does a single Newton iteration
	// on the current horizon (which might be either genuine or dummy)
	//
        bool do_return = false;
	for (int iteration = 1 ; ; ++iteration)
	{
	if (verbose_info.print_algorithm_debug)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "beginning iteration %d (horizon_is_genuine=%d)",
			   iteration, int(horizon_is_genuine));

	//
	// evaluate the expansion Theta(h) and its norms
	// *** this is a synchronous operation across all processors ***
	//
	if (horizon_is_genuine
	    && solver_info.debugging_output_at_each_Newton_iteration)
	   then output_gridfn(*ps_ptr, gfns::gfn__h,
                              "h",
			      IO_info, IO_info.h_base_file_name,
                              IO_info.h_min_digits,
			      hn, verbose_info.print_algorithm_highlights,
			      iteration);

        // calculate the norms also for a surface a bit further out,
        // so that we know how they vary in space.
        // perform this calculation first, so that the real values
        // do not have to be saved.
        const fp epsilon = Jacobian_info.perturbation_amplitude;
	jtutil::norm<fp> shifted_Theta_norms;
        jtutil::norm<fp> shifted_expansion_Theta_norms;
        jtutil::norm<fp> shifted_inner_expansion_Theta_norms;
        jtutil::norm<fp> shifted_product_expansion_Theta_norms;
        jtutil::norm<fp> shifted_mean_curvature_Theta_norms;
        fp shifted_area;
	jtutil::norm<fp> shifted2_Theta_norms;
        jtutil::norm<fp> shifted2_expansion_Theta_norms;
        jtutil::norm<fp> shifted2_inner_expansion_Theta_norms;
        jtutil::norm<fp> shifted2_product_expansion_Theta_norms;
        jtutil::norm<fp> shifted2_mean_curvature_Theta_norms;
        fp shifted2_area;
        
        if (solver_info.want_expansion_gradients) {
          if (horizon_is_genuine) {
            ps_ptr->add_to_ghosted_gridfn(epsilon, gfns::gfn__h);
            // ps_ptr->scale_ghosted_gridfn(1.0+epsilon, gfns::gfn__h);
          }

          const enum expansion_status raw_shifted_expansion_status
            = expansion(ps_ptr, compute_info,
                        cgi, gi,
                        error_info, (iteration == 1),
                        false,	// no, we don't want Jacobian coeffs
                        false,
                        &shifted_Theta_norms,
                        &shifted_expansion_Theta_norms,
                        &shifted_inner_expansion_Theta_norms,
                        &shifted_product_expansion_Theta_norms,
                        &shifted_mean_curvature_Theta_norms);
          if (horizon_is_genuine) {
            shifted_area = ps_ptr->integrate_gridfn
              (gfns::gfn__one, true, true, true,
               gfns::gfn__h,
               gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                                   gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                       gfns::gfn__g_dd_33,
               BH_diagnostics_info.integral_method);
            ps_ptr->add_to_ghosted_gridfn(-epsilon, gfns::gfn__h);
            // ps_ptr->scale_ghosted_gridfn(1.0/(1.0+epsilon), gfns::gfn__h);
            
            ps_ptr->add_to_ghosted_gridfn(-epsilon, gfns::gfn__h);
            // ps_ptr->scale_ghosted_gridfn(1.0/(1.0+epsilon), gfns::gfn__h);
          }
          const enum expansion_status raw_shifted2_expansion_status
            = expansion(ps_ptr, compute_info,
                        cgi, gi,
                        error_info, (iteration == 1),
                        false,	// no, we don't want Jacobian coeffs
                        false,
                        &shifted2_Theta_norms,
                        &shifted2_expansion_Theta_norms,
                        &shifted2_inner_expansion_Theta_norms,
                        &shifted2_product_expansion_Theta_norms,
                        &shifted2_mean_curvature_Theta_norms);
          if (horizon_is_genuine) {
            shifted2_area = ps_ptr->integrate_gridfn
              (gfns::gfn__one, true, true, true,
               gfns::gfn__h,
               gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                                   gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                       gfns::gfn__g_dd_33,
               BH_diagnostics_info.integral_method);
            ps_ptr->add_to_ghosted_gridfn(epsilon, gfns::gfn__h);
            // ps_ptr->scale_ghosted_gridfn(1.0+epsilon, gfns::gfn__h);
          }
        } // if want_expansion_gradients

        // now calculate the real values.
	jtutil::norm<fp> Theta_norms;
        jtutil::norm<fp> expansion_Theta_norms;
        jtutil::norm<fp> inner_expansion_Theta_norms;
        jtutil::norm<fp> product_expansion_Theta_norms;
        jtutil::norm<fp> mean_curvature_Theta_norms;
	const enum expansion_status raw_expansion_status
               = expansion(ps_ptr, compute_info,
			   cgi, gi,
			   error_info, (iteration == 1),
			   true,	// yes, we want Jacobian coeffs
			   verbose_info.print_algorithm_details,
                           &Theta_norms,
                           &expansion_Theta_norms,
                           &inner_expansion_Theta_norms,
                           &product_expansion_Theta_norms,
                           &mean_curvature_Theta_norms);
        fp area;
        
        if (horizon_is_genuine)
	        {
	        area = ps_ptr->integrate_gridfn
			(gfns::gfn__one, true, true, true,
			 gfns::gfn__h,
		         gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
					     gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
							         gfns::gfn__g_dd_33,
	                 BH_diagnostics_info.integral_method);
	        }
	const bool Theta_is_ok = (raw_expansion_status == expansion_success);
        const bool norms_are_ok = horizon_is_genuine && Theta_is_ok;

        if (norms_are_ok) {
          const fp this_norm = Theta_norms.infinity_norm();
          if (this_norm > previous_Theta_norm) {
            ++ num_Theta_growth_iterations;
          } else {
            num_Theta_growth_iterations = 0;
          }
          previous_Theta_norm = this_norm;
          
          if (this_norm >= best_Theta_norm) {
            ++ num_Theta_nonshrink_iterations;
          } else {
            num_Theta_nonshrink_iterations = 0;
            best_Theta_norm = this_norm;
          }
        }

        if (I_am_pretracking && norms_are_ok) {
          pretracking_mean_expansion = expansion_Theta_norms.mean();
          pretracking_have_horizon_info = true;
        }

	if (verbose_info.print_algorithm_debug)
	   then {
		CCTK_VInfo(CCTK_THORNSTRING,
			   "   Newton_solve(): Theta_is_ok=%d",
			   Theta_is_ok);
		if (norms_are_ok)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "   Theta rms-norm %.1e, infinity-norm %.1e",
				   double(Theta_norms.rms_norm()),
				   double(Theta_norms.infinity_norm()));
		}

	if (horizon_is_genuine && Theta_is_ok
	    && solver_info.debugging_output_at_each_Newton_iteration)
	   then {
		output_gridfn(*ps_ptr, gfns::gfn__Theta,
                              "Theta",
			      IO_info, IO_info.Theta_base_file_name,
                              IO_info.h_min_digits,
			      hn, verbose_info.print_algorithm_highlights,
			      iteration);
		output_gridfn(*ps_ptr, gfns::gfn__mean_curvature,
                              "mean_curvature",
			      IO_info, IO_info.mean_curvature_base_file_name,
                              IO_info.h_min_digits,
			      hn, verbose_info.print_algorithm_highlights,
			      iteration);
		}


	//
	// have we found this horizon?
	// if so, compute and output BH diagnostics
	//
                
        found_this_horizon
           = (norms_are_ok
              && (I_was_pretracking
                  ? pretracking_was_successful
                  : Theta_norms.infinity_norm() <= solver_info.Theta_norm_for_convergence));

	if (horizon_is_genuine)
	   then AH_data_ptr->found_flag = found_this_horizon;

	if (found_this_horizon && ! I_am_pretracking)
	   then {
		struct BH_diagnostics& BH_diagnostics
			= AH_data_ptr->BH_diagnostics;
                const fp mean_expansion       = expansion_Theta_norms.mean();
                const fp mean_inner_expansion = inner_expansion_Theta_norms.mean();
                const fp mean_product_expansion = product_expansion_Theta_norms.mean();
                const fp mean_mean_curvature  = mean_curvature_Theta_norms.mean();
                // const fp area_gradient                 = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_area                               - area                              ) / epsilon;
                // const fp mean_expansion_gradient       = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_expansion_Theta_norms.mean()       - expansion_Theta_norms.mean()      ) / epsilon;
                // const fp mean_inner_expansion_gradient = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_inner_expansion_Theta_norms.mean() - inner_expansion_Theta_norms.mean()) / epsilon;
                // const fp mean_product_expansion_gradient = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_product_expansion_Theta_norms.mean() - product_expansion_Theta_norms.mean()) / epsilon;
                // const fp mean_mean_curvature_gradient  = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_mean_curvature_Theta_norms.mean()  - mean_curvature_Theta_norms.mean() ) / epsilon;
                const fp area_gradient                 = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_area                               - shifted2_area                              ) / (2*epsilon);
                const fp mean_expansion_gradient       = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_expansion_Theta_norms.mean()       - shifted2_expansion_Theta_norms.mean()      ) / (2*epsilon);
                const fp mean_inner_expansion_gradient = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_inner_expansion_Theta_norms.mean() - shifted2_inner_expansion_Theta_norms.mean()) / (2*epsilon);
                const fp mean_product_expansion_gradient = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_product_expansion_Theta_norms.mean() - shifted2_product_expansion_Theta_norms.mean()) / (2*epsilon);
                const fp mean_mean_curvature_gradient  = ! solver_info.want_expansion_gradients ? 0.0 : (shifted_mean_curvature_Theta_norms.mean()  - shifted2_mean_curvature_Theta_norms.mean() ) / (2*epsilon);

		BH_diagnostics.compute(*ps_ptr,
                                       area,
                                       mean_expansion,
                                       mean_inner_expansion,
                                       mean_product_expansion,
                                       mean_mean_curvature,
                                       area_gradient,
                                       mean_expansion_gradient,
                                       mean_inner_expansion_gradient,
                                       mean_product_expansion_gradient,
                                       mean_mean_curvature_gradient,
                                       BH_diagnostics_info);

		if (IO_info.output_BH_diagnostics)
		   then {
			if (AH_data_ptr->BH_diagnostics_fileptr == NULL)
			   then AH_data_ptr->BH_diagnostics_fileptr
				  = BH_diagnostics.setup_output_file
					(IO_info, N_horizons, hn);

			BH_diagnostics.output(AH_data_ptr
					      ->BH_diagnostics_fileptr,
					      IO_info);
			}
		}


	//
	// see if the expansion is too big
	// (if so, we'll give up on this horizon)
	//
	const bool expansion_is_too_large
		= norms_are_ok
                  && (   Theta_norms.infinity_norm() > solver_info.max_allowable_Theta
                      || (   solver_info.max_allowable_Theta_growth_iterations > 0
                          && num_Theta_growth_iterations > solver_info.max_allowable_Theta_growth_iterations)
                      || (   solver_info.max_allowable_Theta_nonshrink_iterations > 0
                          && num_Theta_nonshrink_iterations > solver_info.max_allowable_Theta_nonshrink_iterations)
                         );


	//
	// compute the mean horizon radius, and if it's too large,
	// then pretend expansion() returned a "surface too large" error status
	//
	jtutil::norm<fp> h_norms;
	if (horizon_is_genuine) {
          ps_ptr->ghosted_gridfn_norms(gfns::gfn__h, h_norms);
        }
	const fp mean_horizon_radius
		= horizon_is_genuine ? h_norms.mean() : 0.0;
	const bool horizon_is_too_large
		= (mean_horizon_radius > solver_info
					 .max_allowable_horizon_radius[hn]);

	const enum expansion_status effective_expansion_status
		= horizon_is_too_large ? expansion_failure__surface_too_large
				       : raw_expansion_status;


	//
	// see if we need more iterations (either on this or another horizon)
	//

	// does *this* horizon need more iterations?
	// i.e. has this horizon's Newton iteration not yet converged?
        const bool this_horizon_needs_more_iterations
	   = horizon_is_genuine && Theta_is_ok
	     && !found_this_horizon
	     && !expansion_is_too_large
	     && !horizon_is_too_large
	     && (iteration < max_iterations);

	// do I (this processor) need to do more iterations
	// on this or a following horizon?
	const bool I_need_more_iterations
	   = this_horizon_needs_more_iterations
	     || there_is_another_genuine_horizon
             || I_am_pretracking;

	if (verbose_info.print_algorithm_details)
	   then {
		CCTK_VInfo(CCTK_THORNSTRING,
			   "   flags: found_this_horizon=%d",
			   int(found_this_horizon));
		CCTK_VInfo(CCTK_THORNSTRING,
			   "          this_horizon_needs_more_iterations=%d",
			   int(this_horizon_needs_more_iterations));
		CCTK_VInfo(CCTK_THORNSTRING,
			   "          I_need_more_iterations=%d",
			   int(I_need_more_iterations));
		}

	//
	// broadcast iteration status from each active processor
	// to all processors, and inclusive-or the "we need more iterations"
	// flags to see if *any* (active) processor needs more iterations
	//
        
	const bool any_proc_needs_more_iterations
	  = broadcast_status(
			     N_procs, N_active_procs,
			     my_proc, my_active_flag,
			     hn, iteration, effective_expansion_status,
			     mean_horizon_radius,
			     (norms_are_ok ? Theta_norms.infinity_norm() : 0.0),
			     found_this_horizon, I_need_more_iterations,
			     isb);
	// set found-this-horizon flags
	// for all active processors' non-dummy horizons
		  {
		for (int found_proc = 0 ;
		     found_proc < N_active_procs ;
		     ++found_proc)
		{
		const int found_hn = isb.hn_buffer[found_proc];
		if (found_hn > 0)
		   then AH_data_array[found_hn]->found_flag
				= isb.found_horizon_buffer[found_proc];
		}
		  }
	if (verbose_info.print_algorithm_details)
	   then {
		CCTK_VInfo(CCTK_THORNSTRING,
			   "          ==> any_proc_needs_more_iterations=%d",
			   int(any_proc_needs_more_iterations));
		}


	//
	// print the iteration status info
	//
	if ((my_proc == 0) && verbose_info.print_algorithm_highlights)
	   then print_status(N_active_procs, isb);


	//
	// for each processor which found a horizon,
	// broadcast its horizon info to all processors
	// and print the BH diagnostics on processor 0
	//
		  {
		for (int found_proc = 0 ;
		     found_proc < N_active_procs ;
		     ++found_proc)
		{

		if (! isb.found_horizon_buffer[found_proc])
		   then continue;
		const int found_hn = isb.hn_buffer[found_proc];
		struct AH_data& found_AH_data = *AH_data_array[found_hn];

		if (verbose_info.print_algorithm_details)
		   then CCTK_VInfo(CCTK_THORNSTRING,
			   "   broadcasting proc %d/horizon %d diagnostics%s",
				   found_proc, found_hn,
				   (broadcast_horizon_shape ? "+shape" : ""));

		broadcast_horizon_data(
				       my_proc == found_proc,
				       broadcast_horizon_shape,
                                       found_AH_data,
				       found_AH_data.BH_diagnostics,
				       *found_AH_data.ps_ptr,
				       found_AH_data.horizon_buffers);

		if ((my_proc == 0) && verbose_info.print_physics_details)
		   then found_AH_data.BH_diagnostics
				     .print(N_horizons, found_hn);
		}
		  }


	//
	// if we found our horizon, maybe output the horizon shape?
	//
	if (found_this_horizon && ! I_am_pretracking)
	   then {
		// printf("will output h/Th/mc: %d/%d/%d\n", IO_info.output_h, IO_info.output_Theta, IO_info.output_mean_curvature); //xxxxxxxxxxxx
		if (IO_info.output_h)
		   then {
			// if this is the first time we've output h for this
			// horizon, maybe output an OpenDX control file?
			if (!AH_data_ptr->h_files_written)
			   then setup_h_files(*ps_ptr, IO_info, hn);
			output_gridfn(*ps_ptr, gfns::gfn__h,
                                      "h",
				      IO_info, IO_info.h_base_file_name,
                                      IO_info.h_min_digits,
				      hn, verbose_info
					  .print_algorithm_highlights);
			}
		if (IO_info.output_Theta)
		   then output_gridfn(*ps_ptr, gfns::gfn__Theta,
                                      "Theta",
				      IO_info, IO_info.Theta_base_file_name,
                                      IO_info.h_min_digits,
				      hn, verbose_info
					  .print_algorithm_highlights);
		if (IO_info.output_mean_curvature)
		   then output_gridfn(*ps_ptr, gfns::gfn__mean_curvature,
                                      "mean_curvature",
				      IO_info, IO_info.mean_curvature_base_file_name,
                                      IO_info.h_min_digits,
				      hn, verbose_info
					  .print_algorithm_highlights);
		}


	//
	// are all processors done with all their genuine horizons?
	// or if this is a single-processor run, are we done with this horizon?
	//
	if (! any_proc_needs_more_iterations)
	   then {
		if (verbose_info.print_algorithm_details)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "==> all processors are finished!");
                do_return = true;
		break;					// *** LOOP EXIT ***
		}
	if ((N_procs == 1) && !this_horizon_needs_more_iterations)
	   then {
		if (verbose_info.print_algorithm_debug)
		   then CCTK_VInfo(CCTK_THORNSTRING,
			   "==> [single processor] Skipping to next horizon");
		break;					// *** LOOP EXIT ***
		}


	//
	// compute the Jacobian matrix
	// *** this is a synchronous operation across all processors ***
	//
	if (verbose_info.print_algorithm_debug)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "   computing Jacobian: genuine/dummy flag %d",
			   int(this_horizon_needs_more_iterations));
	const enum expansion_status
	  Jacobian_status
	     = expansion_Jacobian
		 (this_horizon_needs_more_iterations ? ps_ptr : NULL,
		  this_horizon_needs_more_iterations ? Jac_ptr : NULL,
		  compute_info,
		  cgi, gi, Jacobian_info,
		  error_info, (iteration == 1),
		  verbose_info.print_algorithm_details);
	const bool Jacobian_is_ok = (Jacobian_status == expansion_success);


	//
	// skip to the next horizon unless
	// this is a genuine Jacobian computation, and it went ok
	//
	if (! (this_horizon_needs_more_iterations && Jacobian_is_ok))
	   then {
		if (verbose_info.print_algorithm_debug)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "   skipping to next horizon");
		break;					// *** LOOP EXIT ***
		}


	//
	// compute the Newton step
	//
	if (verbose_info.print_algorithm_details)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "solving linear system for Delta_h (%d unknowns)",
			   Jac_ptr->N_rows());
	const fp rcond
	   = Jac_ptr->solve_linear_system(gfns::gfn__Theta, gfns::gfn__Delta_h,
					  solver_info.linear_solver_pars,
					  verbose_info.print_algorithm_details);
	if	((rcond >= 0.0) && (rcond < 100.0*FP_EPSILON))
	   then {
		CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "Newton_solve: Jacobian matrix is numerically singular!");
		// give up on this horizon
		break;				// *** LOOP CONTROL ***
		}
	if (verbose_info.print_algorithm_details)
	   then {
		if (rcond > 0.0)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "done solving linear system (rcond=%.1e)",
				   double(rcond));
		   else CCTK_VInfo(CCTK_THORNSTRING,
				   "done solving linear system");
		}

	if (solver_info.debugging_output_at_each_Newton_iteration)
	   then output_gridfn(*ps_ptr, gfns::gfn__Delta_h,
                              "Delta_h", 
			      IO_info, IO_info.Delta_h_base_file_name,
                              IO_info.h_min_digits,
			      hn, verbose_info.print_algorithm_details,
			      iteration);


	//
	// take the Newton step (scaled if need be)
	//
	Newton_step(*ps_ptr,
		    mean_horizon_radius, solver_info
					 .max_allowable_Delta_h_over_h,
		    verbose_info);

	// end of this Newton iteration
	}

        if (! I_am_pretracking) {
          if (horizon_is_genuine) {
            if (! AH_data_ptr->initial_guess_info.reset_horizon_after_not_finding) {
              if (! found_this_horizon) {
                // the surface failed; backtrack and continue
                AH_data_ptr->ps_ptr->ghosted_gridfn_copy
                  (gfns::gfn__save_h, gfns::gfn__h);
              }
            }
          }
        }

        // exit
        if (do_return) return;				// *** NORMAL RETURN ***

        // break out of the pretracking loop if we are not pretracking
        if (! I_am_pretracking) break;
        
        if (! found_this_horizon) {
          // the surface failed; backtrack and continue
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Pretracking: solving failed; backtracking");
          AH_data_ptr->ps_ptr->ghosted_gridfn_copy (gfns::gfn__save_h, gfns::gfn__h);
          if (pretracking_have_lower_bound) {
            assert (AH_data_ptr->pretracking_value >= pretracking_lower_bound - 1.0e-10 * fabs(pretracking_lower_bound));
          }
          pretracking_have_lower_bound = true;
          pretracking_lower_bound = AH_data_ptr->pretracking_value;
          if (pretracking_have_upper_bound) {
            AH_data_ptr->pretracking_delta = 0.5 * (pretracking_upper_bound - pretracking_lower_bound);
            AH_data_ptr->pretracking_value = pretracking_lower_bound + 0.5 * (pretracking_upper_bound - pretracking_lower_bound);
          } else {
            AH_data_ptr->pretracking_delta = fabs(AH_data_ptr->pretracking_delta);
            AH_data_ptr->pretracking_delta *= 2.0;
            if (AH_data_ptr->pretracking_delta > AH_data_ptr->pretracking_maximum_delta) {
              AH_data_ptr->pretracking_delta = AH_data_ptr->pretracking_maximum_delta;
            }
            AH_data_ptr->pretracking_value += AH_data_ptr->pretracking_delta;
            if (AH_data_ptr->pretracking_value > AH_data_ptr->pretracking_maximum_value) {
              AH_data_ptr->pretracking_value = AH_data_ptr->pretracking_maximum_value;
            }
          }
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Pretracking: new value %g, delta %g",
                     double(AH_data_ptr->pretracking_value),
                     double(AH_data_ptr->pretracking_delta));
          if (pretracking_lower_bound >= AH_data_ptr->pretracking_maximum_value * 0.9999999999) {
            // give up
            CCTK_VInfo(CCTK_THORNSTRING,
                       "Pretracking: upper bound reached, giving up.");
            I_am_pretracking = false;
            I_was_pretracking = true;
            pretracking_was_successful = false;
            // restore old pretracking goal
            AH_data_ptr->pretracking_value = old_pretracking_value;
          } else if (AH_data_ptr->pretracking_delta < AH_data_ptr->pretracking_minimum_delta) {
            // give up
            CCTK_VInfo(CCTK_THORNSTRING,
                       "Pretracking: step size too small, giving up.");
            I_am_pretracking = false;
            I_was_pretracking = true;
            pretracking_was_successful = true;
          }
        } else {
          // the surface was okay
          // get mean expansion
          assert (pretracking_have_horizon_info);
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Pretracking: solving succeeded; expansion is now %g",
                     double(pretracking_mean_expansion));
//           if (fabs(AH_data_ptr->pretracking_value) > 2.0 * AH_data_ptr->pretracking_minimum_delta) {
          if (AH_data_ptr->pretracking_value > AH_data_ptr->pretracking_minimum_value + 1.0e-10 * AH_data_ptr->pretracking_minimum_delta) {
            // it is not yet a horizon
            if (pretracking_have_upper_bound) {
              assert (AH_data_ptr->pretracking_value <= pretracking_upper_bound + 1.0e-10 * fabs(pretracking_upper_bound));
            }
            pretracking_have_upper_bound = true;
            pretracking_upper_bound = AH_data_ptr->pretracking_value;
            if (pretracking_have_lower_bound) {
#if 1
              // TODO
              // move lower bound further down
              pretracking_lower_bound -= pretracking_upper_bound - pretracking_lower_bound;
              if (pretracking_lower_bound < AH_data_ptr->pretracking_minimum_value) pretracking_lower_bound = AH_data_ptr->pretracking_minimum_value;
#endif
              AH_data_ptr->pretracking_delta = 0.5 * (pretracking_lower_bound - pretracking_upper_bound);
              AH_data_ptr->pretracking_value = pretracking_lower_bound + 0.5 * (pretracking_upper_bound - pretracking_lower_bound);
            } else {
              AH_data_ptr->pretracking_delta = - fabs(AH_data_ptr->pretracking_delta);
              AH_data_ptr->pretracking_delta *= 2.0;
              if (- AH_data_ptr->pretracking_delta > AH_data_ptr->pretracking_maximum_delta) {
                AH_data_ptr->pretracking_delta = - AH_data_ptr->pretracking_maximum_delta;
              }
              AH_data_ptr->pretracking_value += AH_data_ptr->pretracking_delta;
              if (AH_data_ptr->pretracking_value < AH_data_ptr->pretracking_minimum_value) {
                AH_data_ptr->pretracking_value = AH_data_ptr->pretracking_minimum_value;
              }
            }
            CCTK_VInfo(CCTK_THORNSTRING,
                       "Pretracking: new value radius %g, delta %g",
                       double(AH_data_ptr->pretracking_value),
                       double(AH_data_ptr->pretracking_delta));
            if (pretracking_upper_bound <= AH_data_ptr->pretracking_minimum_value * 1.00000000001) {
              // give up
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Pretracking: lower bound reached, giving up.");
              I_am_pretracking = false;
              I_was_pretracking = true;
              pretracking_was_successful = false;
              // restore old pretracking goal
              AH_data_ptr->pretracking_value = old_pretracking_value;
            } else if (- AH_data_ptr->pretracking_delta < AH_data_ptr->pretracking_minimum_delta) {
              // give up
              CCTK_VInfo(CCTK_THORNSTRING,
                         "Pretracking: step size too small, giving up.");
              I_am_pretracking = false;
              I_was_pretracking = true;
              pretracking_was_successful = true;
            }
          } else {
            // a true horizon was found; we are done
            CCTK_VInfo(CCTK_THORNSTRING,
                       "Pretracking: done.");
            I_am_pretracking = false;
            I_was_pretracking = true;
            pretracking_was_successful = true;
          }
        } // if surface found

    // end of pretracking loop
    }
    I_am_pretracking = false;

    // end of this horizon
    }

// we should never get to here
assert( false );
}

//******************************************************************************

//
// This function is called (via the magic of function aliasing) by
// other thorns to find out our local coordinate origin for a given AH.
//
// Results:
// This function returns 0 for ok, or -1 if the horizon number is invalid.
//
CCTK_INT Horizon::AHFinderDirect_local_coordinate_origin
    (CCTK_INT horizon_number,
     CCTK_REAL* p_origin_x, CCTK_REAL* p_origin_y, CCTK_REAL* p_origin_z)
{
  const struct verbose_info& verbose_info = state.verbose_info;

  if (!  ((horizon_number >= 1) && (horizon_number <= state.N_horizons))  )
    then {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "AHFinderDirect_local_coordinate_origin():\n"
                 "        horizon_number=%d must be in the range [1,N_horizons=%d]!\n"
                 ,
                 int(horizon_number), state.N_horizons);
      return -1;					// *** ERROR RETURN ***
    }

  assert(state.AH_data_array[horizon_number] != NULL);
  const struct AH_data& AH_data = *state.AH_data_array[horizon_number];

  assert(AH_data.ps_ptr != NULL);
  const patch_system& ps = *AH_data.ps_ptr;

  assert(p_origin_x != NULL);
  assert(p_origin_y != NULL);
  assert(p_origin_z != NULL);
  *p_origin_x = ps.origin_x();
  *p_origin_y = ps.origin_y();
  *p_origin_z = ps.origin_z();

  return 0;						// *** NORMAL RETURN ***
}

//******************************************************************************

//
// This function is called (via the Cactus flesh function-aliasing mechanism)
// by other thorns to query whether or not the specified horizon was found
// the last time we searched for it.
//
// Results:
// This function returns
//  1 if the horizon was found
//  0 if the horizon was not found
//  -1 if the horizon number is invalid.
//
CCTK_INT Horizon::AHFinderDirect_horizon_was_found(CCTK_INT horizon_number)
{
const struct verbose_info& verbose_info = state.verbose_info;

if (!  ((horizon_number >= 1) && (horizon_number <= state.N_horizons))  )
   then {
	CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"AHFinderDirect_horizon_was_found():\n"
"        horizon_number=%d must be in the range [1,N_horizons=%d]!\n"
		   ,
		   int(horizon_number), state.N_horizons);
	return -1;					// *** ERROR RETURN ***
	}

assert(state.AH_data_array[horizon_number] != NULL);
const struct AH_data& AH_data = *state.AH_data_array[horizon_number];

return AH_data.found_flag ? 1 : 0;
}

//******************************************************************************

//
// This function is called (via the Cactus flesh function-aliasing mechanism)
// by other thorns to query whether or not the specified horizon was found
// the last time we searched for it, and if so, to determine the horizon
// centroid.
//
// Results:
// This function returns:
//  1 if the horizon was found; in this case  *centroid_[xyz]_ptr
//    are set to the centroid position 
//  0 if the horizon was not found; in this case  *centroid_[xyz]_ptr
//    set to zeros
//  negative for an error
//
CCTK_INT Horizon::AHFinderDirect_horizon_centroid
    (CCTK_INT horizon_number,
     CCTK_REAL* p_centroid_x, CCTK_REAL* p_centroid_y, CCTK_REAL* p_centroid_z)
{
const struct verbose_info& verbose_info = state.verbose_info;

if (!  ((horizon_number >= 1) && (horizon_number <= state.N_horizons))  )
   then {
	CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
"AHFinderDirect_horizon_centroid():\n"
"        horizon_number=%d must be in the range [1,N_horizons=%d]!\n"
		   ,
		   int(horizon_number), state.N_horizons);
	return -1;					// *** ERROR RETURN ***
	}

assert(state.AH_data_array[horizon_number] != NULL);
const struct AH_data& AH_data = *state.AH_data_array[horizon_number];
const struct BH_diagnostics& BH_diagnostics = AH_data.BH_diagnostics;

assert(p_centroid_x != NULL);
assert(p_centroid_y != NULL);
assert(p_centroid_z != NULL);
if (AH_data.found_flag)
   then {
	*p_centroid_x = BH_diagnostics.centroid_x;
	*p_centroid_y = BH_diagnostics.centroid_y;
	*p_centroid_z = BH_diagnostics.centroid_z;
	return 1;
	}
   else {
	*p_centroid_x = 0.0;
	*p_centroid_y = 0.0;
	*p_centroid_z = 0.0;
	return 0;
	}
}

//******************************************************************************

//
// This function is called (via the Cactus flesh function-aliasing mechanism)
// by other thorns to find out a given AH's radius in the direction from
// its local coordinate origin to a given (x,y,z) coordinate or coordinates.
//
// Arguments:
// horizon_number = must be in the range 1 to N_horizons
// N_points = should be >= 0
// x[], y[], z[] = these give the (x,y,z) coordinates
// radius[] = this is set to the horizon radius values (Euclidean distance
//	      from the local coordinate origin), or to all -1.0 if we didn't
//	      find this horizon the last time we looked for it
//
// Results:
// This function returns 0 for ok, or -1 if the horizon number is invalid.
//
CCTK_INT Horizon::AHFinderDirect_radius_in_direction
    (CCTK_INT horizon_number,
     CCTK_INT N_points,
     const CCTK_REAL* const x, const CCTK_REAL* const y, const CCTK_REAL* const z,
     CCTK_REAL* const radius)
{
  const struct verbose_info& verbose_info = state.verbose_info;
  
  if (! ((horizon_number >= 1) && (horizon_number <= state.N_horizons)) ) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "AHFinderDirect_distance_outside_thorn():\n"
               "        horizon_number=%d must be in the range [1,N_horizons=%d]!\n"
               ,
               int(horizon_number), state.N_horizons);
    return -1;					// *** ERROR RETURN ***
  }
  
  assert(state.AH_data_array[horizon_number] != NULL);
  const struct AH_data& AH_data = *state.AH_data_array[horizon_number];
  
  if (AH_data.found_flag) {
    
    assert(AH_data.ps_ptr != NULL);
    const patch_system& ps = *AH_data.ps_ptr;
    
    std::vector<fp> local_xs(N_points);
    std::vector<fp> local_ys(N_points);
    std::vector<fp> local_zs(N_points);
    
    for (int point = 0 ; point < N_points ; ++point) {
      
      local_xs.at(point) = x[point] - ps.origin_x();
      local_ys.at(point) = y[point] - ps.origin_y();
      local_zs.at(point) = z[point] - ps.origin_z();
      
    }
    
    ps.radii_in_local_xyz_directions (gfns::gfn__h,
                                      N_points,
                                      & local_xs.front(),
                                      & local_ys.front(),
                                      & local_zs.front(),
                                      radius);
    
  } else {
    // if not found
    
    for (int point = 0 ; point < N_points ; ++point) {
      radius[point] = -1.0;
    }
    
  } // if not found

return 0;						// *** NORMAL RETURN ***
}


//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function sets up the global xyz positions of the grid points
// in the gridfns global_[xyz].  These will be used by interplate_geometry().
void Horizon::setup_xyz_posns(patch_system& ps, const bool print_msg_flag)
{
  if (print_msg_flag)
    then CCTK_VInfo(CCTK_THORNSTRING,
                    "      xyz positions and derivative coefficients");

#ifdef GEOMETRY_INTERP_DEBUG2
  printf ("AH exp\n");
#endif
  for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
  {
    patch& p = ps.ith_patch(pn);

    for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
    {
      for (int isigma = p.min_isigma() ;
           isigma <= p.max_isigma() ;
           ++isigma)
      {
        const fp r = p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
        const fp rho = p.rho_of_irho(irho);
        const fp sigma = p.sigma_of_isigma(isigma);
        fp local_x, local_y, local_z;
        p.xyz_of_r_rho_sigma(r,rho,sigma, local_x,local_y,local_z);
#ifdef GEOMETRY_INTERP_DEBUG2
        printf ("   pn=%d irho=%d isigma=%d   x=%g y=%g z=%g\n",
                pn, irho, isigma, local_x, local_y, local_z);
#endif

        const fp global_x = ps.origin_x() + local_x;
        const fp global_y = ps.origin_y() + local_y;
        const fp global_z = ps.origin_z() + local_z;

        p.gridfn(gfns::gfn__global_x, irho,isigma) = global_x;
        p.gridfn(gfns::gfn__global_y, irho,isigma) = global_y;
        p.gridfn(gfns::gfn__global_z, irho,isigma) = global_z;

        const fp global_xx = global_x * global_x;
        const fp global_xy = global_x * global_y;
        const fp global_xz = global_x * global_z;
        const fp global_yy = global_y * global_y;
        const fp global_yz = global_y * global_z;
        const fp global_zz = global_z * global_z;

        p.gridfn(gfns::gfn__global_xx, irho,isigma) = global_xx;
        p.gridfn(gfns::gfn__global_xy, irho,isigma) = global_xy;
        p.gridfn(gfns::gfn__global_xz, irho,isigma) = global_xz;
        p.gridfn(gfns::gfn__global_yy, irho,isigma) = global_yy;
        p.gridfn(gfns::gfn__global_yz, irho,isigma) = global_yz;
        p.gridfn(gfns::gfn__global_zz, irho,isigma) = global_zz;
      }
    }
  }
}

//******************************************************************************

//
// If ps_ptr != NULL, this function interpolates the Cactus gridfns
//	gxx...gzz
//	kxx...kzz
//	psi			# optional
// to determine the nominal-grid angular gridfns
//	g_dd_ij
//	partial_d_g_dd_kij
//	K_dd_ij
//	psi			# optional
//	partial_d_psi_k		# optional
// at the nominal-grid trial horizon surface positions given by the
// global_(x,y,z) angular gridfns in the patch system *ps_ptr.  The psi
// interpolation is only done if the cgi.use_Cactus_conformal_metric flag
// is set.  Note that this function ignores the physical-vs-conformal
// semantics of the gridfns; it just interpolates and takes derivatives
// of the stored gridfn values.
//
// If ps_ptr == NULL, this function does (only) the parameter-table
// setup and a a dummy interpolator call, as described in the comments
// to  expansion()  above.
//
// The interpolation is done via  CCTK_InterpGridArrays() .  This has the
// option to return both an overall interpolation status, and a "local"
// status which gives the results of interpolating only the points requested
// on *this* processor; if the local status is available we use it, otherwise
// we fall back to the overall status.
//
// Inputs (angular gridfns, all on the nominal grid):
//	global_[xyz]			# xyz positions of grid points
//
// Inputs (Cactus 3-D gridfns):
//      ahmask                          # excision mask
//	gxx,gxy,gxz,gyy,gyz,gzz		# 3-metric $g_{ij}$
//					# (may be either physical or conformal)
//	kxx,kxy,kxz,kyy,kyz,kzz		# extrinsic curvature $K_{ij}$
//	psi				# optional conformal factor $\psi$
//
// Outputs (angular gridfns, all on the nominal grid):
//      mask                                    # excision mask
//      partial_d_mask_[123]                    # derivatives of the mask
//	g_dd_{11,12,13,22,23,33}		# $\stored{g}_{ij}$
//	K_dd_{11,12,13,22,23,33}		# $K_{ij}$
//	partial_d_g_dd_[123]{11,12,13,22,23,33}	# $\partial_k \stored{g}_{ij}$
//	psi					# (optional) $\psi$
//	partial_d_psi_[123]			# (optional) $\partial_k \psi$
//
// This function may also modify the interpolator parameter table.
//
// Results:
// This function returns a status code indicating whether the computation
// succeeded or failed, and if the latter, what caused the failure.  Possible
// failure codes are
// * expansion_failure__surface_outside_grid
// * expansion_failure__surface_in_excised_region	// not implemented yet
//


enum expansion_status
  Horizon::interpolate_geometry(patch_system* ps_ptr,
			  const struct cactus_grid_info& cgi,
			  const struct geometry_info& gi,
			  const struct error_info& error_info, bool initial_flag,
			  bool print_msg_flag)
{
const bool active_flag = (ps_ptr != NULL);
const bool psi_flag = cgi.use_Cactus_conformal_metric;

//
// Implementation Notes:
//
// To handle the optional interpolation of psi, we set up all the data
// type and pointer arrays to include psi, but with the psi entries at
// the end.  We then choose the array sizes passed to the interpolator
// to either include or exclude the psi entries as appropriate.
//
// We remember whether or not psi was interpolated on the previous call,
// and only modify the interpolator parameter table if this changes (or
// if this is our first call).
//

if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "      interpolating %s from Cactus grid",
		   (psi_flag ? "{g_ij, K_ij, psi}" : "{g_ij, K_ij}"));

int status;

#define CAST_PTR_OR_NULL(type_,ptr_)	\
	(ps_ptr == NULL) ? NULL : static_cast<type_>(ptr_)


//
// ***** interpolation points *****
//
const int N_interp_points = (ps_ptr == NULL) ? 0 : ps_ptr->N_grid_points();
const int interp_coords_type_code = CCTK_VARIABLE_REAL;
const void* const interp_coords[N_GRID_DIMS]
  = {
    CAST_PTR_OR_NULL(const void*, ps_ptr->gridfn_data(gfns::gfn__global_x)),
    CAST_PTR_OR_NULL(const void*, ps_ptr->gridfn_data(gfns::gfn__global_y)),
    CAST_PTR_OR_NULL(const void*, ps_ptr->gridfn_data(gfns::gfn__global_z)),
    };


//
// ***** input arrays *****
//

const CCTK_INT input_array_variable_indices[]
	= {
          cgi.mask_varindex,
	  cgi.g_dd_11_varindex, cgi.g_dd_12_varindex, cgi.g_dd_13_varindex,
				cgi.g_dd_22_varindex, cgi.g_dd_23_varindex,
						      cgi.g_dd_33_varindex,
	  cgi.K_dd_11_varindex, cgi.K_dd_12_varindex, cgi.K_dd_13_varindex,
				cgi.K_dd_22_varindex, cgi.K_dd_23_varindex,
						      cgi.K_dd_33_varindex,
	  cgi.psi_varindex,
	  };
const int N_input_arrays_for_psi = 1;
const int N_input_arrays_dim =   sizeof(input_array_variable_indices)
			       / sizeof(input_array_variable_indices[0]);
const int N_input_arrays_use
	= psi_flag ? N_input_arrays_dim
		   : N_input_arrays_dim - N_input_arrays_for_psi;


//
// ***** output arrays *****
//

const CCTK_INT output_array_type_codes[]
	= {
 // mask             $\partial_x$ mask   $\partial_y$ mask   $\partial_z$ mask
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 // $g_{ij}$         $\partial_x g_{ij}$ $\partial_y g_{ij}$ $\partial_z g_{ij}$
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
 // $K_{ij}$
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
		     CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
					 CCTK_VARIABLE_REAL,
 // $\psi$           $\partial_x \psi$   $\partial_y \psi$   $\partial_z \psi$
 CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
	  };

const CCTK_INT operand_indices[]
	= {
	  0, 0, 0, 0,		// mask, partial_[xyz] mask
	  1, 1, 1, 1,		// g_dd_11, partial_[xyz] g_dd_11
	  2, 2, 2, 2,		// g_dd_12, partial_[xyz] g_dd_12
	  3, 3, 3, 3,		// g_dd_13, partial_[xyz] g_dd_13
	  4, 4, 4, 4,		// g_dd_22, partial_[xyz] g_dd_22
	  5, 5, 5, 5,		// g_dd_23, partial_[xyz] g_dd_23
	  6, 6, 6, 6,		// g_dd_33, partial_[xyz] g_dd_33
	  7, 8, 9, 10, 11, 12,	// K_dd_{11,12,13,22,23,33}
	  13, 13, 13, 13,	// psi, partial_[xyz] psi
	  };
#define DERIV(x)	x
const CCTK_INT operation_codes[]
  = {
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // mask, partial_[xyz] mask
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // g_dd_11, partial_[xyz] g_dd_11
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // g_dd_12, partial_[xyz] g_dd_12
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // g_dd_13, partial_[xyz] g_dd_13
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // g_dd_22, partial_[xyz] g_dd_22
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // g_dd_23, partial_[xyz] g_dd_23
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // g_dd_33, partial_[xyz] g_dd_33
    DERIV(0), DERIV(0), DERIV(0), DERIV(0), DERIV(0), DERIV(0),
					    // K_dd_{11,12,13,22,23,33}
    DERIV(0), DERIV(1), DERIV(2), DERIV(3), // psi, partial_[xyz] psi
    };

void* const output_arrays[]
  = {
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__mask)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_mask_1)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_mask_2)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_mask_3)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__g_dd_11)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_111)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_211)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_311)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__g_dd_12)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_112)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_212)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_312)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__g_dd_13)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_113)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_213)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_313)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__g_dd_22)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_122)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_222)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_322)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__g_dd_23)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_123)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_223)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_323)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__g_dd_33)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_133)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_233)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_333)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__K_dd_11)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__K_dd_12)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__K_dd_13)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__K_dd_22)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__K_dd_23)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__K_dd_33)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__psi)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_psi_1)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_psi_2)),
    CAST_PTR_OR_NULL(void*, ps_ptr->gridfn_data(gfns::gfn__partial_d_psi_3)),
    };

const int N_output_arrays_for_psi = 4;
const int N_output_arrays_dim
	= sizeof(output_arrays) / sizeof(output_arrays[0]);
const int N_output_arrays_use
	= psi_flag ? N_output_arrays_dim
		   : N_output_arrays_dim - N_output_arrays_for_psi;

//
// ***** parameter table *****
//

// this flag is true if and only if the parameter table already has the
//	suppress_warnings
//	operand_indices
//	operand_codes
// entries for  psi_flag == par_table_psi_flag .
static bool par_table_setup = false;

// if  par_table_setup,
//    this flag is the value of  psi_flag  for the parameter table entries
// otherwise this flag is ignored
static bool par_table_psi_flag = false;

// TODOMARKS
// do need this setting anyway
if (par_table_setup && (psi_flag == par_table_psi_flag))
   then {
	// parameter table is already set to just what we need
	// ==> no-op here
	}
   else {
	// store derivative info in interpolator parameter table
	if (print_msg_flag)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "         setting up interpolator derivative info");


	par_table_setup = true;
	par_table_psi_flag = psi_flag;
	}


//
// ***** the actual interpolation *****
//
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "         calling geometry interpolator (%s%d points)",
		   (active_flag ? "" : "dummy: "), N_interp_points);

#ifdef GEOMETRY_INTERP_DEBUG2
	  {
	printf("AHFinderDirect:: proc %d: CCTK_InterpGridArrays() coordinates are:\n",
	       int(CCTK_MyProc(cgi.GH)));
	for (int pt = 0 ; pt < N_interp_points ; ++ pt)
	{
	printf("   pt=%d   [x,y,z]=[%g,%g,%g]\n",
               pt,
               double(((const CCTK_REAL*)interp_coords[0])[pt]),
               double(((const CCTK_REAL*)interp_coords[1])[pt]),
               double(((const CCTK_REAL*)interp_coords[2])[pt]));
	}
	  }
#endif	/* GEOMETRY_INTERP_DEBUG2 */

#ifdef GEOMETRY_INTERP_DEBUG
printf("AHFinderDirect:: proc %d: initializing interpolator outputs to 999.999\n",
       int(CCTK_MyProc(cgi.GH)));
	  {
	for (int pt = 0 ; pt < N_interp_points ; ++pt)
	{
		for (int out = 0 ; out < N_output_arrays_use ; ++out)
		{
		CCTK_REAL* const out_ptr
			= static_cast<CCTK_REAL*>(output_arrays[out]);
		out_ptr[pt] = 999.999;
		}
	}
	  }
#endif


          hierarchy->getMPI().Barrier();


 // all coordinate information is stored in interp_coords
 status = CCTK_InterpGridArrays(N_GRID_DIMS,
        		       gi.operator_handle, gi.param_table_handle,
        		       cgi.coord_system_handle,
        		       N_interp_points,
        			  interp_coords_type_code,
        			  interp_coords,
        		       N_input_arrays_use,
        			  input_array_variable_indices,
        		       N_output_arrays_use,
        			  output_array_type_codes,
        			  output_arrays);


#ifdef GEOMETRY_INTERP_DEBUG2
	  {
	for (int pt = 0 ; pt < N_interp_points ; pt = 2*pt + (pt == 0))
	{
	printf("AHFinderDirect:: proc %d: CCTK_InterpGridArrays() results for pt=%d at [x,y,z]=[%g,%g,%g]:\n",
	       int(CCTK_MyProc(cgi.GH)), pt,
               double(((const CCTK_REAL*)interp_coords[0])[pt]),
               double(((const CCTK_REAL*)interp_coords[1])[pt]),
               double(((const CCTK_REAL*)interp_coords[2])[pt]));
		for (int out = 0 ; out < N_output_arrays_use ; ++out)
		{
		const CCTK_REAL* const out_ptr
			= static_cast<const CCTK_REAL*>(output_arrays[out]);
		printf("   out=%d   result=%g\n", out, double(out_ptr[pt]));
		}
	}
	  }
#endif	/* GEOMETRY_INTERP_DEBUG2 */



if (status < 0)
   then error_exit(ERROR_EXIT,
"***** interpolate_geometry(): error return %d from interpolator!\n",
		   status);					/*NOTREACHED*/
// TODOMARKS
// excision mask has not been supported yet
//
// ***** check the interpolated excision mask *****
//
// {
// bool did_use_excised_gridpoint = false;
// if (active_flag)
//    then {
//         for (int pn = 0 ; pn < ps_ptr->N_patches() ; ++pn)
//             {
//             patch& p = ps_ptr->ith_patch(pn);
        
//             for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
//             for (int isigma = p.min_isigma() ; isigma <= p.max_isigma() ; ++isigma)
//         	{
//                 const fp m = p.gridfn(gfns::gfn__mask, irho,isigma);
//                 const fp m1 = p.gridfn(gfns::gfn__partial_d_mask_1, irho,isigma);
//                 const fp m2 = p.gridfn(gfns::gfn__partial_d_mask_2, irho,isigma);
//                 const fp m3 = p.gridfn(gfns::gfn__partial_d_mask_3, irho,isigma);
//                 if (fabs(m) > 1.0e-12
//                     || fabs(m1) > 1.0e-12 || fabs(m2) > 1.0e-12 || fabs(m3) > 1.0e-12)
//                    then did_use_excised_gridpoint = true;
//                 }
//             }
//         }
// if (gi.mask_is_noshrink && did_use_excised_gridpoint)
//    then {
// 	if (print_msg_flag)
// 	   then {
// 		// see if we can get further info
// 		const int warn_level
// 		   = initial_flag
// 		     ? error_info.warn_level__point_outside__initial
// 		     : error_info.warn_level__point_outside__subsequent;


// 		CCTK_VWarn(warn_level, __LINE__, __FILE__, CCTK_THORNSTRING,
// "interpolate_geometry():\n"
// "        one or more points on the trial horizon surface point\n"
// "        is/are in an excised region (or too close to the excision boundary)\n");
// 		}

// 	return expansion_failure__surface_in_excised_region;	// *** ERROR RETURN ***
// 	}
// }

return expansion_success;				// *** NORMAL RETURN ***
}


//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// If ps_ptr != NULL, this function computes the LHS function Theta(h),
// and optionally also its Jacobian coefficients (from which the Jacobian
// matrix may be computed later).
//
// If ps_ptr == NULL, this function does a dummy computation, described
// below.
//
// Inputs (angular gridfns, on ghosted grid):
// ... defined on ghosted grid
// ... only values on nominal grid are actually used as input
//	h				# shape of trial surface
//
// Inputs (Cactus 3-D gridfns):
//	gxx,gxy,gxz,gyy,gyz,gzz		# 3-metric $g_{ij}$
//	kxx,kxy,kxz,kyy,kyz,kzz		# extrinsic curvature $K_{ij}$
//
// Outputs (temporaries computed at each grid point)
//	## computed by hand-written code
//	global_[xyz]			# xyz positions of grid points
//	X_ud_*, X_udd_*			# xyz derivative coefficients
//	## computed by Maple-generated code
//	g_uu_{11,12,13,22,23,33}	# $g^{ij}$
//	K				# $K$
//	K_dd_{11,12,13,22,23,33}	# $K^{ij}$
//	partial_d_ln_sqrt_g_d		# $\partial_i \ln \sqrt{g}$
//	partial_d_g_uu_{1,2,3}{11,12,13,22,23,33}	# $\partial_k g^{ij}$
//
// Outputs (angular gridfns, all on the nominal grid):
//	## interpolated from 3-D Cactus grid
//	g_dd_{11,12,13,22,23,33}			# $g_{ij}$
//	K_dd_{11,12,13,22,23,33}			# $K_{ij}$
//	partial_d_g_dd_{1,2,3}{11,12,13,22,23,33}	# $\partial_k g_{ij}$
//	## computed at the nominal grid points
//	Theta						# $\Theta = \Theta(h)$
//
// Arguments:
// ps_ptr --> The patch system, or == NULL to do (only) a dummy computation
//	      in which only the parameter-table setup and a dummy geometry
//	      interpolator call are done, the latter with the number of
//	      interpolation points is set to 0 and all the output array
//	      pointers set to NULL.
// add_to_expansion = A real number which is added to the computed expansion
//		      at each grid point.
// initial_flag = true if this is the first evaluation of  expansion()
//		       for this horizon,
//		  false otherwise;
//		  this is used (only) to select which elements of  error_info
//		  are relevant
// Jacobian_flag = true to compute the Jacobian coefficients,
//		   false to skip this.
// print_msg_flag = true to print status messages,
//		    false to skip this.
// Theta_norms_ptr = (out) If this pointer is non-NULL, the norm object it
//			   points to is updated with all the Theta values
//			   in the grid.  This norm object can then be used
//			   to compute various (gridwise) norms of Theta.
//
// Results:
// This function returns a status code indicating whether the computation
// succeeded or failed, and if the latter, what caused the failure.
//
enum expansion_status Horizon::
  expansion(patch_system* ps_ptr,
            const struct what_to_compute& compute_info,
	    const struct cactus_grid_info& cgi,
	    const struct geometry_info& gi,
	    const struct error_info& error_info, bool initial_flag,
	    bool Jacobian_flag  /*= false*/,
	    bool print_msg_flag  /*= false*/ ,
	    jtutil::norm<fp>* Theta_norms_ptr /* = NULL */,
	    jtutil::norm<fp>* expansion_Theta_norms_ptr /* = NULL */,
	    jtutil::norm<fp>* inner_expansion_Theta_norms_ptr /* = NULL */,
	    jtutil::norm<fp>* product_expansion_Theta_norms_ptr /* = NULL */,
	    jtutil::norm<fp>* mean_curvature_Theta_norms_ptr /* = NULL */)
{
  const bool active_flag = (ps_ptr != NULL);
  if (print_msg_flag)
    then CCTK_VInfo(CCTK_THORNSTRING,
                    "   %s expansion",
                    active_flag ? "" : "dummy ");

  if (active_flag)
    then {
      //
      // normal computation
      //

      // fill in values of all ghosted gridfns in ghost zones
      ps_ptr->synchronize();

      if (gi.check_that_h_is_finite && !h_is_finite(*ps_ptr,
                                                    error_info, initial_flag,
                                                    print_msg_flag))
        then return expansion_failure__surface_nonfinite;
      // *** ERROR RETURN ***

      // set up xyz positions of grid points
      setup_xyz_posns(*ps_ptr, print_msg_flag);
    }

  // compute the "geometry" g_ij, K_ij, and partial_k g_ij
  // if (gi.hardwire_Schwarzschild_EF_geometry)
  //   then {
  //     if (active_flag)
  //       then Schwarzschild_EF_geometry(*ps_ptr,
  //                                      gi,
  //                                      print_msg_flag);
  //   }
  {
    // this is the only function we call unconditionally; it looks at
    // ps_ptr (non-NULL vs NULL) to choose a normal vs dummy computation
    const enum expansion_status status
      = interpolate_geometry(ps_ptr,
                             cgi, gi,
                             error_info, initial_flag,
                             print_msg_flag);
    if (status != expansion_success)
      then return status;				// *** ERROR RETURN ***
    // if (active_flag && cgi.use_Cactus_conformal_metric)
    //   then convert_conformal_to_physical(*ps_ptr,
    //                                      print_msg_flag);
  }

  if (active_flag)
    then {


      if (gi.check_that_geometry_is_finite
          && !geometry_is_finite(*ps_ptr,
                                 error_info, initial_flag,
                                 print_msg_flag))
        then return expansion_failure__geometry_nonfinite;
      // *** ERROR RETURN ***

      // Ensure that there is a norm object
      const bool want_norms = Theta_norms_ptr;
      jtutil::norm<fp> norms;
      if (compute_info.surface_selection != selection_definition)
        then if (! Theta_norms_ptr) Theta_norms_ptr = &norms;


      // compute remaining gridfns --> $\Theta$
      // and optionally also the Jacobian coefficients
      // by algebraic ops and angular finite differencing
      what_to_compute this_compute_info (compute_info);
      this_compute_info.surface_selection = selection_definition;
      if (!compute_Theta(*ps_ptr, this_compute_info,
                         Jacobian_flag, Theta_norms_ptr,
                         expansion_Theta_norms_ptr,
                         inner_expansion_Theta_norms_ptr,
                         product_expansion_Theta_norms_ptr,
                         mean_curvature_Theta_norms_ptr,
                         error_info, initial_flag,
                         print_msg_flag))
        then return expansion_failure__gij_not_positive_definite;
      // *** ERROR RETURN ***

      if (compute_info.surface_selection != selection_definition) {
        //
        // Apply correction to find a surface by its areal radius
        //
        // get mean expansion
        fp mean_expansion;
        fp areal_radius;
        switch (compute_info.surface_selection) {
        case selection_mean_coordinate_radius: {
          const int np = ps_ptr->N_grid_points();
          fp sum_expansion = 0;
          fp sum_radius = 0;
          for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
            patch& p = ps_ptr->ith_patch(pn);
            for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
              for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
                sum_expansion += p.gridfn(gfns::gfn__Theta, irho,isigma);
                sum_radius += p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
              }
            }
          }
          mean_expansion = sum_expansion / np;
          areal_radius = sum_radius / np;
          break;
        }
        case selection_areal_radius: {
          // get surface area
          const fp area = ps_ptr->integrate_gridfn
            (gfns::gfn__one, true, true, true,
             gfns::gfn__h,
             gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
             gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
             gfns::gfn__g_dd_33,
             patch::integration_method__automatic_choice);
          mean_expansion = Theta_norms_ptr->mean();
          areal_radius = sqrt(area / (4.0*PI));
          break;
        }
        case selection_expansion_mean_coordinate_radius: {
          const int np = ps_ptr->N_grid_points();
          fp sum_expansion = 0;
          fp sum_radius = 0;
          for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
            patch& p = ps_ptr->ith_patch(pn);
            for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
              for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
                sum_expansion += p.gridfn(gfns::gfn__Theta, irho,isigma);
                sum_radius += p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
              }
            }
          }
          mean_expansion = sum_expansion / np;
          areal_radius = mean_expansion * sum_radius / np;
          break;
        }
        case selection_expansion_areal_radius: {
          // get surface area
          const fp area = ps_ptr->integrate_gridfn
            (gfns::gfn__one, true, true, true,
             gfns::gfn__h,
             gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
             gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
             gfns::gfn__g_dd_33,
             patch::integration_method__automatic_choice);
          mean_expansion = Theta_norms_ptr->mean();
          areal_radius = mean_expansion * sqrt(area / (4.0*PI));
          break;
        }
        default:
          assert (0);
        } // switch areal_radius_definition
          
        if (! ps_ptr->N_additional_points()) {
          // calculate correction
          const fp correction
            = (- mean_expansion
               + areal_radius - compute_info.desired_value);
          // apply correction
          ps_ptr->add_to_gridfn(-correction, gfns::gfn__Theta);
        } else {
          const int np = ps_ptr->N_grid_points();
          const int gnp = ps_ptr->ghosted_N_grid_points();
          // apply correction
          const fp correction
            = ps_ptr->ghosted_gridfn_data(gfns::gfn__h)[gnp];
          ps_ptr->add_to_gridfn(-correction, gfns::gfn__Theta);
          ps_ptr->gridfn_data(gfns::gfn__Theta)[np]
            = (mean_expansion - correction
               - areal_radius + compute_info.desired_value);
        }
        if (want_norms) {
          // recalculate norms
          ps_ptr->gridfn_norms (gfns::gfn__Theta, *Theta_norms_ptr);
        } // if want_norms
      } else {
        //
        // do not apply correction
        //
        if (ps_ptr->N_additional_points()) {
          const int np = ps_ptr->N_grid_points();
          const int gnp = ps_ptr->ghosted_N_grid_points();
          ps_ptr->gridfn_data(gfns::gfn__Theta)[np]
            = ps_ptr->ghosted_gridfn_data(gfns::gfn__h)[gnp];
        }
      } // if use-areal-radius
    }

  return expansion_success;				// *** NORMAL RETURN ***
}

bool Horizon::h_is_finite(patch_system& ps,
		 const struct error_info& error_info, bool initial_flag,
		 bool print_msg_flag)
{
  if (print_msg_flag)
    then CCTK_VInfo(CCTK_THORNSTRING, "      checking that h is finite");

#ifdef HAVE_FINITE
  for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
  {
    patch& p = ps.ith_patch(pn);

    for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
    {
      for (int isigma = p.min_isigma() ;
           isigma <= p.max_isigma() ;
           ++isigma)
      {
        const fp h = p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
        if (!finite(h))
          then {
            const fp rho = p.rho_of_irho(irho);
            const fp sigma = p.sigma_of_isigma(isigma);
            const fp drho   = jtutil::degrees_of_radians(rho);
            const fp dsigma = jtutil::degrees_of_radians(sigma);
            CCTK_VWarn(error_info.warn_level__nonfinite_geometry,
                       __LINE__, __FILE__, CCTK_THORNSTRING,
                       "\n"
                       "   h=%g isn't finite!\n"
                       "   %s patch (rho,sigma)=(%g,%g) (drho,dsigma)=(%g,%g)\n"
                       ,
                       double(h),
                       p.name(), double(rho), double(sigma),
                       double(drho), double(dsigma));
            return false;			// *** found a NaN ***
          }
      }
    }
  }
  return true;					// *** all values finite ***
#else
  CCTK_VWarn(error_info.warn_level__skipping_finite_check,
             __LINE__, __FILE__, CCTK_THORNSTRING,
             "      no finite() fn ==>  skipping is-h-finite check");
  return true;					// *** no check possible ***
#endif
}


bool Horizon::geometry_is_finite(patch_system& ps,
			const struct error_info& error_info, bool initial_flag,
			bool print_msg_flag)
{
  if (print_msg_flag)
    then CCTK_VInfo(CCTK_THORNSTRING, "      checking that geometry is finite");

#ifdef HAVE_FINITE
  for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
  {
    patch& p = ps.ith_patch(pn);

    for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
    {
      for (int isigma = p.min_isigma() ;
           isigma <= p.max_isigma() ;
           ++isigma)
      {
	const fp g_dd_11 = p.gridfn(gfns::gfn__g_dd_11, irho,isigma);
	const fp g_dd_12 = p.gridfn(gfns::gfn__g_dd_12, irho,isigma);
	const fp g_dd_13 = p.gridfn(gfns::gfn__g_dd_13, irho,isigma);
	const fp g_dd_22 = p.gridfn(gfns::gfn__g_dd_22, irho,isigma);
	const fp g_dd_23 = p.gridfn(gfns::gfn__g_dd_23, irho,isigma);
	const fp g_dd_33 = p.gridfn(gfns::gfn__g_dd_33, irho,isigma);

	const fp K_dd_11 = p.gridfn(gfns::gfn__K_dd_11, irho,isigma);
	const fp K_dd_12 = p.gridfn(gfns::gfn__K_dd_12, irho,isigma);
	const fp K_dd_13 = p.gridfn(gfns::gfn__K_dd_13, irho,isigma);
	const fp K_dd_22 = p.gridfn(gfns::gfn__K_dd_22, irho,isigma);
	const fp K_dd_23 = p.gridfn(gfns::gfn__K_dd_23, irho,isigma);
	const fp K_dd_33 = p.gridfn(gfns::gfn__K_dd_33, irho,isigma);

	const fp partial_d_g_dd_111
          = p.gridfn(gfns::gfn__partial_d_g_dd_111, irho,isigma);
	const fp partial_d_g_dd_112
          = p.gridfn(gfns::gfn__partial_d_g_dd_112, irho,isigma);
	const fp partial_d_g_dd_113
          = p.gridfn(gfns::gfn__partial_d_g_dd_113, irho,isigma);
	const fp partial_d_g_dd_122
          = p.gridfn(gfns::gfn__partial_d_g_dd_122, irho,isigma);
	const fp partial_d_g_dd_123
          = p.gridfn(gfns::gfn__partial_d_g_dd_123, irho,isigma);
	const fp partial_d_g_dd_133
          = p.gridfn(gfns::gfn__partial_d_g_dd_133, irho,isigma);
	const fp partial_d_g_dd_211
          = p.gridfn(gfns::gfn__partial_d_g_dd_211, irho,isigma);
	const fp partial_d_g_dd_212
          = p.gridfn(gfns::gfn__partial_d_g_dd_212, irho,isigma);
	const fp partial_d_g_dd_213
          = p.gridfn(gfns::gfn__partial_d_g_dd_213, irho,isigma);
	const fp partial_d_g_dd_222
          = p.gridfn(gfns::gfn__partial_d_g_dd_222, irho,isigma);
	const fp partial_d_g_dd_223
          = p.gridfn(gfns::gfn__partial_d_g_dd_223, irho,isigma);
	const fp partial_d_g_dd_233
          = p.gridfn(gfns::gfn__partial_d_g_dd_233, irho,isigma);
	const fp partial_d_g_dd_311
          = p.gridfn(gfns::gfn__partial_d_g_dd_311, irho,isigma);
	const fp partial_d_g_dd_312
          = p.gridfn(gfns::gfn__partial_d_g_dd_312, irho,isigma);
	const fp partial_d_g_dd_313
          = p.gridfn(gfns::gfn__partial_d_g_dd_313, irho,isigma);
	const fp partial_d_g_dd_322
          = p.gridfn(gfns::gfn__partial_d_g_dd_322, irho,isigma);
	const fp partial_d_g_dd_323
          = p.gridfn(gfns::gfn__partial_d_g_dd_323, irho,isigma);
	const fp partial_d_g_dd_333
          = p.gridfn(gfns::gfn__partial_d_g_dd_333, irho,isigma);

	if (    !finite(g_dd_11) || !finite(g_dd_12) || !finite(g_dd_13)
                || !finite(g_dd_22) || !finite(g_dd_23)
                || !finite(g_dd_33)
                || !finite(K_dd_11) || !finite(K_dd_12) || !finite(K_dd_13)
                || !finite(K_dd_22) || !finite(K_dd_23)
                || !finite(K_dd_33)
                || !finite(partial_d_g_dd_111)
                || !finite(partial_d_g_dd_112)
                || !finite(partial_d_g_dd_113)
                || !finite(partial_d_g_dd_122)
                || !finite(partial_d_g_dd_123)
                || !finite(partial_d_g_dd_133)
                || !finite(partial_d_g_dd_211)
                || !finite(partial_d_g_dd_212)
                || !finite(partial_d_g_dd_213)
                || !finite(partial_d_g_dd_222)
                || !finite(partial_d_g_dd_223)
                || !finite(partial_d_g_dd_233)
                || !finite(partial_d_g_dd_311)
                || !finite(partial_d_g_dd_312)
                || !finite(partial_d_g_dd_313)
                || !finite(partial_d_g_dd_322)
                || !finite(partial_d_g_dd_323)
                || !finite(partial_d_g_dd_333)    )
          then {
            const fp h = p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
            const fp rho = p.rho_of_irho(irho);
            const fp sigma = p.sigma_of_isigma(isigma);
            const fp drho   = jtutil::degrees_of_radians(rho);
            const fp dsigma = jtutil::degrees_of_radians(sigma);
            fp local_x, local_y, local_z;
            p.xyz_of_r_rho_sigma(h,rho,sigma, local_x,local_y,local_z);
            const fp global_x = ps.origin_x() + local_x;
            const fp global_y = ps.origin_y() + local_y;
            const fp global_z = ps.origin_z() + local_z;
            return false;			// *** found a NaN ***
          }
      }
    }
  }
  return true;						// *** no NaNs found ***
#else
  CCTK_VWarn(error_info.warn_level__skipping_finite_check,
             __LINE__, __FILE__, CCTK_THORNSTRING,
             "      no finite() ==>  skipping is-geometry-finite check");
  return true;					// *** no check possible ***
#endif
}

bool Horizon::compute_Theta(patch_system& ps, const struct what_to_compute& compute_info,
		   bool Jacobian_flag,
                   jtutil::norm<fp>* Theta_norms_ptr,
                   jtutil::norm<fp>* expansion_Theta_norms_ptr,
                   jtutil::norm<fp>* inner_expansion_Theta_norms_ptr,
                   jtutil::norm<fp>* product_expansion_Theta_norms_ptr,
                   jtutil::norm<fp>* mean_curvature_Theta_norms_ptr,
		   const struct error_info& error_info, bool initial_flag,
		   bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING, "      computing Theta(h)");

#if 0
 fp mean_radius, areal_radius;
 switch (compute_info.surface_modification) {
 case modification_none:
 case modification_radius:
 case modification_radius2:
   // do nothing
   break;
 case modification_mean_radius: {
   // get average coordinate radius
   const int np = ps.N_grid_points();
   fp sum_radius = 0;
   for (int pn = 0; pn < ps.N_patches(); ++pn) {
     patch& p = ps.ith_patch(pn);
     for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
       for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
         sum_radius += p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
       }
     }
   }
   mean_radius = sum_radius / np;
   break;
 }
 case modification_areal_radius: {
   // get surface area
   const fp area = ps.integrate_gridfn
     (gfns::gfn__one, true, true, true,
      gfns::gfn__h,
      gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                          gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                              gfns::gfn__g_dd_33,
      patch::integration_method__automatic_choice);
   areal_radius = sqrt(area / (4.0*PI));
   break;
 }
 default:
   assert (0);
 }
#endif

	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		//
		// compute the X_ud and X_udd derivative coefficients
		// ... n.b. this uses the *local* (x,y,z) coordinates
		//
		const fp r = p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
		const fp rho = p.rho_of_irho(irho);
		const fp sigma = p.sigma_of_isigma(isigma);
		fp xx, yy, zz;
		p.xyz_of_r_rho_sigma(r,rho,sigma, xx, yy, zz);
		#ifdef COMPUTE_THETA_DEBUG
		if (    ((irho == 0) && (isigma == 0))
		     || ((irho == 5) && (isigma == 5))    )
		   then {
			printf("AHFinderDirect:: computing at %s patch (%d,%d) r=%g ==> local xyz=(%g,%g,%g)\n",
			       p.name(), irho,isigma, double(r),
			       double(xx), double(yy), double(zz));
			printf("AHFinderDirect:: got g_dd_11=%g g_dd_12=%g g_dd_13=%g\n",
			       p.gridfn(gfns::gfn__g_dd_11, irho,isigma),
			       p.gridfn(gfns::gfn__g_dd_12, irho,isigma),
			       p.gridfn(gfns::gfn__g_dd_13, irho,isigma));
			printf("                                g_dd_22=%g g_dd_23=%g\n",
			       p.gridfn(gfns::gfn__g_dd_22, irho,isigma),
			       p.gridfn(gfns::gfn__g_dd_23, irho,isigma));
			printf("                                           g_dd_33=%g\n",
			       p.gridfn(gfns::gfn__g_dd_33, irho,isigma));
			fflush(stdout);
			}
		#endif

		// 1st derivative coefficients X_ud
		const fp X_ud_11 = p.partial_rho_wrt_x(xx, yy, zz);
		const fp X_ud_12 = p.partial_rho_wrt_y(xx, yy, zz);
		const fp X_ud_13 = p.partial_rho_wrt_z(xx, yy, zz);
		const fp X_ud_21 = p.partial_sigma_wrt_x(xx, yy, zz);
		const fp X_ud_22 = p.partial_sigma_wrt_y(xx, yy, zz);
		const fp X_ud_23 = p.partial_sigma_wrt_z(xx, yy, zz);

		// 2nd derivative coefficient gridfns X_udd
		const fp X_udd_111 = p.partial2_rho_wrt_xx(xx, yy, zz);
		const fp X_udd_112 = p.partial2_rho_wrt_xy(xx, yy, zz);
		const fp X_udd_113 = p.partial2_rho_wrt_xz(xx, yy, zz);
		const fp X_udd_122 = p.partial2_rho_wrt_yy(xx, yy, zz);
		const fp X_udd_123 = p.partial2_rho_wrt_yz(xx, yy, zz);
		const fp X_udd_133 = p.partial2_rho_wrt_zz(xx, yy, zz);
		const fp X_udd_211 = p.partial2_sigma_wrt_xx(xx, yy, zz);
		const fp X_udd_212 = p.partial2_sigma_wrt_xy(xx, yy, zz);
		const fp X_udd_213 = p.partial2_sigma_wrt_xz(xx, yy, zz);
		const fp X_udd_222 = p.partial2_sigma_wrt_yy(xx, yy, zz);
		const fp X_udd_223 = p.partial2_sigma_wrt_yz(xx, yy, zz);
		const fp X_udd_233 = p.partial2_sigma_wrt_zz(xx, yy, zz);

		//
		// "call" the Maple-generated code
		// ... each cg/*.c file has a separate set of temp variables,
		//     and so must be inside its own set of { } braces
		//

		// gridfn #defines
		#include "gr/cg.hh"

		  {
		// g_uu
		#include "gr.cg/inverse_metric.c"
		  }

		  {
		// K, K_uu
		#include "gr.cg/extrinsic_curvature_trace_raise.c"
		  }

		  {
		// partial_d_g_uu
		#include "gr.cg/inverse_metric_gradient.c"
		  }

		  {
		// partial_d_ln_sqrt_g
		#include "gr.cg/metric_det_gradient.c"
		  }

		  {
		// Theta_A, Theta_B, Theta_C, Theta_D
		#include "gr.cg/expansion.c"
		  }

		if (Theta_D <= 0)
		   then {
const int warn_level
  = initial_flag ? error_info.warn_level__gij_not_positive_definite__initial
		 : error_info.warn_level__gij_not_positive_definite__subsequent;
CCTK_VWarn(warn_level, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   compute_Theta(): Theta_D = $g^{ij} s_i s_j$ = %g <= 0\n"
"                    at %s patch rho=%g sigma=%g!\n"
"                    (i.e. the interpolated g_ij isn't positive definite)",
	   double(Theta_D),
	   p.name(), double(rho), double(sigma));
			return false;			// *** ERROR RETURN ***
			}

                assert (compute_info.surface_selection == selection_definition);

		// compute H via equation (14) of my 1996 horizon finding paper
		const fp sqrt_Theta_D = sqrt(Theta_D);

		const fp Theta_X = + Theta_A/(Theta_D*sqrt_Theta_D)
                                   + Theta_B/sqrt_Theta_D;
                const fp Theta_Y = + Theta_C/Theta_D
                                   - K;

#define mean_curvature	p.gridfn(gfns::gfn__mean_curvature, irho,isigma)
		mean_curvature = Theta_X;
#undef mean_curvature

                switch (compute_info.surface_definition) {
                case definition_expansion:
                  Theta = + Theta_X + Theta_Y;
                  break;
                case definition_inner_expansion:
                  Theta = - Theta_X + Theta_Y;
                  break;
                case definition_mean_curvature:
                  Theta = + Theta_X;
                  break;
                case definition_expansion_product:
                  Theta = (+ Theta_X + Theta_Y) * (- Theta_X + Theta_Y);
                  break;
                default:
                  assert (0);
                }

                switch (compute_info.surface_modification) {
                case modification_none:
                  // do nothing
                  break;
                case modification_radius:
                  // multiply by radius
                  Theta *= r;
                  break;
                case modification_radius2:
                  // multiply by radius^2
                  Theta *= pw2(r);
                  break;
#if 0
                case modification_mean_radius:
                  // multiply by average coordinate radius
                  Theta *= mean_radius;
                  break;
                case modification_areal_radius:
                  // multiply by areal radius
                  Theta *= areal_radius;
                  break;
#endif
                default:
                  assert (0);
                }

                Theta -= compute_info.desired_value;
                
		// update running norms of Theta(h) function
		if (Theta_norms_ptr != NULL)
		   then Theta_norms_ptr->data(Theta);

		if (expansion_Theta_norms_ptr != NULL)
                   then expansion_Theta_norms_ptr->data(+ Theta_X + Theta_Y);

		if (inner_expansion_Theta_norms_ptr != NULL)
                   then inner_expansion_Theta_norms_ptr->data(- Theta_X + Theta_Y);

		if (product_expansion_Theta_norms_ptr != NULL)
                   then product_expansion_Theta_norms_ptr->data((+ Theta_X + Theta_Y) * (- Theta_X + Theta_Y));

		if (mean_curvature_Theta_norms_ptr != NULL)
                   then mean_curvature_Theta_norms_ptr->data(+ Theta_X);

                fp partial_Theta_X_wrt_partial_d_h_1;
                fp partial_Theta_X_wrt_partial_d_h_2;
                fp partial_Theta_X_wrt_partial_dd_h_11;
                fp partial_Theta_X_wrt_partial_dd_h_12;
                fp partial_Theta_X_wrt_partial_dd_h_22;
                fp partial_Theta_Y_wrt_partial_d_h_1;
                fp partial_Theta_Y_wrt_partial_d_h_2;
                fp partial_Theta_Y_wrt_partial_dd_h_11;
                fp partial_Theta_Y_wrt_partial_dd_h_12;
                fp partial_Theta_Y_wrt_partial_dd_h_22;

		if (Jacobian_flag)
		   then {
			// partial_Theta_wrt_partial_d_h,
			// partial_Theta_wrt_partial_dd_h
			#include "gr.cg/expansion_Jacobian.c"
			}

		if (Jacobian_flag) {
                  switch (compute_info.surface_definition) {
                    
                  case definition_expansion:
                    partial_Theta_wrt_partial_d_h_1
                      = (+ partial_Theta_X_wrt_partial_d_h_1
                         + partial_Theta_Y_wrt_partial_d_h_1);
                    partial_Theta_wrt_partial_d_h_2
                      = (+ partial_Theta_X_wrt_partial_d_h_2
                         + partial_Theta_Y_wrt_partial_d_h_2);
                    partial_Theta_wrt_partial_dd_h_11
                      = (+ partial_Theta_X_wrt_partial_dd_h_11
                         + partial_Theta_Y_wrt_partial_dd_h_11);
                    partial_Theta_wrt_partial_dd_h_12
                      = (+ partial_Theta_X_wrt_partial_dd_h_12
                         + partial_Theta_Y_wrt_partial_dd_h_12);
                    partial_Theta_wrt_partial_dd_h_22
                      = (+ partial_Theta_X_wrt_partial_dd_h_22
                         + partial_Theta_Y_wrt_partial_dd_h_22);
                    break;
                    
                  case definition_inner_expansion:
                    partial_Theta_wrt_partial_d_h_1
                      = (- partial_Theta_X_wrt_partial_d_h_1
                         + partial_Theta_Y_wrt_partial_d_h_1);
                    partial_Theta_wrt_partial_d_h_2
                      = (- partial_Theta_X_wrt_partial_d_h_2
                         + partial_Theta_Y_wrt_partial_d_h_2);
                    partial_Theta_wrt_partial_dd_h_11
                      = (- partial_Theta_X_wrt_partial_dd_h_11
                         + partial_Theta_Y_wrt_partial_dd_h_11);
                    partial_Theta_wrt_partial_dd_h_12
                      = (- partial_Theta_X_wrt_partial_dd_h_12
                         + partial_Theta_Y_wrt_partial_dd_h_12);
                    partial_Theta_wrt_partial_dd_h_22
                      = (- partial_Theta_X_wrt_partial_dd_h_22
                         + partial_Theta_Y_wrt_partial_dd_h_22);
                    break;
                    
                  case definition_mean_curvature:
                    partial_Theta_wrt_partial_d_h_1
                      = + partial_Theta_X_wrt_partial_d_h_1;
                    partial_Theta_wrt_partial_d_h_2
                      = + partial_Theta_X_wrt_partial_d_h_2;
                    partial_Theta_wrt_partial_dd_h_11
                      = + partial_Theta_X_wrt_partial_dd_h_11;
                    partial_Theta_wrt_partial_dd_h_12
                      = + partial_Theta_X_wrt_partial_dd_h_12;
                    partial_Theta_wrt_partial_dd_h_22
                      = + partial_Theta_X_wrt_partial_dd_h_22;
                    break;
                    
                  case definition_expansion_product: {
#define f(x,y,dx,dy)  (- x*x + y*y)
#define df(x,y,dx,dy) (- 2*x*dx + 2*y*dy)
                    partial_Theta_wrt_partial_d_h_1   = df(Theta_X, Theta_Y, partial_Theta_X_wrt_partial_d_h_1  , partial_Theta_Y_wrt_partial_d_h_1  );
                    partial_Theta_wrt_partial_d_h_2   = df(Theta_X, Theta_Y, partial_Theta_X_wrt_partial_d_h_2  , partial_Theta_Y_wrt_partial_d_h_2  );
                    partial_Theta_wrt_partial_dd_h_11 = df(Theta_X, Theta_Y, partial_Theta_X_wrt_partial_dd_h_11, partial_Theta_Y_wrt_partial_dd_h_11);
                    partial_Theta_wrt_partial_dd_h_12 = df(Theta_X, Theta_Y, partial_Theta_X_wrt_partial_dd_h_12, partial_Theta_Y_wrt_partial_dd_h_12);
                    partial_Theta_wrt_partial_dd_h_22 = df(Theta_X, Theta_Y, partial_Theta_X_wrt_partial_dd_h_22, partial_Theta_Y_wrt_partial_dd_h_22);
#undef f
#undef df
                    break;
                  }
                    
                  default:
                    assert (0);
                  }
                }

		}
		}
	}
		#include "gr/uncg.hh"
return true;						// *** NORMAL RETURN ***
}

//******************************************************************************

//
// If ps_ptr != NULL and Jac_ptr != NULL, this function computes the
// Jacobian matrix J[Theta(h)] of the expansion Theta(h).  We assume
// that Theta(h) has already been computed.
//
// If ps_ptr == NULL and Jac_ptr == NULL, this function does a dummy
// computation, in which only any expansion() (and hence geometry
// interpolator) calls are done, these with the number of interpolation
// points set to 0 and all the output array pointers set to NULL.
//
// It's illegal for one but not both of ps_ptr and Jac_ptr to be NULL.
//
// Only some values of  Jacobian_info.Jacobian_compute_method  support
// the dummy computation.
//
// Arguments:
// ps_ptr --> The patch system, or == NULL to do (only) a dummy computation.
// Jac_ptr --> The Jacobian, or == NULL to do (only) a dummy computation.
// add_to_expansion = A real number to add to the expansion.
//
// Results:
// This function returns a status code indicating whether the computation
// succeeded or failed, and if the latter, what caused the failure.
//
enum expansion_status
  Horizon::expansion_Jacobian(patch_system* ps_ptr, Jacobian* Jac_ptr,
                     const struct what_to_compute& compute_info,
		     const struct cactus_grid_info& cgi,
		     const struct geometry_info& gi,
		     const struct Jacobian_info& Jacobian_info,
		     const struct error_info& error_info, bool initial_flag,
		     bool print_msg_flag /* = false */)
{
const bool active_flag = (ps_ptr != NULL) && (Jac_ptr != NULL);
enum expansion_status status;

switch	(Jacobian_info.Jacobian_compute_method)
	{
case Jacobian__numerical_perturbation:
	if (active_flag)
	   then {
		status = expansion_Jacobian_NP(*ps_ptr, *Jac_ptr,
					       compute_info,
					       cgi, gi, Jacobian_info,
					       error_info, initial_flag,
					       print_msg_flag);
		if (status != expansion_success)
		   then return status;			// *** ERROR RETURN ***
		break;
		}
	   else error_exit(ERROR_EXIT,
"***** expansion_Jacobian():\n"
                           "        dummy computation isn't supported for\n"
"        Jacobian_compute_method = \"numerical perturbation\"!\n");
								/*NOTREACHED*/

case Jacobian__symbolic_diff:
	error_exit(ERROR_EXIT,
"***** expansion_Jacobian():\n"
                   "        Jacobian_compute_method == \"symbolic differentiation\"\n"
"        isn't implemented (yet)!\n");				/*NOTREACHED*/

case Jacobian__symbolic_diff_with_FD_dr:
	if (active_flag)
	   then expansion_Jacobian_partial_SD(*ps_ptr, *Jac_ptr,
					      cgi, gi, Jacobian_info,
					      print_msg_flag);
	// this function looks at ps_ptr and Jac_ptr (non-NULL vs NULL)
	// to choose a normal vs dummy computation
	  {
	status = expansion_Jacobian_dr_FD(ps_ptr, Jac_ptr,
                                          compute_info,
					  cgi, gi, Jacobian_info,
					  error_info, initial_flag,
					  print_msg_flag);
	if (status != expansion_success)
	   then return status;				// *** ERROR RETURN ***
	  }
	break;

default:
	error_exit(PANIC_EXIT,
"***** expansion_Jacobian():\n"
"        unknown Jacobian_info.Jacobian_compute_method=(int)%d!\n"
"        (this should never happen!)\n"
,
		   int(Jacobian_info.Jacobian_compute_method));	/*NOTREACHED*/
	}

 if (active_flag) {

   switch (compute_info.surface_modification) {

   case modification_none:
     // do nothing
     break;

   case modification_radius: {
     // multiply with the coordinate radius
     // H_{(r)i} = H_i h_i
     // J_{(r)ij} = J_{ij} h_i + H_i \delta_{ij}
     const int np = ps_ptr->N_grid_points();
     for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
       patch& p = ps_ptr->ith_patch(pn);
       for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
         for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
           const int i = ps_ptr->gpn_of_patch_irho_isigma(p, irho,isigma);
           const fp radius = p.ghosted_gridfn(gfns::gfn__h, irho, isigma);
           for (int j=0; j<np; ++j) {
             if (Jac_ptr->is_explicitly_stored (i, j)) {
               const fp val = Jac_ptr->element (i, j);
               Jac_ptr->set_element (i, j, val * radius);
             }
           }
           const fp Theta = (p.gridfn(gfns::gfn__Theta, irho, isigma)
                             + compute_info.desired_value) / radius;
           Jac_ptr->sum_into_element (i, i, Theta);
         }
       }
     }
     break;
   }

   case modification_radius2: {
     // multiply with the square of the coordinate radius
     // H_{(r2)i} = H_i h_i^2
     // J_{(r2)ij} = J_{ij} h_i^2 + 2 H_i h_i \delta_{ij}
     const int np = ps_ptr->N_grid_points();
     for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
       patch& p = ps_ptr->ith_patch(pn);
       for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
         for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
           const int i = ps_ptr->gpn_of_patch_irho_isigma(p, irho,isigma);
           const fp radius = p.ghosted_gridfn(gfns::gfn__h, irho, isigma);
           const fp radius2 = radius * radius;
           for (int j=0; j<np; ++j) {
             if (Jac_ptr->is_explicitly_stored (i, j)) {
               const fp val = Jac_ptr->element (i, j);
               Jac_ptr->set_element (i, j, val * radius2);
             }
           }
           const fp Theta = (p.gridfn(gfns::gfn__Theta, irho, isigma)
                             + compute_info.desired_value) / radius2;
           Jac_ptr->sum_into_element (i, i, 2 * Theta * radius);
         }
       }
     }
     break;
   }

#if 0
   case modification_mean_radius: {
     // multiply with the average coordinate radius
     // H_{(\bar r)i} = H_i \bar r
     // J_{(\bar r)ij} = J_{ij} \bar r + H_i / N
     // calculate average coordinate radius
     const int np = ps_ptr->N_grid_points();
     fp sum_radius = 0;
     for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
       patch& p = ps_ptr->ith_patch(pn);
       for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
         for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
           sum_radius += p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
         }
       }
     }
     mean_radius = sum_radius / np;
     // correct Jacobian
     const int np = ps_ptr->N_grid_points();
     for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
       patch& p = ps_ptr->ith_patch(pn);
       for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
         for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
           const int i = ps_ptr->gpn_of_patch_irho_isigma(p, irho,isigma);
           for (int j=0; j<np; ++j) {
             if (Jac_ptr->is_explicitly_stored (i, j)) {
               const fp val = Jac_ptr->element (i, j);
               Jac_ptr->set_element (i, j, val * mean_radius);
             }
           }
#error "unfinished"
           const fp Theta = (p.gridfn(gfns::gfn__Theta, irho, isigma)
                             + compute_info.desired_value) / areal_radius;
           const fp dRdh = 0.5 * areal_radius;
           Jac_ptr->sum_into_element (i, i, Theta * dRdh);
         }
       }
     }
     break;
   }

   case modification_areal_radius: {
     // multiply with the areal radius
     // H_{(R)i} = H_i R
     // J_{(R)ij} = J_{ij} R + H_i dR/dh_j
     // get surface area
     const fp area = ps_ptr->integrate_gridfn
       (gfns::gfn__one, true, true, true,
        gfns::gfn__h,
        gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                            gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                gfns::gfn__g_dd_33,
        patch::integration_method__automatic_choice);
     const fp areal_radius = sqrt(area / (4.0*PI));
     // correct Jacobian
     const int np = ps_ptr->N_grid_points();
     for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
       patch& p = ps_ptr->ith_patch(pn);
       for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
         for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
           const int i = ps_ptr->gpn_of_patch_irho_isigma(p, irho,isigma);
           for (int j=0; j<np; ++j) {
             if (Jac_ptr->is_explicitly_stored (i, j)) {
               const fp val = Jac_ptr->element (i, j);
               Jac_ptr->set_element (i, j, val * areal_radius);
             }
           }
           const fp Theta = (p.gridfn(gfns::gfn__Theta, irho, isigma)
                             + compute_info.desired_value) / areal_radius;
           const fp dRdh = 0.5 * areal_radius;
           Jac_ptr->sum_into_element (i, i, Theta * dRdh);
         }
       }
     }
     break;
   }
#endif

   default:
     assert (0);
   } // switch surface_modification
   
   if (ps_ptr->N_additional_points()) {
     switch (compute_info.surface_selection) {

     case selection_definition: {
       // we want nothing special
       const int np = ps_ptr->N_grid_points();
       for (int i=0; i<np; ++i) {
         Jac_ptr->set_element (i, np, 0.0);
       }
       for (int j=0; j<np; ++j) {
         Jac_ptr->set_element (np, j, 0.0);
       }
       Jac_ptr->set_element (np, np, 1.0);
       break;
     }

     case selection_mean_coordinate_radius: {
       // Jac_ptr->set_element (II, JJ, x) == dTheta(II)/dh(JJ)
       // \frac{\partial R}{\partial h_j} = 1 / N
       const int np = ps_ptr->N_grid_points();
       for (int i=0; i<np; ++i) {
         Jac_ptr->set_element (i, np, -1.0);
       }
       for (int j=0; j<np; ++j) {
         fp val = 0;
         for (int k=0; k<np; ++k) {
           val += Jac_ptr->element (k, j) / np;
         }
         val -= 1.0 / np;
         Jac_ptr->set_element (np, j, val);
       }
       Jac_ptr->set_element (np, np, -1.0);
       break;
     }

     case selection_areal_radius: {
       // \frac{\partial R_a}{\partial h_j}
       //    = \sqrt{1 / 16 \pi A} \sum_k \sqrt{q_k} dS_k
       // The "trapezoid" method is faster
//        const enum patch::integration_method method
//          = patch::integration_method__automatic_choice;
       const enum patch::integration_method method
         = patch::integration_method__trapezoid;
       const fp area = ps_ptr->integrate_gridfn
         (gfns::gfn__one, true, true, true,
          gfns::gfn__h,
          gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
			      gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
						  gfns::gfn__g_dd_33,
          method);
       const int np = ps_ptr->N_grid_points();
       for (int i=0; i<np; ++i) {
         Jac_ptr->set_element (i, np, -1.0);
       }
       for (int j=0; j<np; ++j) {
         fp val = 0;
         for (int k=0; k<np; ++k) {
           val += Jac_ptr->element (k, j) / np;
         }
         Jac_ptr->set_element (np, j, val);
       }
       for (int jpn = 0; jpn < ps_ptr->N_patches(); ++jpn) {
         patch& jp = ps_ptr->ith_patch(jpn);
         for (int jrho = jp.min_irho(); jrho <= jp.max_irho(); ++jrho) {
           for (int jsigma = jp.min_isigma(); jsigma <= jp.max_isigma(); ++jsigma) {
             const int j = ps_ptr->gpn_of_patch_irho_isigma(jp, jrho,jsigma);
             // const fp radius = jp.ghosted_gridfn(gfns::gfn__h, jrho, jsigma);
             const fp epsilon = Jacobian_info.perturbation_amplitude;
             fp val1, val2;
#if 0
             // Re-calculate all points
             // (this is slow, but it works)
             val1 = area;
             jp.ghosted_gridfn(gfns::gfn__h, jrho, jsigma) += epsilon;
             val2 = ps_ptr->integrate_gridfn
               (gfns::gfn__one, true, true, true,
                gfns::gfn__h,
                gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                                    gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                        gfns::gfn__g_dd_33,
                method);
             jp.ghosted_gridfn(gfns::gfn__h, jrho, jsigma) -= epsilon;
#else
             // Re-calculate all points with non-zero Jacobian entries
             jp.ghosted_gridfn(gfns::gfn__h, jrho, jsigma) -= epsilon/2;
             val1 = 0;
             for (int ipn = 0; ipn < ps_ptr->N_patches(); ++ipn) {
               patch& ip = ps_ptr->ith_patch(ipn);
               for (int irho = ip.min_irho(); irho <= ip.max_irho(); ++irho) {
                 for (int isigma = ip.min_isigma(); isigma <= ip.max_isigma(); ++isigma) {
                   const int i = ps_ptr->gpn_of_patch_irho_isigma(ip, irho,isigma);
                   if (Jac_ptr->is_explicitly_stored (i, j)) {
                     val1 += ps_ptr->integrate_gridpoint
                       (gfns::gfn__one,
                        gfns::gfn__h,
                        gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                                            gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                                gfns::gfn__g_dd_33,
                        method,
                        ipn, irho, isigma);
                   }
                 }
               }
             }
             jp.ghosted_gridfn(gfns::gfn__h, jrho, jsigma) += epsilon;
             val2 = 0;
             for (int ipn = 0; ipn < ps_ptr->N_patches(); ++ipn) {
               patch& ip = ps_ptr->ith_patch(ipn);
               for (int irho = ip.min_irho(); irho <= ip.max_irho(); ++irho) {
                 for (int isigma = ip.min_isigma(); isigma <= ip.max_isigma(); ++isigma) {
                   const int i = ps_ptr->gpn_of_patch_irho_isigma(ip, irho,isigma);
                   if (Jac_ptr->is_explicitly_stored (i, j)) {
                     val2 += ps_ptr->integrate_gridpoint
                       (gfns::gfn__one,
                        gfns::gfn__h,
                        gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
                                            gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
                                                                gfns::gfn__g_dd_33,
                        method,
                        ipn, irho, isigma);
                   }
                 }
               }
             }
             jp.ghosted_gridfn(gfns::gfn__h, jrho, jsigma) -= epsilon/2;
#endif
             const fp val = 1 / sqrt(16*PI*area) * ps_ptr->integrate_correction(true, true, true) * (val2 - val1) / epsilon;
             Jac_ptr->sum_into_element (np, j, -val);
           }
         }
       }
       Jac_ptr->set_element (np, np, -1.0);
       break;
     }

     case selection_expansion_mean_coordinate_radius: {
       // Jac_ptr->set_element (II, JJ, x) == dTheta(II)/dh(JJ)
       const int np = ps_ptr->N_grid_points();
       fp sum_expansion = 0;
       fp sum_radius = 0;
       for (int pn = 0; pn < ps_ptr->N_patches(); ++pn) {
         patch& p = ps_ptr->ith_patch(pn);
         for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho) {
           for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma) {
             sum_expansion += p.gridfn(gfns::gfn__Theta, irho,isigma);
             sum_radius += p.ghosted_gridfn(gfns::gfn__h, irho,isigma);
           }
         }
       }
       for (int i=0; i<np; ++i) {
         Jac_ptr->set_element (i, np, -1.0);
       }
       for (int j=0; j<np; ++j) {
         fp val = 0;
         for (int k=0; k<np; ++k) {
           val += (Jac_ptr->element (k, j) / np) * (1.0 - sum_radius / np);
         }
         val -= (sum_expansion / np) / np;
         Jac_ptr->set_element (np, j, val);
       }
       Jac_ptr->set_element (np, np, -1.0);
       break;
     }

     case selection_expansion_areal_radius: {
       CCTK_WARN (0, "selection_expansion_areal_radius not implemented");
       break;
     }

     default:
       assert (0);
     } // switch surface_selection
   } else {
     assert (compute_info.surface_selection == selection_definition);
   }
   
 } // if active

return expansion_success;				// *** NORMAL RETURN ***
}



//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function computes the Jacobian matrix of the expansion Theta(h)
// by numerical perturbation of the Theta(h) function.  The algorithm is
// as follows:
//
// we assume that Theta = Theta(h) has already been evaluated
// save_Theta = Theta
//	for each point (y,JJ)
//	{
//	const fp save_h_y = h at y;
//	h at y += perturbation_amplitude;
//	evaluate Theta(h) (silently)
//		for each point (x,II)
//		{
//		Jac(II,JJ) = (Theta(II) - save_Theta(II))
//			     / perturbation_amplitude;
//		}
//	h at y = save_h_y;
//	}
// Theta = save_Theta
//
// Inputs (angular gridfns, on ghosted grid):
//	h			# shape of trial surface
//	Theta			# Theta(h) assumed to already be computed
//
// Outputs:
//	The Jacobian matrix is stored in the Jacobian object Jac.
//	As implied by the above algorithm, it's traversed by columns.
//
// Results:
// This function returns a status code indicating whether the computation
// succeeded or failed, and if the latter, what caused the failure.
//

enum expansion_status
  Horizon::expansion_Jacobian_NP
        (patch_system& ps, Jacobian& Jac,
         const struct what_to_compute& compute_info,
	 const struct cactus_grid_info& cgi,
	 const struct geometry_info& gi,
	 const struct Jacobian_info& Jacobian_info,
	 const struct error_info& error_info, bool initial_flag,
	 bool print_msg_flag)
{
  if (print_msg_flag)
    then CCTK_VInfo(CCTK_THORNSTRING,
                    "   horizon Jacobian (numerical perturbation)");
  const fp epsilon = Jacobian_info.perturbation_amplitude;

  ps.gridfn_copy(gfns::gfn__Theta, gfns::gfn__save_Theta);
  ps.gridfn_copy(gfns::gfn__mean_curvature, gfns::gfn__save_mean_curvature);

  for (int ypn = 0 ; ypn < ps.N_patches() ; ++ypn)
  {
    patch& yp = ps.ith_patch(ypn);
    if (print_msg_flag)
      then CCTK_VInfo(CCTK_THORNSTRING,
                      "      perturbing in %s patch",
                      yp.name());

    for (int y_irho = yp.min_irho() ; y_irho <= yp.max_irho() ; ++y_irho)
    {
      for (int y_isigma = yp.min_isigma() ;
           y_isigma <= yp.max_isigma() ;
           ++y_isigma)
      {
	const int JJ = ps.gpn_of_patch_irho_isigma(yp, y_irho,y_isigma);

	const fp save_h_y = yp.ghosted_gridfn(gfns::gfn__h, y_irho,y_isigma);
	yp.ghosted_gridfn(gfns::gfn__h, y_irho,y_isigma) += epsilon;
	const
	  enum expansion_status status = expansion(&ps,
                                                   compute_info,
						   cgi, gi,
						   error_info, initial_flag);
	if (status != expansion_success)
          then return status;				// *** ERROR RETURN ***

        for (int xpn = 0 ; xpn < ps.N_patches() ; ++xpn)
        {
          patch& xp = ps.ith_patch(xpn);

          for (int x_irho = xp.min_irho() ;
               x_irho <= xp.max_irho() ;
               ++x_irho)
          {
            for (int x_isigma = xp.min_isigma() ;
                 x_isigma <= xp.max_isigma() ;
                 ++x_isigma)
            {
              const int II = ps.gpn_of_patch_irho_isigma(xp, x_irho,x_isigma);
              const fp old_Theta = xp.gridfn(gfns::gfn__save_Theta,
                                             x_irho,x_isigma);
              const fp new_Theta = xp.gridfn(gfns::gfn__Theta,
                                             x_irho,x_isigma);
              Jac.set_element(II,JJ, (new_Theta - old_Theta) / epsilon);
            }
          }
        }

	if (ps.N_additional_points())
        {
          const int np = ps.N_grid_points();
          Jac.set_element(np,JJ, 0.0);	// insert dummy value
        }

	yp.ghosted_gridfn(gfns::gfn__h, y_irho,y_isigma) = save_h_y;
      }
    }
  } 

  ps.gridfn_copy(gfns::gfn__save_Theta, gfns::gfn__Theta);
  ps.gridfn_copy(gfns::gfn__save_mean_curvature, gfns::gfn__mean_curvature);
  return expansion_success;				// *** NORMAL RETURN ***
}

//******************************************************************************

//
// This function computes the partial derivative terms in the Jacobian
// matrix of the expansion Theta(h), by symbolic differentiation from
// the Jacobian coefficient (angular) gridfns.  The Jacobian is traversed
// by rows, using equation (25) of my 1996 apparent horizon finding paper.
//
// Inputs (angular gridfns, on ghosted grid):
//	h			# shape of trial surface
//	Theta			# Theta(h) assumed to already be computed
//	partial_Theta_wrt_partial_d_h	# Jacobian coefficients
//	partial_Theta_wrt_partial_dd_h	# (also assumed to already be computed)
//
// Outputs:
//	The Jacobian matrix is stored in the Jacobian object Jac.

 void Horizon::expansion_Jacobian_partial_SD(patch_system& ps, Jacobian& Jac,
				   const struct cactus_grid_info& cgi,
				   const struct geometry_info& gi,
				   const struct Jacobian_info& Jacobian_info,
				   bool print_msg_flag)
{
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   horizon Jacobian: partial-deriv terms (symbolic diff)");

Jac.zero_matrix();
ps.compute_synchronize_Jacobian();

    for (int xpn = 0 ; xpn < ps.N_patches() ; ++xpn)
    {
    patch& xp = ps.ith_patch(xpn);

	for (int x_irho = xp.min_irho() ; x_irho <= xp.max_irho() ; ++x_irho)
	{
	for (int x_isigma = xp.min_isigma() ;
	x_isigma <= xp.max_isigma() ;
	++x_isigma)
	{
	//
	// compute the main Jacobian terms for this grid point, i.e.
	//	partial Theta(this point x, Jacobian row II)
	//	---------------------------------------------
	//	partial h(other points y, Jacobian column JJ)
	//

	// Jacobian row index
	const int II = ps.gpn_of_patch_irho_isigma(xp, x_irho, x_isigma);

	// Jacobian coefficients for this point
	const fp Jacobian_coeff_rho
	   = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_d_h_1,
		       x_irho, x_isigma);
	const fp Jacobian_coeff_sigma
	   = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_d_h_2,
		       x_irho, x_isigma);
	const fp Jacobian_coeff_rho_rho
	   = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_11,
		       x_irho, x_isigma);
	const fp Jacobian_coeff_rho_sigma
	   = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_12,
		       x_irho, x_isigma);
	const fp Jacobian_coeff_sigma_sigma
	   = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_22,
		       x_irho, x_isigma);

	// partial_rho, partial_rho_rho
	      {
	    for (int m_irho = xp.molecule_min_m() ;
		 m_irho <= xp.molecule_max_m() ;
		 ++m_irho)
	    {
	    const int xm_irho = x_irho + m_irho;
	    const fp Jac_rho     = Jacobian_coeff_rho
				   * xp.partial_rho_coeff(m_irho);
	    const fp Jac_rho_rho = Jacobian_coeff_rho_rho
				   * xp.partial_rho_rho_coeff(m_irho);
	    const fp Jac_sum = Jac_rho + Jac_rho_rho;
	    if (xp.is_in_nominal_grid(xm_irho, x_isigma))
	       then {
		    const int xm_JJ
		       = Jac.II_of_patch_irho_isigma(xp,xm_irho,x_isigma);
		    Jac.sum_into_element(II, xm_JJ, Jac_sum);
		    }
	       else add_ghost_zone_Jacobian
			(ps, Jac,
			 Jac_sum,
			 xp, xp.minmax_rho_ghost_zone(m_irho < 0),
			 II, xm_irho, x_isigma);
	    }
	      }

	// partial_sigma, partial_sigma_sigma
	      {
	    for (int m_isigma = xp.molecule_min_m() ;
		 m_isigma <= xp.molecule_max_m() ;
		 ++m_isigma)
	    {
	    const int xm_isigma = x_isigma + m_isigma;
	    const fp Jac_sigma       = Jacobian_coeff_sigma
				       * xp.partial_sigma_coeff(m_isigma);
	    const fp Jac_sigma_sigma = Jacobian_coeff_sigma_sigma
				       * xp.partial_sigma_sigma_coeff(m_isigma);
	    const fp Jac_sum = Jac_sigma + Jac_sigma_sigma;
	    if (xp.is_in_nominal_grid(x_irho, xm_isigma))
	       then {
		    const int xm_JJ
		       = Jac.II_of_patch_irho_isigma(xp, x_irho, xm_isigma);
		    Jac.sum_into_element(II, xm_JJ, Jac_sum);
		    }
	       else add_ghost_zone_Jacobian
			(ps, Jac,
			 Jac_sum,
			 xp, xp.minmax_sigma_ghost_zone(m_isigma < 0),
			 II, x_irho, xm_isigma);
	    }
	      }

	// partial_rho_sigma
	      {
	    for (int m_irho = xp.molecule_min_m() ;
		 m_irho <= xp.molecule_max_m() ;
		 ++m_irho)
	    {
	    for (int m_isigma = xp.molecule_min_m() ;
		 m_isigma <= xp.molecule_max_m() ;
		 ++m_isigma)
	    {
	    const int xm_irho   = x_irho   + m_irho;
	    const int xm_isigma = x_isigma + m_isigma;
	    const fp Jac_rho_sigma
	       = Jacobian_coeff_rho_sigma
		 * xp.partial_rho_sigma_coeff(m_irho, m_isigma);
	    if (xp.is_in_nominal_grid(xm_irho, xm_isigma))
	       then {
		    const int xm_JJ
		       = Jac.II_of_patch_irho_isigma(xp, xm_irho, xm_isigma);
		    Jac.sum_into_element(II, xm_JJ, Jac_rho_sigma);
		    }
	       else {
		    const ghost_zone& xmgz
		       = xp.corner_ghost_zone_containing_point
				(m_irho < 0, m_isigma < 0,
				 xm_irho, xm_isigma);
		    add_ghost_zone_Jacobian(ps, Jac,
					    Jac_rho_sigma,
					    xp, xmgz,
					    II, xm_irho, xm_isigma);
		    }
	    }
	    }
	      }

	if (ps.N_additional_points())
		{
		const int np = ps.N_grid_points();
		Jac.set_element(II,np, 0.0);	// insert dummy value
		}

	}
	}
    }
}

//******************************************************************************

//
// This function adds the ghost-zone Jacobian dependency contributions
// for a single ghost-zone point, to a Jacobian matrix.
//
// Arguments:
// ps = The patch system.
// Jac = (out) The Jacobian matrix.
// mol = The molecule coefficient.
// xp = The patch containing the center point of the molecule.
// xmgz = If the x+m point is in a ghost zone, this must be that ghost zone.
//	  If the x+m point is not in a ghost zone, this argument is ignored.
// x_II = The Jacobian row of the x point.
// xm_(irho,isigma) = The coordinates (in xp) of the x+m point of the molecule.
//
void Horizon::add_ghost_zone_Jacobian(const patch_system& ps,
			     Jacobian& Jac,
			     fp mol,
			     const patch& xp, const ghost_zone& xmgz,
			     int x_II,
			     int xm_irho, int xm_isigma)
{
const patch_edge& xme = xmgz.my_edge();
const int xm_iperp = xme.iperp_of_irho_isigma(xm_irho, xm_isigma);
const int xm_ipar  = xme. ipar_of_irho_isigma(xm_irho, xm_isigma);

// FIXME: this won't change from one call to another
//        ==> it would be more efficient to reuse the same buffer
//            across multiple calls on this function
int global_min_ym, global_max_ym;
ps.synchronize_Jacobian_global_minmax_ym(global_min_ym, global_max_ym);
jtutil::array1d<fp> Jacobian_buffer(global_min_ym, global_max_ym);

// on what other points y does this molecule point xm depend
// via the patch_system::synchronize() operation?
int y_iperp;
int y_posn, min_ym, max_ym;
const patch_edge& ye = ps.synchronize_Jacobian(xmgz,
					       xm_iperp, xm_ipar,
					       y_iperp,
					       y_posn, min_ym, max_ym,
					       Jacobian_buffer);
patch& yp = ye.my_patch();

// add the Jacobian contributions from the ym points
	for (int ym = min_ym ; ym <= max_ym ; ++ym)
	{
	const int y_ipar = y_posn + ym;
	const int y_irho   = ye.  irho_of_iperp_ipar(y_iperp,y_ipar);
	const int y_isigma = ye.isigma_of_iperp_ipar(y_iperp,y_ipar);
	const int y_JJ = Jac.II_of_patch_irho_isigma(yp, y_irho, y_isigma);
	Jac.sum_into_element(x_II, y_JJ, mol*Jacobian_buffer(ym));
	}
}

//******************************************************************************

//
// If ps_ptr != NULL and Jac_ptr != NULL, this function sums the d/dr
// terms into the Jacobian matrix of the expansion Theta(h), computing
// those terms by finite differencing.
//
// If ps_ptr == NULL and Jac_ptr == NULL, this function does a dummy
// computation, in which only any expansion() (and hence geometry
// interpolator) calls are done, these with the number of interpolation
// points set to 0 and all the output array pointers set to NULL.
//
// It's illegal for one but not both of ps_ptr and Jac_ptr to be NULL.
//
// The basic algorithm is that
//	Jac += diag[ (Theta(h+epsilon/2) - Theta(h-epsilon/2)) / epsilon ]
//
// Inputs (angular gridfns, on ghosted grid):
//	h			# shape of trial surface
//	Theta			# Theta(h) assumed to already be computed
//				# (saved and restored, but not used)
//
// Outputs:
//	Jac += d/dr terms
//
// Results:
// This function returns a status code indicating whether the computation
// succeeded or failed, and if the latter, what caused the failure.
//

enum expansion_status
  Horizon::expansion_Jacobian_dr_FD
	(patch_system* ps_ptr, Jacobian* Jac_ptr,
         const struct what_to_compute& compute_info,
	 const struct cactus_grid_info& cgi,
	 const struct geometry_info& gi,
	 const struct Jacobian_info& Jacobian_info,
	 const struct error_info& error_info, bool initial_flag,
	 bool print_msg_flag)
{
const bool active_flag = (ps_ptr != NULL) && (Jac_ptr != NULL);
if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   horizon Jacobian: %sd/dr terms (finite diff)",
		   active_flag ? "" : "dummy ");

const fp epsilon = Jacobian_info.perturbation_amplitude;

what_to_compute this_compute_info (compute_info);
this_compute_info.surface_modification = modification_none;
this_compute_info.surface_selection = selection_definition;
this_compute_info.desired_value = 0.0;

fp additional_save_Theta;

// compute Theta(h-epsilon/2)
if (active_flag)
   then {
	ps_ptr->gridfn_copy(gfns::gfn__Theta, gfns::gfn__save_Theta);
	ps_ptr->gridfn_copy(gfns::gfn__mean_curvature, gfns::gfn__save_mean_curvature);
	if (ps_ptr->N_additional_points())
	   then {
		const int np = ps_ptr->N_grid_points();
		additional_save_Theta = ps_ptr->gridfn_data(gfns::gfn__Theta)[np];
		}
	ps_ptr->add_to_ghosted_gridfn(-epsilon/2, gfns::gfn__h);
	}
const
  enum expansion_status status = expansion(ps_ptr,
                                           this_compute_info,
					   cgi, gi,
					   error_info, initial_flag);
if (status != expansion_success)
   then {
        expansion(NULL,
                  this_compute_info,
                  cgi, gi,
                  error_info, false);
        return status;					// *** ERROR RETURN ***
  }

// compute Theta(h+epsilon/2)
if (active_flag)
   then {
	ps_ptr->gridfn_copy(gfns::gfn__Theta, gfns::gfn__old_Theta);
	ps_ptr->add_to_ghosted_gridfn(epsilon, gfns::gfn__h);
	}
const
  enum expansion_status status2 = expansion(ps_ptr,
                                            this_compute_info,
					    cgi, gi,
					    error_info, initial_flag);
if (status2 != expansion_success)
   then return status2;					// *** ERROR RETURN ***

if (active_flag)
   then {
	    for (int pn = 0 ; pn < ps_ptr->N_patches() ; ++pn)
	    {
	    patch& p = ps_ptr->ith_patch(pn);
		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		const int II = ps_ptr->gpn_of_patch_irho_isigma(p, irho,isigma);
		const fp old_Theta = p.gridfn(gfns::gfn__old_Theta,
					      irho,isigma);
		const fp new_Theta = p.gridfn(gfns::gfn__Theta,
					      irho,isigma);
		const fp d_dr_term = (new_Theta - old_Theta) / epsilon;
		Jac_ptr->sum_into_element(II,II, d_dr_term);
		}
		}
	    }

	// restore h and Theta
	ps_ptr->add_to_ghosted_gridfn(-epsilon/2, gfns::gfn__h);
	ps_ptr->gridfn_copy(gfns::gfn__save_Theta, gfns::gfn__Theta);
	ps_ptr->gridfn_copy(gfns::gfn__save_mean_curvature, gfns::gfn__mean_curvature);
	if (ps_ptr->N_additional_points())
	   then {
		const int np = ps_ptr->N_grid_points();
		ps_ptr->gridfn_data(gfns::gfn__Theta)[np] = additional_save_Theta;
		}
	}

return expansion_success;				// *** NORMAL RETURN ***
}

//******************************************************************************

//
// This function decodes the  Jacobian_compute_method  parameter (string) into
// an internal enum for future use.
//
enum Jacobian_compute_method
  Horizon::decode_Jacobian_compute_method(const char Jacobian_compute_method_string[])
{
if	(STRING_EQUAL(Jacobian_compute_method_string,
		      "numerical perturbation"))
   then return Jacobian__numerical_perturbation;
else if (STRING_EQUAL(Jacobian_compute_method_string,
		      "symbolic differentiation with finite diff d/dr"))
   then return Jacobian__symbolic_diff_with_FD_dr;
else if (STRING_EQUAL(Jacobian_compute_method_string,
		      "symbolic differentiation"))
   then return Jacobian__symbolic_diff;
else	CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   decode_Jacobian_compute_method():\n"
"        unknown Jacobian_compute_method_string=\"%s\"!",
		   Jacobian_compute_method_string);		/*NOTREACHED*/
}

//******************************************************************************

//
// This function returns (a pointer to) a C-style string describing
// an  expansion_status  value.
//
const char* Horizon::expansion_status_string(enum expansion_status status)
{
switch	(status)
	{
case expansion_success:
	return "success";
	break;
case expansion_failure__surface_nonfinite:
	return "infinity/NaN in surface shape!";
	break;
case expansion_failure__surface_too_large:
	return "surface too large";
	break;
case expansion_failure__surface_outside_grid:
	return "surface outside grid";
	break;
case expansion_failure__surface_in_excised_region:
	return "surface in excised region";
	break;
case expansion_failure__geometry_nonfinite:
	return "infinity/NaN in 3-geometry!";
	break;
case expansion_failure__gij_not_positive_definite:
	return "g_ij not positive definite!";
	break;
default:
	CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
		    "expansion_status_string(): unknown status=(int)%d!",
		    status);					/*NOTREACHED*/
	}
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function inputs a gridfn from a data file, with the file name
// chosen automagically.
//
// We assume that this gridfn is ghosted, but the ghost zones are *not*
// present in the data file.
//
//  void Horizon::input_gridfn(patch_system& ps, int unknown_gfn,
// 		  const struct IO_info& IO_info, const char base_file_name[],
// 		  int min_digits,
// 		  int hn, bool print_msg_flag, int AHF_iteration /* = 0 */)
// {
// const char* file_name = io_ASCII_file_name(IO_info, base_file_name, min_digits,
//                                            hn, AHF_iteration);

// input_gridfn__explicit_name(ps, unknown_gfn,
// 			    IO_info, base_file_name, print_msg_flag);
// }

//******************************************************************************

//
// This function inputs a gridfn from a data file, with the file name
// specified explicitly.
//
// We assume that this gridfn is ghosted, but the ghost zones are *not*
// present in the data file.
//
//  void Horizon::input_gridfn__explicit_name(patch_system& ps, int unknown_gfn,
// 				 const struct IO_info& IO_info,
// 				 const char file_name[], bool print_msg_flag)
// {
// if (print_msg_flag)
//    then {
// 	if (unknown_gfn == gfns::gfn__h)
// 	   then CCTK_VInfo(CCTK_THORNSTRING,
// 			   "   reading initial guess from \"%s\"", file_name);
// 	}

// if (IO_info.output_ASCII_files)
//    then	{
// 	ps.read_ghosted_gridfn(unknown_gfn,
// 			       file_name,
// 			       false);		// no ghost zones in data file
// 	}

// else if (IO_info.output_HDF5_files)
//    then {
// 	CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
// "input_gridfn__explicit_name(): reading from HDF5 data files not implemented yet!");
// 	}

// else
// 	{
// 	CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
// "\n"
// "   input_gridfn__explicit_name():\n"
// "        (this should never happen!)");				/*NOTREACHED*/
// 	}
// }

//******************************************************************************

//
// This function encapsulates our file-naming conventions for angular-gridfn
// output files (those used for h, H, and other angular grid functions).
//
// Arguments:
// base_file_name[] = from the parameter file
// hn = the horizon number
// AHF_iteration = the apparent horizon finder's internal iteration
//		   number (>= 1) if this is an intermediate iterate,
//		   or the default (0) if this is a final computed
//		   horizon position
//
// Results:
// This function returns (a pointer to) the file name.  The returned
// result points into a private static buffer; the usual caveats apply.

 const char* Horizon::io_HDF5_file_name(const struct IO_info& IO_info,
                              const char base_file_name[], int min_digits,
                              int hn, int AHF_iteration /* = 0 */)
{
static char file_name_buffer[IO_info::file_name_buffer_size];

const char* file_name_extension
	= IO_info.HDF5_file_name_extension;

if (AHF_iteration == 0)
   then snprintf(file_name_buffer, IO_info::file_name_buffer_size,
		 "%s/%s.ah%d.%s",
		 IO_info.h_directory, base_file_name,
		 hn,
		 file_name_extension);
   else snprintf(file_name_buffer, IO_info::file_name_buffer_size,
		 "%s/%s.ah%d.it%d.%s",
		 IO_info.h_directory, base_file_name,
		 hn, AHF_iteration,
		 file_name_extension);

return file_name_buffer;
}

const char* Horizon::io_ASCII_file_name(const struct IO_info& IO_info,
                               const char base_file_name[], int min_digits,
                               int hn, int AHF_iteration /* = 0 */)
{
static char file_name_buffer[IO_info::file_name_buffer_size];

const char* file_name_extension
	= IO_info.ASCII_gnuplot_file_name_extension;

if (AHF_iteration == 0)
   then snprintf(file_name_buffer, IO_info::file_name_buffer_size,
		 "%s/%s.t%0*d.ah%d.%s",
		 IO_info.h_directory, base_file_name,
		 min_digits, IO_info.time_iteration, hn,
		 file_name_extension);
   else snprintf(file_name_buffer, IO_info::file_name_buffer_size,
		 "%s/%s.t%0*d.ah%d.it%d.%s",
		 IO_info.h_directory, base_file_name,
		 min_digits, IO_info.time_iteration, hn, AHF_iteration,
		 file_name_extension);

return file_name_buffer;
}


//******************************************************************************

//
// This function outputs a gridfn from a data file.
//
// If the gridfn is h, then we also write out the xyz positions of the
// horizon surface points.
//
// FIXME: if the gridfn is not h, we assume that it's nominal-grid.
//
void Horizon::output_gridfn(patch_system& ps, int unknown_gfn,
                   const char gfn_name[],
		   const struct IO_info& IO_info, const char base_file_name[],
                   int min_digits,
		   int hn, bool print_msg_flag, int AHF_iteration /* = 0 */)
{
  
if (IO_info.output_ASCII_files)
   then	{
	const char* file_name
		= io_ASCII_file_name(IO_info, base_file_name, min_digits,
				     hn, AHF_iteration);
	// create the output directory (if it doesn't already exist)
	create_h_directory(IO_info);

	switch	(unknown_gfn)
		{
	case gfns::gfn__h:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing h to \"%s\"", file_name);
		ps.print_ghosted_gridfn_with_xyz
		      (unknown_gfn,
		       true, gfns::gfn__h,
		       file_name,
		       IO_info.output_ghost_zones_for_h); // should we include
							  // ghost zones?
		break;
	case gfns::gfn__Theta:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing Theta to \"%s\"", file_name);
		ps.print_gridfn(unknown_gfn, file_name);
		break;
	case gfns::gfn__mean_curvature:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing mean curvature to \"%s\"", file_name);
		ps.print_gridfn(unknown_gfn, file_name);
		break;
	case gfns::gfn__Delta_h:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing Delta_h to \"%s\"", file_name);
		ps.print_gridfn(unknown_gfn, file_name);
		break;
	default:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing gfn=%d to \"%s\"",
				   unknown_gfn, file_name);
		ps.print_gridfn(unknown_gfn, file_name);
		break;
		}
	}

  
if (IO_info.output_HDF5_files)
   then	{
	const char* file_name
		= io_HDF5_file_name(IO_info, base_file_name, min_digits,
				     hn, AHF_iteration);
	// create the output directory (if it doesn't already exist)
	create_h_directory(IO_info);
	switch	(unknown_gfn)
		{
	case gfns::gfn__h:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing h to \"%s\"", file_name);
		ps.output_ghosted_gridfn_with_xyz
		      (unknown_gfn, gfn_name,
		       true, gfns::gfn__h,
		       file_name,
		       IO_info.output_ghost_zones_for_h); // should we include
							  // ghost zones?
		break;
	case gfns::gfn__Theta:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing Theta to \"%s\"", file_name);
		ps.output_gridfn(unknown_gfn, gfn_name,
                                 file_name);
		break;
	case gfns::gfn__Delta_h:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing Delta_h to \"%s\"", file_name);
		ps.output_gridfn(unknown_gfn, gfn_name,
                                 file_name);
		break;
	default:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing gfn=%d to \"%s\"",
				   unknown_gfn, file_name);
		ps.output_gridfn(unknown_gfn, gfn_name,
                                 file_name);
		break;
		}
	}
}

void Horizon::do_evaluate_expansions(int my_proc, int N_horizons,
                                     horizon_sequence& hs,
                                     struct AH_data* const AH_data_array[],
                                     const struct cactus_grid_info& cgi,
                                     const struct geometry_info& gi,
                                     const struct IO_info& IO_info,
                                     const struct error_info& error_info,
                                     const struct verbose_info& verbose_info,
                                     int timer_handle)
{
const bool active_flag = (my_proc == 0);

if (active_flag)
   then {
	assert( hs.N_horizons() == N_horizons );
	assert( hs.my_N_horizons() == N_horizons );

		for (int hn = hs.init_hn() ;
		     hs.is_genuine() ;
		     hn = hs.next_hn())
		{
		assert( AH_data_array[hn] != NULL );
		struct AH_data& AH_data = *AH_data_array[hn];
		patch_system& ps = *AH_data.ps_ptr;

		jtutil::norm<fp> Theta_norms;
		const bool Theta_ok = expansion(&ps,
                                                AH_data.compute_info,
						cgi, gi,
						error_info, true,// initial eval
						false,	// no Jacobian coeffs
						true,	// yes, print msgs
						&Theta_norms);

		if (IO_info.output_h)
		   then output_gridfn(ps, gfns::gfn__h,
                                      "h",
				      IO_info, IO_info.h_base_file_name,
                                      IO_info.h_min_digits,
				      hn, verbose_info.print_algorithm_details);

		if (Theta_ok)
		   then {
			CCTK_VInfo(CCTK_THORNSTRING,
			   "   Theta(h) rms-norm %.2e, infinity-norm %.2e",
			   Theta_norms.rms_norm(), Theta_norms.infinity_norm());
			if (IO_info.output_Theta)
			   then output_gridfn(ps, gfns::gfn__Theta,
                                              "Theta",
					      IO_info, IO_info
						       .Theta_base_file_name,
                                              IO_info.h_min_digits,
					      hn, verbose_info
						  .print_algorithm_details);
			if (IO_info.output_mean_curvature)
			   then output_gridfn(ps, gfns::gfn__mean_curvature,
                                              "mean_curvature",
					      IO_info, IO_info
						       .mean_curvature_base_file_name,
                                              IO_info.h_min_digits,
					      hn, verbose_info
						  .print_algorithm_details);
			}
		}
	}
   else {
                struct what_to_compute new_compute_info;
		for (int i = 0 ; i < N_horizons ; ++i)
		{
                expansion(NULL, new_compute_info,
			  cgi, gi,
			  error_info, true);	// initial evaluation
		}
	}
}

// void Horizon::do_test_expansion_Jacobians(int my_proc, int N_horizons,
// 				 struct AH_data* const AH_data_array[],
// 				 const struct cactus_grid_info& cgi,
// 				 const struct geometry_info& gi,
// 				       struct Jacobian_info& Jac_info,
// 				 bool test_all_Jacobian_compute_methods,
// 				 const struct IO_info& IO_info,
// 				 const struct error_info& error_info,
// 				 const struct verbose_info& verbose_info,
// 				 int timer_handle)
// {
// const bool active_flag = (my_proc == 0);
// assert(N_horizons >= 1);

// const bool print_msg_flag = true;
// const int hn = 1;

// struct AH_data* const AH_data_ptr = active_flag ? AH_data_array[hn]   : NULL;
// patch_system*   const      ps_ptr = active_flag ? AH_data_ptr->ps_ptr : NULL;
// struct what_to_compute dummy_compute_info;
// struct what_to_compute & compute_info =
// 	active_flag
// 	? AH_data_ptr->compute_info
// 	: dummy_compute_info;

// //
// // numerical-perturbation Jacobian
// //
// Jacobian* Jac_NP_ptr = active_flag ? AH_data_ptr->Jac_ptr : NULL;
// expansion(ps_ptr, compute_info,
// 	  cgi, gi,
// 	  error_info, true);		// initial evaluation
// Jac_info.Jacobian_compute_method = Jacobian__numerical_perturbation;
// expansion_Jacobian(ps_ptr, Jac_NP_ptr,
//                    compute_info,
// 		   cgi, gi, Jac_info,
// 		   error_info, true,	// initial evaluation
// 		   print_msg_flag);

// Jacobian* Jac_SD_FDdr_ptr = NULL;
// if (test_all_Jacobian_compute_methods)
//    then {
// 	// symbolic differentiation with finite diff d/dr
// 	Jac_SD_FDdr_ptr = active_flag
// 			  ? new_Jacobian(Jac_info.Jacobian_store_solve_method,
// 					 *ps_ptr,
// 					 verbose_info.print_algorithm_details)
// 			  : NULL;
// 	expansion(ps_ptr, compute_info,
// 		  cgi, gi,
// 		  error_info, true,	// initial evaluation
// 		  true);		// compute SD Jacobian coeffs
// 	Jac_info.Jacobian_compute_method = Jacobian__symbolic_diff_with_FD_dr;
// 	expansion_Jacobian(ps_ptr, Jac_SD_FDdr_ptr,
//                            compute_info,
// 			   cgi, gi, Jac_info,
// 			   error_info, true,	// initial evaluation
// 			   print_msg_flag);
// 	}

// if (active_flag)
//    then output_Jacobians(*ps_ptr,
// 			 Jac_NP_ptr, Jac_SD_FDdr_ptr,
// 			 IO_info, IO_info.Jacobian_base_file_name,
// 			 IO_info.h_min_digits,
// 			 hn, print_msg_flag);
// }

void Horizon::AHFinderDirect_find_horizons(int cctk_iteration, double cctk_time)
{

// determine whether a horizon should be found at this iteration
bool find_any = false;
for (int hn = 1; hn <= state.my_hs->N_horizons(); ++ hn)
{
  // only try to find horizons every  find_every  time steps
  const int my_find_after = find_after_individual[hn];
  const int my_dont_find_after = dont_find_after_individual[hn];
  const fp my_find_after_time = find_after_individual_time[hn];
  const fp my_dont_find_after_time = dont_find_after_individual_time[hn];
  const int my_find_every = (find_every_individual[hn] >= 0
                             ? find_every_individual[hn]
                             : find_every);
  const bool find_this = cctk_iteration >= my_find_after
                         && (my_dont_find_after < 0
                             ? true
                             : cctk_iteration <= my_dont_find_after)
                         && cctk_time >= my_find_after_time
                         && (my_dont_find_after_time <= my_find_after_time
                             ? true
                             : cctk_time <= my_dont_find_after_time)
                         && my_find_every > 0
                         && cctk_iteration % my_find_every == 0
                         && ! disable_horizon[hn];
  struct AH_data& AH_data = *state.AH_data_array[hn];
  AH_data.search_flag = find_this;
  find_any = find_any || find_this;
}
if (! find_any) return;

const int my_proc = state.my_proc;
horizon_sequence& hs = *state.my_hs;
const bool active_flag = hs.has_genuine_horizons();
const bool broadcast_horizon_shape = true;

      struct cactus_grid_info&          cgi = state.cgi;
const struct    geometry_info&           gi = state.gi;
      struct    Jacobian_info&     Jac_info = state.Jac_info;
      struct          IO_info&      IO_info = state.IO_info;
const struct       error_info&   error_info = state.error_info;
const struct     verbose_info& verbose_info = state.verbose_info;

// what are the semantics of the Cactus gxx variables? (these may
// change from one call to another, so we have to re-check each time)
// setting whether use conformal metric here 
   cgi.use_Cactus_conformal_metric = false;

// update parameters
IO_info.output_ASCII_files = (output_ASCII_files != 0);
IO_info.output_HDF5_files = (output_HDF5_files != 0);
IO_info.output_initial_guess = (output_initial_guess != 0);
IO_info.output_h_every     = output_h_every;
IO_info.output_Theta_every = output_Theta_every;
IO_info.output_mean_curvature_every = output_mean_curvature_every;
IO_info.output_h     = false;	// dummy value
IO_info.output_Theta = false;	// dummy value
IO_info.output_mean_curvature = false;	// dummy value

IO_info.output_BH_diagnostics              = (output_BH_diagnostics != 0);
IO_info.BH_diagnostics_directory
	= (strlen(BH_diagnostics_directory) == 0)
	  ? /* IO:: */ cur_directory
	  : BH_diagnostics_directory;
IO_info.BH_diagnostics_base_file_name      = BH_diagnostics_base_file_name;
IO_info.BH_diagnostics_file_name_extension = BH_diagnostics_file_name_extension;

IO_info.output_ghost_zones_for_h  = (output_ghost_zones_for_h != 0);
IO_info.ASCII_gnuplot_file_name_extension = ASCII_gnuplot_file_name_extension;
IO_info.HDF5_file_name_extension          = HDF5_file_name_extension;
IO_info.h_directory
	= (strlen(h_directory) == 0)
	  ? /* IO:: */ cur_directory
	  : h_directory;
IO_info.h_base_file_name         = h_base_file_name;
IO_info.Theta_base_file_name     = Theta_base_file_name;
IO_info.mean_curvature_base_file_name     = mean_curvature_base_file_name;
IO_info.Delta_h_base_file_name   = Delta_h_base_file_name;
IO_info.h_min_digits             = h_min_digits;
IO_info.Jacobian_base_file_name  = Jacobian_base_file_name;
IO_info.output_OpenDX_control_files  = (output_OpenDX_control_files != 0);
IO_info.OpenDX_control_file_name_extension = OpenDX_control_file_name_extension;
IO_info.time_iteration = 0;
IO_info.time           = 0.0;

// get the Cactus time step and decide if we want to output h and/or Theta now
IO_info.time_iteration = cctk_iteration;
IO_info.time           = cctk_time;
IO_info.output_h
   = (IO_info.output_h_every > 0)
     && ((IO_info.time_iteration % IO_info.output_h_every) == 0);
IO_info.output_Theta
   = (IO_info.output_Theta_every > 0)
     && ((IO_info.time_iteration % IO_info.output_Theta_every) == 0);
IO_info.output_mean_curvature
   = (IO_info.output_mean_curvature_every > 0)
     && ((IO_info.time_iteration % IO_info.output_mean_curvature_every) == 0);

// set initial guess for any (genuine) horizons that need it,
// i.e. for any (genuine) horizons where we didn't find the horizon previously
	for (int hn = hs.init_hn() ; hs.is_genuine() ; hn = hs.next_hn())
	{
	assert( state.AH_data_array[hn] != NULL );
	struct AH_data& AH_data = *state.AH_data_array[hn];
        if (verbose_info.print_algorithm_details) {
          printf ("AHF find_horizons[%d] initial_find_flag=%d\n", hn, (int) AH_data.initial_find_flag);
          printf ("AHF find_horizons[%d] really_initial_find_flag=%d\n", hn, (int) AH_data.really_initial_find_flag);
          printf ("AHF find_horizons[%d] search_flag=%d\n", hn, (int) AH_data.search_flag);
          printf ("AHF find_horizons[%d] found_flag=%d\n", hn, (int) AH_data.found_flag);
        }
	if (AH_data.found_flag)
           then {
                AH_data.initial_find_flag = false;
                AH_data.really_initial_find_flag = false;
                }
	   else {
                if (AH_data.really_initial_find_flag
                    || AH_data.initial_guess_info.reset_horizon_after_not_finding)
                   then {
		        patch_system& ps = *AH_data.ps_ptr;
                        if (verbose_info.print_algorithm_details) {
                          printf ("AHF find_horizons[%d] setup_initial_guess\n", hn);
                        }
                        if (track_origin_from_grid_scalar[hn] && state.method == method__find_horizons) {
                           track_origin(ps, &AH_data, hn, verbose_info.print_algorithm_details);
                           set_initial_guess_parameters(AH_data, hn, 
                                                        ps.origin_x(), ps.origin_y(), ps.origin_z());
                        }
        		setup_initial_guess(ps,
        		          	    AH_data.initial_guess_info,
        				    IO_info,
        				    hn, N_horizons, verbose_info);
        		if (active_flag && IO_info.output_initial_guess)
        		   then output_gridfn(ps, gfns::gfn__h,
                                              "h",
        				      IO_info, IO_info.h_base_file_name,
                                              IO_info.h_min_digits,
        				      hn, verbose_info
        				      .print_algorithm_highlights);
        		AH_data.initial_find_flag = true;
        		}
                }
	}

//
// now the main horizon finding (or other computation)
//
switch	(state.method)
	{
case method__evaluate_expansions:
	do_evaluate_expansions(my_proc, N_horizons,
			       *state.my_hs, state.AH_data_array,
			       cgi, gi, IO_info,
			       error_info, verbose_info,
			       state.timer_handle);
	break;

case method__test_expansion_Jacobians:
	// do_test_expansion_Jacobians(my_proc, N_horizons,
	// 			    state.AH_data_array,
	// 			    cgi, gi, Jac_info,
	// 			    (test_all_Jacobian_compute_methods != 0),
	// 			    IO_info, error_info, verbose_info,
	// 			    state.timer_handle);
     error_exit(ERROR_EXIT,
                "not support test expansion jacobian method");					/*NOTREACHED*/

	break;

case method__find_horizons:
	  {
	Newton(
	       state.N_procs, state.N_active_procs, my_proc,
	       *state.my_hs, state.AH_data_array,
	       cgi, gi, Jac_info, state.solver_info,
	       IO_info, state.BH_diagnostics_info, broadcast_horizon_shape,
	       error_info, verbose_info,
	       state.isb);
	break;
	  }

default:
	CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   find_horizons(): unknown method=(int)%d!\n"
"                    (this should never happen!)"
		   ,
		   int(state.method));				/*NOTREACHED*/
	}

if (state.timer_handle >= 0)
   then {
	// CCTK_VInfo(CCTK_THORNSTRING,
	// 	   "timer stats for computation:");
	}
}

void Horizon::setup_initial_guess(patch_system& ps,
                                  const struct initial_guess_info& igi,
                                  const struct IO_info& IO_info,
                                  int hn, int N_horizons,
                                  const struct verbose_info& verbose_info)
{
if (verbose_info.print_algorithm_highlights)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "setting initial guess for horizon %d/%d",
		   hn, N_horizons);

switch	(igi.method)
	{
// case initial_guess__read_from_named_file:
// 	input_gridfn__explicit_name(ps, gfns::gfn__h,
// 				    IO_info,
// 				    igi.read_from_named_file_info.file_name,
// 				    verbose_info.print_algorithm_highlights);
// 	break;

// case initial_guess__read_from_h_file:
// 	input_gridfn(ps, gfns::gfn__h,
// 		     IO_info, IO_info.h_base_file_name, IO_info.h_min_digits,
// 		     hn, verbose_info.print_algorithm_highlights);
// 	break;

// case initial_guess__Kerr_Kerr:
// 	setup_Kerr_horizon(ps,
// 			   igi.Kerr_Kerr_info.x_posn,
// 			   igi.Kerr_Kerr_info.y_posn,
// 			   igi.Kerr_Kerr_info.z_posn,
// 			   igi.Kerr_Kerr_info.mass,
// 			   igi.Kerr_Kerr_info.spin,
// 			   false,		// Kerr coordinates
// 			   verbose_info);
// 	break;

// case initial_guess__Kerr_KerrSchild:
// 	setup_Kerr_horizon(ps,
// 			   igi.Kerr_KerrSchild_info.x_posn,
// 			   igi.Kerr_KerrSchild_info.y_posn,
// 			   igi.Kerr_KerrSchild_info.z_posn,
// 			   igi.Kerr_KerrSchild_info.mass,
// 			   igi.Kerr_KerrSchild_info.spin,
// 			   true,		// Kerr-Schild coordinates
// 			   verbose_info);
// 	break;

case initial_guess__coord_sphere:
	setup_coord_ellipsoid(ps,
			      igi.coord_sphere_info.x_center,
			      igi.coord_sphere_info.y_center,
			      igi.coord_sphere_info.z_center,
			      igi.coord_sphere_info.radius,
			      igi.coord_sphere_info.radius,
			      igi.coord_sphere_info.radius,
			      verbose_info.print_algorithm_highlights);
	break;

case initial_guess__coord_ellipsoid:
	setup_coord_ellipsoid(ps,
			      igi.coord_ellipsoid_info.x_center,
			      igi.coord_ellipsoid_info.y_center,
			      igi.coord_ellipsoid_info.z_center,
			      igi.coord_ellipsoid_info.x_radius,
			      igi.coord_ellipsoid_info.y_radius,
			      igi.coord_ellipsoid_info.z_radius,
			      verbose_info.print_algorithm_highlights);
	break;

default:
	CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "unknown initial guess method=(int)%d!",
		   int(igi.method));				/*NOTREACHED*/
	}
}

void Horizon::setup_coord_ellipsoid(patch_system& ps,
			   fp x_center, fp y_center, fp z_center,
			   fp x_radius, fp y_radius, fp z_radius,
			   bool print_msg_flag)
{
if (print_msg_flag)
   then {
	CCTK_VInfo(CCTK_THORNSTRING,
		   "   setting ellipsoid: center=(%g,%g,%g)",
		   double(x_center), double(y_center), double(z_center));
	CCTK_VInfo(CCTK_THORNSTRING,
		   "                      radius=(%g,%g,%g)",
		   double(x_radius), double(y_radius), double(z_radius));
	}

	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	patch& p = ps.ith_patch(pn);

		for (int irho = p.min_irho() ; irho <= p.max_irho() ; ++irho)
		{
		for (int isigma = p.min_isigma() ;
		     isigma <= p.max_isigma() ;
		     ++isigma)
		{
		const fp rho = p.rho_of_irho(irho);
		const fp sigma = p.sigma_of_isigma(isigma);
		fp xcos, ycos, zcos;
		p.xyzcos_of_rho_sigma(rho,sigma, xcos,ycos,zcos);

		// set up variables used by Maple-generated code
		const fp AU = x_center - ps.origin_x();
		const fp BV = y_center - ps.origin_y();
		const fp CW = z_center - ps.origin_z();
		const fp a = x_radius;
		const fp b = y_radius;
		const fp c = z_radius;

		// compute the solutions r_plus and r_minus
		fp r_plus, r_minus;
		#include "driver/ellipsoid.c"

		// exactly one of the solutions (call it r) should be positive
		fp r;
		if      ((r_plus > 0.0) && (r_minus < 0.0))
		   then r = r_plus;
		else if ((r_plus < 0.0) && (r_minus > 0.0))
		   then r = r_minus;
		else    CCTK_VWarn(FATAL_ERROR,
				   __LINE__, __FILE__, CCTK_THORNSTRING,
				   "\n"
"   setup_coord_ellipsoid():\n"
"        expected exactly one r>0 solution to quadratic, got 0 or 2!\n"
"        %s patch (irho,isigma)=(%d,%d) ==> (rho,sigma)=(%g,%g)\n"
"        direction cosines (xcos,ycos,zcos)=(%g,%g,%g)\n"
"        r_plus=%g r_minus=%g\n"
"        ==> this probably means the initial guess surface doesn't contain\n"
"            the local origin point, or more generally that the initial\n"
"            guess surface isn't a Strahlkoerper (\"star-shaped region\")\n"
"            with respect to the local origin point\n"
				   ,
				   p.name(), irho, isigma,
				   double(rho), double(sigma),
				   double(xcos), double(ycos), double(zcos),
				   double(r_plus), double(r_minus));
		   						/*NOTREACHED*/

		// r = horizon radius at this grid point
		p.ghosted_gridfn(gfns::gfn__h, irho,isigma) = r;
		}
		}
	}
}

void Horizon::output_OpenDX_control_file(const patch_system& ps,
				const struct IO_info& IO_info,
				int hn)
{
if (! IO_info.output_OpenDX_control_files) return;
static char file_name_buffer[IO_info::file_name_buffer_size];
snprintf(file_name_buffer, IO_info::file_name_buffer_size,
	 "%s/%s.ah%d.%s",
	 IO_info.h_directory, IO_info.h_base_file_name,
         hn, IO_info.OpenDX_control_file_name_extension);

FILE *fileptr = fopen(file_name_buffer, "w");
if (fileptr == NULL)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
		"output_OpenDX_control_file(): can't open output file \"%s\"!",
		   file_name_buffer);				/*NOTREACHED*/

fprintf(fileptr, "# list the size of each patch (N_rho x N_sigma)\n");
	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	const patch& p = ps.ith_patch(pn);
	fprintf(fileptr, "object \"%s patch\" class array ", p.name());
	fprintf(fileptr, "type int rank 1 shape 2 items 1 data follows %d %d\n",
		p.effective_N_irho(IO_info.output_ghost_zones_for_h),
		p.effective_N_isigma(IO_info.output_ghost_zones_for_h));
	}
fprintf(fileptr, "\n");

fprintf(fileptr, "# collect all patch sizes into a single OpenDX group\n");
fprintf(fileptr, "# for the ImportAHFinderDirectGnuplot macro to read\n");
fprintf(fileptr, "object \"patchsizes\" class group\n");
	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	const patch& p = ps.ith_patch(pn);
	fprintf(fileptr, "member %d value \"%s patch\"\n", pn, p.name());
	}

fclose(fileptr);
}

void Horizon::setup_h_files(patch_system& ps, const struct IO_info& IO_info, int hn)
{
// create the output directory (if it doesn't already exist)
create_h_directory(IO_info);
output_OpenDX_control_file(ps, IO_info, hn);
}

void Horizon::create_h_directory(const struct IO_info& IO_info)
{
// create the output directory (if it doesn't already exist)
const int status = CCTK_CreateDirectory(IO_info.default_directory_permission,
					IO_info.h_directory);
if (status < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   setup_h_files():\n"
"        error %d trying to create output directory\n"
"        \"%s\"!"
		   ,
		   status,
		   IO_info.h_directory);			/*NOTREACHED*/
}


int Horizon::CCTK_InterpGridArrays(
  int N_dims,
  int param_table_handle,
  int local_interp_handle,
  int coord_system_handle,
  int N_interp_points,
  int interp_coords_type,
  const void *const interp_coords[],
  int N_input_arrays,
  const CCTK_INT input_array_indices[],
  int N_output_arrays,
  const CCTK_INT output_array_types[],
  void *const output_arrays[])
{

  return 1;
}

}

