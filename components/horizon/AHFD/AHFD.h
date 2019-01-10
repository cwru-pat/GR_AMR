#ifndef AHFD_H
#define AHFD_H

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>

#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/hier/LocalId.h"

#include "../../bssn/bssn.h"

#include "AHFD_macros.h"
#include "driver/BH_diagnostics.hh"
#include "jtutil/util.hh"
#include "jtutil/array.hh"
#include "jtutil/util_String.h"
#include "jtutil/cpm_map.hh"
#include "jtutil/linear_map.hh"
#include "jtutil/interpolator/InterpLocalUniform.h"

#include "patch/grid.hh"
#include "patch/fd_grid.hh"
#include "patch/patch.hh"
#include "patch/patch_edge.hh"
#include "patch/patch_interp.hh"
#include "patch/ghost_zone.hh"
#include "patch/patch_system.hh"

#include "elliptic/Jacobian.hh"

#include "gr/gfns.hh"
#include "gr/gr.hh"

#include "driver/horizon_sequence.hh"
#include "AHFD_types.h"

using namespace SAMRAI;
using namespace cosmo;

#define AHFD_VEC_INIT(field, value)             \
  for(int i = 0; i <= 100 ; i++) field[i] = value

#define AHFD_DEFINE_TEMP_MI(I, J)          \
  double tempmi##I##J = 0;

#define AHFD_DEFINE_TEMP_MJ(I, J)          \
  double tempmj##I##J = 0;


#define AHFD_DEFINE_TEMP_KI(I, J)          \
  double tempKi##I##J = 0;

#define AHFD_DEFINE_TEMP_KJ(I, J)          \
  double tempKj##I##J = 0;


#define AHFD_DEFINE_TEMP_DMI(K, I, J)            \
  double tempdmi##K##I##J = 0;

#define AHFD_DEFINE_TEMP_DMJ(K, I, J)            \
  double tempdmj##K##I##J = 0;

#define AHFD_INTERPOLATE_DM_1(K, I, J) \
  tempdmi##K##I##J += 1.0 / PW2(bd.chi) * \
    ( -2.0 * bd.d##K##chi * bd.gamma##I##J / bd.chi + bd.d##K##g##I##J ) *     \
    ( (dx[0] - (x0 - x ) * (2.0 * (i - i0) -1.0)) / dx[0])

#define AHFD_INTERPOLATE_DM_2(K, I, J)           \
  tempdmj##K##I##J += tempdmi##K##I##J * \
    ( (dx[1] - (y0 - y ) * (2.0 * (j - j0) -1.0))/ dx[1])

#define AHFD_INTERPOLATE_DM_3(K, I, J)           \
  d##K##g##I##J += tempdmj##K##I##J * \
    ((dx[2] - (z0 - z ) * (2.0 * (k - k0) -1.0)) / dx[2])


#define AHFD_INTERPOLATE_M_1(I, J)  \
  tempmi##I##J += 1.0 / PW2(bd.chi) * bd.gamma##I##J *        \
    ( (dx[0] - (x0 - x ) * (2.0 * (i - i0) -1.0)) / dx[0])

#define AHFD_INTERPOLATE_M_2(I, J)  \
  tempmj##I##J += tempmi##I##J * \
    ( (dx[1] - (y0 - y ) * (2.0 * (j - j0) -1.0))/ dx[1])

#define AHFD_INTERPOLATE_M_3(I, J)  \
  g##I##J += tempmj##I##J * \
    ((dx[2] - (z0 - z ) * (2.0 * (k - k0) -1.0)) / dx[2])


#define AHFD_INTERPOLATE_K_1(I, J)  \
  tempKi##I##J += 1.0 / PW2(bd.chi) * (bd.A##I##J +  bd.gamma##I##J * bd.K / 3.0) * \
    ( (dx[0] - (x0 - x ) * (2.0 * (i - i0) -1.0)) / dx[0])

#define AHFD_INTERPOLATE_K_2(I, J)  \
  tempKj##I##J += tempKi##I##J * \
    ( (dx[1] - (y0 - y ) * (2.0 * (j - j0) -1.0))/ dx[1])


#define AHFD_INTERPOLATE_K_3(I, J)  \
  K##I##J += tempKj##I##J * \
    ((dx[2] - (z0 - z ) * (2.0 * (k - k0) -1.0)) / dx[2])




namespace AHFinderDirect
{

//
// (A single copy of) this struct holds all of our information about
// a single apparent horizon.
//
struct	AH_data
	{
	//
	// Any given horizon is allocated to a single processor.
	// On that processor (where we actually find the horizon)
	// we keep a "full-fledged" patch system and a Jacobian.
	// On other processors (where we only use the horizon shape
	// to (optionally) set the BH excision mask), we keep only a
	// "skeletal" patch system, and no Jacobian.
	//
	patch_system* ps_ptr;
	Jacobian* Jac_ptr;
        what_to_compute compute_info;

        bool move_origins;

        bool use_pretracking;
        int pretracking_max_iterations;

        fp pretracking_value;
        fp pretracking_minimum_value;
        fp pretracking_maximum_value;
        fp pretracking_delta;
        fp pretracking_minimum_delta;
        fp pretracking_maximum_delta;

        int depends_on;
        fp desired_value_factor;
        fp desired_value_offset;

        fp shiftout_factor;
        fp smoothing_factor;

	// are we finding this horizon "from scratch"?
	// ... true if this is the first time we've tried to find it,
	//          or if we've tried before but failed to find it the
	//          last time
	//     false if we've tried before and succeeded the last time
	//           (so we have that position as a very good initial guess)
	// ... also false if we're not finding this horizon on this processor
	bool initial_find_flag;
        // is this really the first time in this simulation that we are
        // trying to find a horizon?
        // ... true initially if this is a genuine horizon,
        //     false after the first time that a horizon has been found
	bool really_initial_find_flag;

	// used only if we're finding this horizon on this processor
	struct initial_guess_info initial_guess_info;

	bool search_flag;	// did we search for this horizon
	bool found_flag;	// did we find this horizon (successfully)
	bool h_files_written;	// have we written horizon-shape or similar
				// files for this horizon yet?

	struct BH_diagnostics BH_diagnostics;
	FILE *BH_diagnostics_fileptr;

	// interprocessor-communication buffers
	// for this horizon's BH diagnostics and (optionally) horizon shape
	struct horizon_buffers horizon_buffers;
	};


struct	state
{
  enum method method;
  struct error_info error_info;
  struct verbose_info verbose_info;
  int timer_handle;
  int N_procs;			// total number of processors
  int my_proc;			// processor number of this processor
  // (0 to N_procs-1)

  int N_horizons;			// total number of genuine horizons
					// being searched for
  int N_active_procs;		// total number of active processors
  // (the active processors are processor
  //  numbers 0 to N_active_procs-1)

  struct cactus_grid_info cgi;
  struct geometry_info gi;
  struct Jacobian_info Jac_info;
  struct solver_info solver_info;
  struct IO_info IO_info;

  struct BH_diagnostics_info BH_diagnostics_info;
  struct mask_info mask_info;

  // interprocessor-communication buffers for broadcasting
  // Newton-iteration status from active processors to all processors
  struct iteration_status_buffers isb;

  horizon_sequence *my_hs;	// --> new-allocated object describing
  //     the sequence of genuine horizons
  //     assigned to this processor

  // horizon numbers ("hn") run from 1 to N_horizons inclusive
  struct AH_data** AH_data_array;	// --> new[]-allocated array of size
					//     N_horizons+1, subscripted by hn,
					//     of --> info
};



 
  
class Horizon
{
  //
// prototypes for functions visible outside their source files
//
 public:
// setup.cc
// ... called from Cactus Scheduler

  struct state state;
  
  Horizon(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
          BSSN * bssn_in,
          const tbox::Dimension& dim_in,
          std::shared_ptr<tbox::Database> database_in,
          const char * visit_d_name,
          int w_idx_in);

  ~Horizon();
  
  void AHFinderDirect_setup();

  enum method
    decode_method(const char method_string[]);
  enum verbose_level
    decode_verbose_level(const char verbose_level_string[]);

  int allocate_horizons_to_processor(int N_procs, int my_proc,
                                     int N_horizons, bool multiproc_flag,
                                     const CCTK_INT depends_on[],
                                     horizon_sequence& my_hs,
                                     const struct verbose_info& verbose_info);

  enum patch_system::patch_system_type
    choose_patch_system_type(const char grid_domain[],
                             const char grid_bitant_plane[],
                             const char grid_quadrant_direction[],
                             const char grid_rotation_axis[],
                             fp origin_x, fp origin_y, fp origin_z);
  

// find_horizons.cc
// ... called from Cactus Scheduler
  void AHFinderDirect_find_horizons();

// announce.cc
// ... called from Cactus Scheduler
  void AHFinderDirect_announce();

// mask.cc
// ... called from Cactus Scheduler
  void AHFinderDirect_maybe_do_masks();

// aliased_functions.cc
// ... called from other thorns via the Cactus flesh function-aliasing mechanism
  CCTK_INT AHFinderDirect_local_coordinate_origin
    (CCTK_INT horizon_number,
     CCTK_REAL* origin_x_ptr, CCTK_REAL* origin_y_ptr, CCTK_REAL* origin_z_ptr);

  CCTK_INT AHFinderDirect_horizon_was_found(CCTK_INT horizon_number);

   void AHFinderDirect_find_horizons(int cctk_iteration, double cctk_time);
  
  CCTK_INT AHFinderDirect_radius_in_direction
    (CCTK_INT horizon_number,
     CCTK_INT N_points,
     const CCTK_REAL* const x, const CCTK_REAL* const y, const CCTK_REAL* const z,
     CCTK_REAL* const radius);
  // int CCTK_CreateDirectory (int mode, const char *pathname);
  
 private:
// initial_guess.cc
void setup_initial_guess(patch_system& ps,
			 const struct initial_guess_info& igi,
			 const struct IO_info& IO_info,
			 int hn, int N_horizons,
			 const struct verbose_info& verbose_info);
enum initial_guess_method
  decode_initial_guess_method(const char initial_guess_method_string[]);


void set_initial_guess_parameters(struct AH_data& AH_data, const int hn,
                                  const fp origin_x, const fp origin_y, const fp origin_z);

// Newton.cc
// returns true for success, false for failure to converge
void Newton(
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
	    struct iteration_status_buffers& isb);

// Tracks coordinate origin
void track_origin( patch_system& ps, 
                  struct AH_data* const AH_data_ptr, 
                  const int hn, const bool print_algorithm_details);

// io.cc
/* void input_gridfn(patch_system& ps, int unknown_gfn, */
/* 		  const struct IO_info& IO_info, const char base_file_name[], */
/* 		  int min_digits, */
/* 		  int hn, bool print_msg_flag, int AHF_iteration = 0); */
/* void input_gridfn__explicit_name(patch_system& ps, int unknown_gfn, */
/* 				 const struct IO_info& IO_info, */
/* 				 const char file_name[], bool print_msg_flag); */
void setup_h_files(patch_system& ps, const struct IO_info& IO_info, int hn);
 
 void output_OpenDX_control_file(const patch_system& ps,
                                 const struct IO_info& IO_info,
                                 int hn);

void output_gridfn(patch_system& ps, int unknown_gfn,
                   const char gfn_name[],
		   const struct IO_info& IO_info, const char base_file_name[],
                   int min_digits,
		   int hn, bool print_msg_flag, int AHF_iteration = 0);
void output_Jacobians(const patch_system& ps,
		      const Jacobian* Jac_NP_ptr,
		      const Jacobian* Jac_SD_FDdr_ptr,
		      const struct IO_info& IO_info, const char base_file_name[],
                      int min_digits,
		      int hn, bool print_msg_flag, int AHF_iteration = 0);

// misc-driver.cc
int Cactus_gridfn_varindex(const char gridfn_name[]);

bool broadcast_status(
		      int N_procs, int N_active_procs,
		      int my_proc, bool my_active_flag,
		      int hn, int iteration,
		      enum expansion_status expansion_status,
		      fp mean_horizon_radius, fp infinity_norm,
		      bool found_this_horizon, bool I_need_more_iterations,
		      struct iteration_status_buffers& isb);
void broadcast_horizon_data(
			    bool broadcast_flag, bool broadcast_horizon_shape,
                            struct AH_data& AH_data,
			    struct BH_diagnostics& BH_diagnostics,
			    patch_system& ps,
			    struct horizon_buffers& horizon_buffers);

void print_status(int N_active_procs,
		  const struct iteration_status_buffers& isb);
void Newton_step(patch_system& ps,
		 fp mean_horizon_radius, fp max_allowable_Delta_h_over_h,
		 const struct verbose_info& verbose_info);
	

void setup_xyz_posns(patch_system& ps, bool print_msg_flag);
enum expansion_status
  interpolate_geometry(patch_system* ps_ptr,
		       const struct cactus_grid_info& cgi,
		       const struct geometry_info& gi,
		       const struct error_info& error_info, bool initial_flag,
		       bool print_msg_flag);

void convert_conformal_to_physical(patch_system& ps,
				   bool print_msg_flag);
void Schwarzschild_EF_geometry(patch_system& ps,
			       const struct cactus_grid_info& cgi,
			       const struct geometry_info& gi,
			       bool print_msg_flag);

bool h_is_finite(patch_system& ps,
		 const struct error_info& error_info, bool initial_flag,
		 bool print_msg_flag);
bool geometry_is_finite(patch_system& ps,
			const struct error_info& error_info, bool initial_flag,
			bool print_msg_flag);

bool compute_Theta(patch_system& ps, const struct what_to_compute& comput_info,
		   bool Jacobian_flag,
                   jtutil::norm<fp>* Theta_norms_ptr,
                   jtutil::norm<fp>* expansion_Theta_norms_ptr,
                   jtutil::norm<fp>* inner_expansion_Theta_norms_ptr,
                   jtutil::norm<fp>* product_expansion_Theta_norms_ptr,
                   jtutil::norm<fp>* mean_curvature_Theta_norms_ptr,
		   const struct error_info& error_info, bool initial_flag,
		   bool print_msg_flag);

enum expansion_status
  expansion_Jacobian(patch_system* ps_ptr, Jacobian* Jac_ptr,
                     const struct what_to_compute& compute_info,
		     const struct cactus_grid_info& cgi,
		     const struct geometry_info& gi,
		     const struct Jacobian_info& Jacobian_info,
		     const struct error_info& error_info, bool initial_flag,
                              bool print_msg_flag /* = false */);

 
enum expansion_status
  expansion_Jacobian_NP
	(patch_system& ps, Jacobian& Jac,
         const struct what_to_compute& comput_info,
	 const struct cactus_grid_info& cgi,
	 const struct geometry_info& gi,
	 const struct Jacobian_info& Jacobian_info,
	 const struct error_info& error_info, bool initial_flag,
	 bool print_msg_flag);

void expansion_Jacobian_partial_SD(patch_system& ps, Jacobian& Jac,
				   const struct cactus_grid_info& cgi,
				   const struct geometry_info& gi,
				   const struct Jacobian_info& Jacobian_info,
				   bool print_msg_flag);
void add_ghost_zone_Jacobian(const patch_system& ps,
			     Jacobian& Jac,
			     fp mol,
			     const patch& xp, const ghost_zone& xmgz,
			     int x_II,
			     int xm_irho, int xm_isigma);
enum expansion_status
  expansion_Jacobian_dr_FD
	(patch_system* ps_ptr, Jacobian* Jac_ptr,
         const struct what_to_compute& compute_info,
	 const struct cactus_grid_info& cgi,
	 const struct geometry_info& gi,
	 const struct Jacobian_info& Jacobian_info,
	 const struct error_info& error_info, bool initial_flag,
	 bool print_msg_flag);

const char* expansion_status_string(enum expansion_status status);

enum Jacobian_compute_method
decode_Jacobian_compute_method(const char Jacobian_compute_method_string[]);

const char* io_HDF5_file_name(const struct IO_info& IO_info,
                              const char base_file_name[],
                              int min_digits,
                              int hn, int AHF_iteration = 0);
 
const char* io_ASCII_file_name(const struct IO_info& IO_info,
                                        const char base_file_name[], int min_digits,
                                        int hn, int AHF_iteration /* = 0 */);

void create_h_directory(const struct IO_info& IO_info);

  
CCTK_INT AHFinderDirect_horizon_centroid
    (CCTK_INT horizon_number,
     CCTK_REAL* p_centroid_x, CCTK_REAL* p_centroid_y, CCTK_REAL* p_centroid_z);

enum expansion_status 
  expansion(patch_system* ps_ptr,
            const struct what_to_compute& comput_info,
	    const struct cactus_grid_info& cgi,
	    const struct geometry_info& gi,
	    const struct error_info& error_info, bool initial_flag,
	    bool Jacobian_flag = false,
	    bool print_msg_flag = false,
	    jtutil::norm<fp>* Theta_norms_ptr = NULL,
	    jtutil::norm<fp>* expansion_Theta_norms_ptr = NULL,
	    jtutil::norm<fp>* inner_expansion_Theta_norms_ptr = NULL,
	    jtutil::norm<fp>* product_expansion_Theta_norms_ptr = NULL,
	    jtutil::norm<fp>* mean_curvature_Theta_norms_ptr = NULL);
 
void do_evaluate_expansions(int my_proc, int N_horizons,
			    horizon_sequence& hs,
			       struct AH_data* const AH_data_array[],
			    const struct cactus_grid_info& cgi,
			    const struct geometry_info& gi,
			    const struct IO_info& IO_info,
			    const struct error_info& error_info,
			    const struct verbose_info& verbose_info,
			    int timer_handle);
void do_test_expansion_Jacobians(int my_proc, int N_horizons,
				 struct AH_data* const AH_data_array[],
				 const struct cactus_grid_info& cgi,
				 const struct geometry_info& gi,
				       struct Jacobian_info& Jac_info,
				 bool test_all_Jacobian_compute_methods,
				 const struct IO_info& IO_info,
				 const struct error_info& error_info,
				 const struct verbose_info& verbose_info,
				 int timer_handle);


 void setup_Kerr_horizon(patch_system& ps,
			fp x_posn, fp y_posn, fp z_posn,
			fp m, fp a,
			bool Kerr_Schild_flag,
			const struct verbose_info& verbose_info);
void setup_coord_ellipsoid(patch_system& ps,
			   fp x_center, fp y_center, fp z_center,
			   fp x_radius, fp y_radius, fp z_radius,
			   bool print_msg_flag);

int CCTK_InterpGridArrays(
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
  void *const output_arrays[]);

//*****************************************************************************
 std::shared_ptr<hier::PatchHierarchy> hierarchy;
 BSSN * bssn;
 std::shared_ptr<tbox::Database>& AHFD_db;

 const tbox::Dimension& dim;
 int w_idx;

 double domain_lower[3], domain_upper[3];
 
//******************************************************************************
  CCTK_REAL ILUCG__error_tolerance;
  CCTK_REAL Jacobian_perturbation_amplitude;
  CCTK_REAL Theta_norm_for_convergence;
  CCTK_REAL desired_value[101];
  CCTK_REAL desired_value_factor[101];
  CCTK_REAL desired_value_offset[101];
  CCTK_REAL dont_find_after_individual_time[101];
  CCTK_REAL find_after_individual_time[101];
  CCTK_REAL geometry__Schwarzschild_EF__Delta_xyz;
  CCTK_REAL geometry__Schwarzschild_EF__epsilon;
  CCTK_REAL geometry__Schwarzschild_EF__mass;
  CCTK_REAL geometry__Schwarzschild_EF__x_posn;
  CCTK_REAL geometry__Schwarzschild_EF__y_posn;
  CCTK_REAL geometry__Schwarzschild_EF__z_posn;
  CCTK_REAL initial_guess__Kerr_KerrSchild__mass[101];
  CCTK_REAL initial_guess__Kerr_KerrSchild__spin[101];
  CCTK_REAL initial_guess__Kerr_KerrSchild__x_posn[101];
  CCTK_REAL initial_guess__Kerr_KerrSchild__y_posn[101];
  CCTK_REAL initial_guess__Kerr_KerrSchild__z_posn[101];
  CCTK_REAL initial_guess__Kerr_Kerr__mass[101];
  CCTK_REAL initial_guess__Kerr_Kerr__spin[101];
  CCTK_REAL initial_guess__Kerr_Kerr__x_posn[101];
  CCTK_REAL initial_guess__Kerr_Kerr__y_posn[101];
  CCTK_REAL initial_guess__Kerr_Kerr__z_posn[101];
  CCTK_REAL initial_guess__coord_ellipsoid__x_center[101];
  CCTK_REAL initial_guess__coord_ellipsoid__x_radius[101];
  CCTK_REAL initial_guess__coord_ellipsoid__y_center[101];
  CCTK_REAL initial_guess__coord_ellipsoid__y_radius[101];
  CCTK_REAL initial_guess__coord_ellipsoid__z_center[101];
  CCTK_REAL initial_guess__coord_ellipsoid__z_radius[101];
  CCTK_REAL initial_guess__coord_sphere__radius[101];
  CCTK_REAL initial_guess__coord_sphere__x_center[101];
  CCTK_REAL initial_guess__coord_sphere__y_center[101];
  CCTK_REAL initial_guess__coord_sphere__z_center[101];
  CCTK_REAL mask_buffer_thickness;
  CCTK_REAL mask_radius_multiplier;
  CCTK_REAL mask_radius_offset;
  CCTK_REAL max_allowable_Delta_h_over_h;
  CCTK_REAL max_allowable_Theta;
  CCTK_REAL max_allowable_horizon_radius[101];
  CCTK_REAL min_horizon_radius_points_for_mask;
  CCTK_REAL old_style_mask_buffer_value;
  CCTK_REAL old_style_mask_inside_value;
  CCTK_REAL old_style_mask_outside_value;
  CCTK_REAL origin_x[101];
  CCTK_REAL origin_y[101];
  CCTK_REAL origin_z[101];
  CCTK_REAL pretracking_delta[101];
  CCTK_REAL pretracking_maximum_delta[101];
  CCTK_REAL pretracking_maximum_value[101];
  CCTK_REAL pretracking_minimum_delta[101];
  CCTK_REAL pretracking_minimum_value[101];
  CCTK_REAL pretracking_value[101];
  CCTK_REAL shiftout_factor[101];
  CCTK_REAL smoothing_factor[101];
  const char * ASCII_gnuplot_file_name_extension;
  const char * BH_diagnostics_base_file_name;
  const char * BH_diagnostics_directory;
  const char * BH_diagnostics_file_name_extension;
  const char * Delta_h_base_file_name;
  const char * HDF5_file_name_extension;
  const char * Jacobian_base_file_name;
  const char * Jacobian_compute_method;
  const char * Jacobian_store_solve_method;
  const char * OpenDX_control_file_name_extension;
  const char * Theta_base_file_name;
  const char * coordinate_system_name;
  const char * geometry_interpolator_name;
  const char * geometry_interpolator_pars;
  const char * h_base_file_name;
  const char * h_directory;
  const char * initial_guess__read_from_named_file__file_name[101];
  const char * initial_guess_method[101];
  const char * integral_method;
  const char * interpatch_interpolator_name;
  const char * interpatch_interpolator_pars;
  const char * mean_curvature_base_file_name;
  const char * method;
  const char * new_style_mask_bitfield_name;
  const char * new_style_mask_buffer_value;
  const char * new_style_mask_gridfn_name;
  const char * new_style_mask_inside_value;
  const char * new_style_mask_outside_value;
  const char * old_style_mask_gridfn_name;
  const char * patch_system_type[101];
  const char * surface_definition[101];
  const char * surface_interpolator_name;
  const char * surface_interpolator_pars;
  const char * surface_modification[101];
  const char * surface_selection[101];
  const char * track_origin_source_x[101];
  const char * track_origin_source_y[101];
  const char * track_origin_source_z[101];
  const char * verbose_level;
  const char * which_surface_to_store_info_by_name[101];
  CCTK_INT ILUCG__limit_CG_iterations;
 public:
  CCTK_INT N_horizons;
  CCTK_INT find_every;
 private:
  CCTK_INT N_zones_per_right_angle[101];
  CCTK_INT UMFPACK__N_II_iterations;
  CCTK_INT check_that_geometry_is_finite;
  CCTK_INT check_that_h_is_finite;
  CCTK_INT debugging_output_at_each_Newton_iteration;
  CCTK_INT depends_on[101];
  CCTK_INT disable_horizon[101];
  CCTK_INT dont_find_after_individual[101];
  CCTK_INT find_after_individual[101];
  CCTK_INT find_every_individual[101];
  CCTK_INT ghost_zone_width;
  CCTK_INT h_min_digits;
  CCTK_INT hardwire_Schwarzschild_EF_geometry;
  CCTK_INT mask_is_noshrink;
  CCTK_INT max_N_zones_per_right_angle;
  CCTK_INT max_Newton_iterations__initial;
  CCTK_INT max_Newton_iterations__subsequent;
  CCTK_INT max_allowable_Theta_growth_iterations;
  CCTK_INT max_allowable_Theta_nonshrink_iterations;
  CCTK_INT move_origins;
  CCTK_INT output_ASCII_files;
  CCTK_INT output_BH_diagnostics;
  CCTK_INT output_HDF5_files;
  CCTK_INT output_OpenDX_control_files;
  CCTK_INT output_Theta_every;
  CCTK_INT output_ghost_zones_for_h;
  CCTK_INT output_h_every;
  CCTK_INT output_initial_guess;
  CCTK_INT output_mean_curvature_every;
  CCTK_INT patch_overlap_width;
  CCTK_INT predict_origin_movement;
  CCTK_INT pretracking_max_iterations[101];
  CCTK_INT print_timing_stats;
  CCTK_INT reset_horizon_after_not_finding[101];
  CCTK_INT reshape_while_moving;
  CCTK_INT set_mask_for_all_horizons;
  CCTK_INT set_mask_for_individual_horizon[101];
  CCTK_INT set_new_style_mask;
  CCTK_INT set_old_style_mask;
  CCTK_INT test_all_Jacobian_compute_methods;
  CCTK_INT track_origin_from_grid_scalar[101];
  CCTK_INT use_pretracking[101];
  CCTK_INT want_expansion_gradients;
  CCTK_INT warn_level__gij_not_positive_definite__initial;
  CCTK_INT warn_level__gij_not_positive_definite__subsequent;
  CCTK_INT warn_level__nonfinite_geometry;
  CCTK_INT warn_level__point_outside__initial;
  CCTK_INT warn_level__point_outside__subsequent;
  CCTK_INT warn_level__skipping_finite_check;
  CCTK_INT which_horizon_to_announce_centroid;
  CCTK_INT which_surface_to_store_info[101];

  char cur_directory[FILENAME_MAX];
};

}//ending namespace
#endif
