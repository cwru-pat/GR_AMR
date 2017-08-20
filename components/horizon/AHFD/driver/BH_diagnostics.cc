// BH_diagnostics.cc -- compute/print BH diagnostics
// $Header$
//
// BH_diagnostics::BH_diagnostics - initialize a  struct BH_diagnostics
//
// BH_diagnostics::copy_to_buffer - copy diagnostics to buffer
// BH_diagnostics::copy_from_buffer - copy buffer to diagnostics
//
// BH_diagnostics::compute - compute BH diagnostics after an AH has been found
// BH_diagnostics::surface_integral - integrate gridfn over the 2-sphere
//
// print - print a line or two summarizing the diagnostics
// setup_output_file - create/open output file, write header describing fields
// output - write a (long) line of all the diagnostics
// store - copy the surface and the diagnostics into the SphericalSurface
//         variables
// save - copy the surface and the diagnostics into the Cactus variables
// load - set the surface and the diagnostics from the Cactus variables
//

#include "BH_diagnostics.hh"


#include "../jtutil/util_String.h"

using namespace SAMRAI;
// all the code in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// ***** access to persistent data *****
//
//extern struct state state;

//******************************************************************************

//
// This function initializes a  struct BH_diagnostics  to all zeros.
//
BH_diagnostics::BH_diagnostics()
	: origin_x(0.0), origin_y(0.0), origin_z(0.0),
	  centroid_x(0.0), centroid_y(0.0), centroid_z(0.0),
          quadrupole_xx(0.0), quadrupole_xy(0.0), quadrupole_xz(0.0),
          quadrupole_yy(0.0), quadrupole_yz(0.0), quadrupole_zz(0.0),
	  min_radius(0.0), max_radius(0.0), mean_radius(0.0),
	  min_x(0.0), max_x(0.0),
	  min_y(0.0), max_y(0.0),
	  min_z(0.0), max_z(0.0),
	  circumference_xy(0.0), circumference_xz(0.0), circumference_yz(0.0),
	  area(1.0),
          expansion(0.0),
          inner_expansion(0.0),
          product_expansion(0.0),
          mean_curvature(0.0),
          area_gradient(0.0),
          expansion_gradient(0.0),
          inner_expansion_gradient(0.0),
          product_expansion_gradient(0.0),
          mean_curvature_gradient(0.0),
          mean_curvature_minimum(0.0),
          mean_curvature_maximum(0.0),
          mean_curvature_integral(0.0)
{ }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function copies the diagnostics to a user-supplied buffer.
//
void BH_diagnostics::copy_to_buffer(CCTK_REAL buffer[N_buffer])
	const
{
buffer[posn__origin_x] = this->origin_x;
buffer[posn__origin_y] = this->origin_y;
buffer[posn__origin_z] = this->origin_z;

buffer[posn__centroid_x] = centroid_x;
buffer[posn__centroid_y] = centroid_y;
buffer[posn__centroid_z] = centroid_z;

buffer[posn__quadrupole_xx] = quadrupole_xx;
buffer[posn__quadrupole_xy] = quadrupole_xy;
buffer[posn__quadrupole_xz] = quadrupole_xz;
buffer[posn__quadrupole_yy] = quadrupole_yy;
buffer[posn__quadrupole_xz] = quadrupole_yz;
buffer[posn__quadrupole_zz] = quadrupole_zz;

buffer[posn__min_radius]  = min_radius;
buffer[posn__max_radius]  = max_radius;
buffer[posn__mean_radius] = mean_radius;

buffer[posn__min_x] = min_x;
buffer[posn__max_x] = max_x;
buffer[posn__min_y] = min_y;
buffer[posn__max_y] = max_y;
buffer[posn__min_z] = min_z;
buffer[posn__max_z] = max_z;

buffer[posn__circumference_xy]  = circumference_xy;
buffer[posn__circumference_xz]  = circumference_xz;
buffer[posn__circumference_yz]  = circumference_yz;
buffer[posn__area]              = area;
buffer[posn__expansion]         = expansion;
buffer[posn__inner_expansion]   = inner_expansion;
buffer[posn__product_expansion] = product_expansion;
buffer[posn__mean_curvature]    = mean_curvature;

buffer[posn__area_gradient]              = area_gradient;
buffer[posn__expansion_gradient]         = expansion_gradient;
buffer[posn__inner_expansion_gradient]   = inner_expansion_gradient;
buffer[posn__product_expansion_gradient] = product_expansion_gradient;
buffer[posn__mean_curvature_gradient]    = mean_curvature_gradient;

buffer[posn__mean_curvature_minimum]    = mean_curvature_minimum;
buffer[posn__mean_curvature_maximum]    = mean_curvature_maximum;
buffer[posn__mean_curvature_integral]   = mean_curvature_integral;
}

//******************************************************************************

//
// This function copies a user-supplied buffer to the diagnostics.
//
void BH_diagnostics::copy_from_buffer(const CCTK_REAL buffer[N_buffer])
{
this->origin_x = buffer[posn__origin_x];
this->origin_y = buffer[posn__origin_y];
this->origin_z = buffer[posn__origin_z];

centroid_x = buffer[posn__centroid_x];
centroid_y = buffer[posn__centroid_y];
centroid_z = buffer[posn__centroid_z];

quadrupole_xx = buffer[posn__quadrupole_xx];
quadrupole_xy = buffer[posn__quadrupole_xy];
quadrupole_xz = buffer[posn__quadrupole_xz];
quadrupole_yy = buffer[posn__quadrupole_yy];
quadrupole_yz = buffer[posn__quadrupole_yz];
quadrupole_zz = buffer[posn__quadrupole_zz];

 min_radius = buffer[posn__min_radius];
 max_radius = buffer[posn__max_radius];
mean_radius = buffer[posn__mean_radius];

min_x = buffer[posn__min_x];
max_x = buffer[posn__max_x];
min_y = buffer[posn__min_y];
max_y = buffer[posn__max_y];
min_z = buffer[posn__min_z];
max_z = buffer[posn__max_z];

 circumference_xy = buffer[posn__circumference_xy];
 circumference_xz = buffer[posn__circumference_xz];
 circumference_yz = buffer[posn__circumference_yz];
             area = buffer[posn__area];
        expansion = buffer[posn__expansion];
  inner_expansion = buffer[posn__inner_expansion];
product_expansion = buffer[posn__product_expansion];
   mean_curvature = buffer[posn__mean_curvature];

             area_gradient = buffer[posn__area_gradient];
        expansion_gradient = buffer[posn__expansion_gradient];
  inner_expansion_gradient = buffer[posn__inner_expansion_gradient];
product_expansion_gradient = buffer[posn__product_expansion_gradient];
   mean_curvature_gradient = buffer[posn__mean_curvature_gradient];

   mean_curvature_minimum  = buffer[posn__mean_curvature_minimum];
   mean_curvature_maximum  = buffer[posn__mean_curvature_maximum];
   mean_curvature_integral = buffer[posn__mean_curvature_integral];
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// Given that an apparent horizon has been found, this function computes
// various black hole diagnostics.
//
// Inputs (gridfns)
// h		# ghosted
// one		# nominal
// global_[xyz]	# nominal
//
// Bugs:
// The computation is rather inefficient -- we make many passes over the
// angular grid, instead of doing everything in one pass.
//
void BH_diagnostics::compute
	(const patch_system& ps,
         const fp the_area,
         const fp mean_expansion,
         const fp mean_inner_expansion,
         const fp mean_product_expansion,
         const fp mean_mean_curvature,
         const fp the_area_gradient,
         const fp mean_expansion_gradient,
         const fp mean_inner_expansion_gradient,
         const fp mean_product_expansion_gradient,
         const fp mean_mean_curvature_gradient,
	 const struct BH_diagnostics_info& BH_diagnostics_info)
{
//
// min/max radius of horizon
//
jtutil::norm<fp> h_norms;
ps.ghosted_gridfn_norms(gfns::gfn__h, h_norms);
min_radius = h_norms.min_abs_value();
max_radius = h_norms.max_abs_value();


//
// xyz bounding box of horizon
//

// compute bounding box of nominal grid
// ... this is only the stored part of the horizon if there are symmetries
jtutil::norm<fp> x_norms;
ps.gridfn_norms(gfns::gfn__global_x, x_norms);
min_x = x_norms.min_value();
max_x = x_norms.max_value();

jtutil::norm<fp> y_norms;
ps.gridfn_norms(gfns::gfn__global_y, y_norms);
min_y = y_norms.min_value();
max_y = y_norms.max_value();

jtutil::norm<fp> z_norms;
ps.gridfn_norms(gfns::gfn__global_z, z_norms);
min_z = z_norms.min_value();
max_z = z_norms.max_value();

// adjust the bounding box for the symmetries
#define REFLECT(origin_, max_)	(origin_ - (max_ - origin_))
switch	(ps.type())
	{
case patch_system::patch_system__full_sphere:
	break;
case patch_system::patch_system__plus_z_hemisphere:
	min_z = REFLECT(ps.origin_z(), max_z);
	break;
case patch_system::patch_system__plus_xy_quadrant_mirrored:
case patch_system::patch_system__plus_xy_quadrant_rotating:
	min_x = REFLECT(ps.origin_x(), max_x);
	min_y = REFLECT(ps.origin_y(), max_y);
	break;
case patch_system::patch_system__plus_xz_quadrant_mirrored:
case patch_system::patch_system__plus_xz_quadrant_rotating:
	min_x = REFLECT(ps.origin_x(), max_x);
	min_z = REFLECT(ps.origin_z(), max_z);
	break;
case patch_system::patch_system__plus_xyz_octant_mirrored:
case patch_system::patch_system__plus_xyz_octant_rotating:
	min_x = REFLECT(ps.origin_x(), max_x);
	min_y = REFLECT(ps.origin_y(), max_y);
	min_z = REFLECT(ps.origin_z(), max_z);
	break;
default:
	error_exit(PANIC_EXIT,
"***** BH_diagnostics::compute(): unknown patch system type()=(int)%d!\n"
"                                 (this should never happen!)\n",
		   int(ps.type()));				/*NOTREACHED*/
	}


//
// surface integrals
//
const fp integral_one = surface_integral(ps,
					 gfns::gfn__one, true, true, true,
					 BH_diagnostics_info.integral_method);
const fp integral_h = surface_integral(ps,
				       gfns::gfn__h, true, true, true,
				       BH_diagnostics_info.integral_method);
const fp integral_x = surface_integral(ps,
				       gfns::gfn__global_x, true, true, false,
				       BH_diagnostics_info.integral_method);
const fp integral_y = surface_integral(ps,
				       gfns::gfn__global_y, true, false, true,
				       BH_diagnostics_info.integral_method);
const fp integral_z = surface_integral(ps,
				       gfns::gfn__global_z, false, true, true,
				       BH_diagnostics_info.integral_method);
const fp integral_xx = surface_integral(ps,
                                        gfns::gfn__global_xx, true, true, true,
                                        BH_diagnostics_info.integral_method);
const fp integral_xy = surface_integral(ps,
                                        gfns::gfn__global_xy, true, false, false,
                                        BH_diagnostics_info.integral_method);
const fp integral_xz = surface_integral(ps,
                                        gfns::gfn__global_xz, false, true, false,
                                        BH_diagnostics_info.integral_method);
const fp integral_yy = surface_integral(ps,
                                        gfns::gfn__global_yy, true, true, true,
                                        BH_diagnostics_info.integral_method);
const fp integral_yz = surface_integral(ps,
                                        gfns::gfn__global_yz, false, false, true,
                                        BH_diagnostics_info.integral_method);
const fp integral_zz = surface_integral(ps,
                                        gfns::gfn__global_zz, true, true, true,
                                        BH_diagnostics_info.integral_method);


//
// originds
//
this->origin_x = ps.origin_x();
this->origin_y = ps.origin_y();
this->origin_z = ps.origin_z();


//
// centroids
//
centroid_x = integral_x / integral_one;
centroid_y = integral_y / integral_one;
centroid_z = integral_z / integral_one;


//
// quadrupoles
//
quadrupole_xx = integral_xx / integral_one;
quadrupole_xy = integral_xy / integral_one;
quadrupole_xz = integral_xz / integral_one;
quadrupole_yy = integral_yy / integral_one;
quadrupole_yz = integral_yz / integral_one;
quadrupole_zz = integral_zz / integral_one;


//
// area, mean radius, and mass
//
mean_radius = integral_h / integral_one;


//
// expansion
//
area              = the_area;
expansion         = mean_expansion;
inner_expansion   = mean_inner_expansion;
product_expansion = mean_product_expansion;
mean_curvature    = mean_mean_curvature;

area_gradient              = the_area_gradient;
expansion_gradient         = mean_expansion_gradient;
inner_expansion_gradient   = mean_inner_expansion_gradient;
product_expansion_gradient = mean_product_expansion_gradient;
mean_curvature_gradient    = mean_mean_curvature_gradient;

//
// minimum, maximum and the integral of the mean curvature
//
jtutil::norm<fp> mean_curvature_norms;
ps.gridfn_norms(gfns::gfn__mean_curvature, mean_curvature_norms);
mean_curvature_minimum = mean_curvature_norms.min_value();
mean_curvature_maximum = mean_curvature_norms.max_value();
mean_curvature_integral = surface_integral(ps,
                                           gfns::gfn__mean_curvature,
                                           true, true, true,
                                           BH_diagnostics_info.integral_method);


//
// circumferences
//
circumference_xy
  = ps.circumference("xy", gfns::gfn__h,
		     gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
					 gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
							     gfns::gfn__g_dd_33,
		     BH_diagnostics_info.integral_method);
circumference_xz
  = ps.circumference("xz", gfns::gfn__h,
		     gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
					 gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
							     gfns::gfn__g_dd_33,
		     BH_diagnostics_info.integral_method);
circumference_yz
  = ps.circumference("yz", gfns::gfn__h,
		     gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
					 gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
							     gfns::gfn__g_dd_33,
		     BH_diagnostics_info.integral_method);
}

//******************************************************************************

//
// This function computes the surface integral of a gridfn over the
// horizon.
//
//static
  fp BH_diagnostics::surface_integral
	(const patch_system& ps,
	 int src_gfn, bool src_gfn_is_even_across_xy_plane,
		      bool src_gfn_is_even_across_xz_plane,
		      bool src_gfn_is_even_across_yz_plane,
	 enum patch::integration_method method)
{
return ps.integrate_gridfn
	   (src_gfn, src_gfn_is_even_across_xy_plane,
		     src_gfn_is_even_across_xz_plane,
		     src_gfn_is_even_across_yz_plane,
	    gfns::gfn__h,
	    gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
				gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
						    gfns::gfn__g_dd_33,
	    method);
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function prints a line or two summarizing the diagnostics,
// using CCTK_VInfo().
//
void BH_diagnostics::print(int N_horizons, int hn)
	const
{

const fp m_irreducible = sqrt(area / (16*PI));
CCTK_VInfo(CCTK_THORNSTRING,
	   "AH %d/%d: r=%g at (%f,%f,%f)",
	   hn, N_horizons,
	   double(mean_radius),
	   double(centroid_x), double(centroid_y), double(centroid_z));
CCTK_VInfo(CCTK_THORNSTRING,
	   "AH %d/%d: area=%.10g m_irreducible=%.10g",
	   hn, N_horizons,
	   double(area), double(m_irreducible));
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function creates a BH-diagnostics output file, writes a suitable
// header comment identifying the fields to be written by  output() ,
// flushes the stream (to help in examining the output while Cactus is
// still running), and finally returns a stdio file pointer which can be
// used by  output()  to output data to the file.
//
FILE* BH_diagnostics::setup_output_file(const struct IO_info& IO_info,
					int N_horizons, int hn)
	const
{
char file_name_buffer[IO_info::file_name_buffer_size];

const char* directory = IO_info.BH_diagnostics_directory;
// TODOMARKS
const int status = CCTK_CreateDirectory(IO_info.default_directory_permission,
 					directory);
if (status < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   BH_diagnostics::setup_output_file():\n"
"        error %d trying to create output directory\n"
"        \"%s\"!"
		   ,
		   status,
		   directory);					/*NOTREACHED*/

snprintf(file_name_buffer, IO_info::file_name_buffer_size,
	 "%s/%s.ah%d.%s",
	 directory, IO_info.BH_diagnostics_base_file_name,
	 hn, IO_info.BH_diagnostics_file_name_extension);
const char *openMode;
   openMode = "w";

FILE *fileptr = fopen(file_name_buffer, openMode);
if (fileptr == NULL)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   BH_diagnostics::setup_output_file():\n"
"        can't open BH-diagnostics output file\n"
"        \"%s\"!"
		   ,
		   file_name_buffer);				/*NOTREACHED*/
 
fprintf(fileptr, "# apparent horizon %d/%d\n", hn, N_horizons);
fprintf(fileptr, "#\n");
fprintf(fileptr, "# column  1 = cctk_iteration\n");
fprintf(fileptr, "# column  2 = cctk_time\n");
fprintf(fileptr, "# column  3 = centroid_x\n");
fprintf(fileptr, "# column  4 = centroid_y\n");
fprintf(fileptr, "# column  5 = centroid_z\n");
fprintf(fileptr, "# column  6 = min radius\n");
fprintf(fileptr, "# column  7 = max radius\n");
fprintf(fileptr, "# column  8 = mean radius\n");
fprintf(fileptr, "# column  9 = quadrupole_xx\n");
fprintf(fileptr, "# column 10 = quadrupole_xy\n");
fprintf(fileptr, "# column 11 = quadrupole_xz\n");
fprintf(fileptr, "# column 12 = quadrupole_yy\n");
fprintf(fileptr, "# column 13 = quadrupole_yz\n");
fprintf(fileptr, "# column 14 = quadrupole_zz\n");
fprintf(fileptr, "# column 15 = min x\n");
fprintf(fileptr, "# column 16 = max x\n");
fprintf(fileptr, "# column 17 = min y\n");
fprintf(fileptr, "# column 18 = max y\n");
fprintf(fileptr, "# column 19 = min z\n");
fprintf(fileptr, "# column 20 = max z\n");
fprintf(fileptr, "# column 21 = xy-plane circumference\n");
fprintf(fileptr, "# column 22 = xz-plane circumference\n");
fprintf(fileptr, "# column 23 = yz-plane circumference\n");
fprintf(fileptr, "# column 24 = ratio of xz/xy-plane circumferences\n");
fprintf(fileptr, "# column 25 = ratio of yz/xy-plane circumferences\n");
fprintf(fileptr, "# column 26 = area\n");
fprintf(fileptr, "# column 27 = m_irreducible\n");
fprintf(fileptr, "# column 28 = areal radius\n");
fprintf(fileptr, "# column 29 = expansion Theta_(l)\n");
fprintf(fileptr, "# column 30 = inner expansion Theta_(n)\n");
fprintf(fileptr, "# column 31 = product of the expansions\n");
fprintf(fileptr, "# column 32 = mean curvature\n");
fprintf(fileptr, "# column 33 = gradient of the areal radius\n");
fprintf(fileptr, "# column 34 = gradient of the expansion Theta_(l)\n");
fprintf(fileptr, "# column 35 = gradient of the inner expansion Theta_(n)\n");
fprintf(fileptr, "# column 36 = gradient of the product of the expansions\n");
fprintf(fileptr, "# column 37 = gradient of the mean curvature\n");
fprintf(fileptr, "# column 38 = minimum  of the mean curvature\n");
fprintf(fileptr, "# column 39 = maximum  of the mean curvature\n");
fprintf(fileptr, "# column 40 = integral of the mean curvature\n");
fflush(fileptr);

return fileptr;
}

//******************************************************************************

//
// This function outputs a BH-diagnostics line to a stdio stream, then
// flushes the stream (to help in examining the output while Cactus is
// still running).
//
// Arguments:
// BH_diagnostics = The BH diagnostics to be written
// fileptr = The stdio file pointer to append to
//
void BH_diagnostics::output(FILE*fileptr, const struct IO_info& IO_info)
	const
{
assert(fileptr != NULL);

fprintf(fileptr,
     //  cctk_iteration        min radius      mean radius
     //  ==  cctk_time         ======  max radius
     //  ==  ====  centroid_[xyz]      ======  ======
     //  ==  ====  ==========  ======  ======  ======
	"%d\t%.3f\t%f\t%f\t%f\t%#.10g\t%#.10g\t%#.10g\t",
	IO_info.time_iteration, double(IO_info.time),
	double(centroid_x), double(centroid_y), double(centroid_z),
	double(min_radius), double(max_radius), double(mean_radius));

fprintf(fileptr,
     //  quadrupole_{xx,xy,xz,yy,yz,zz}
     //  ================================================
	"%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t",
	double(quadrupole_xx - centroid_x*centroid_x),
        double(quadrupole_xy - centroid_x*centroid_y),
	double(quadrupole_xz - centroid_x*centroid_z),
        double(quadrupole_yy - centroid_y*centroid_y),
	double(quadrupole_yz - centroid_y*centroid_z),
        double(quadrupole_zz - centroid_z*centroid_z));

fprintf(fileptr,
     //  {min,max}_x     {min,max}_y     {min,max}_z
     //  ==============  ==============  ==============
	"%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t",
	double(min_x), double(max_x),
	double(min_y), double(max_y),
	double(min_z), double(max_z));

fprintf(fileptr,
     //  {xy,xz,yz}-plane         xz/xy  yz/xy
     //  circumferences          circumference
     //                          ratios
     //  ======================  ==============
	"%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t",
	double(circumference_xy),
	double(circumference_xz),
	double(circumference_yz),
	double(circumference_xz / circumference_xy),
	double(circumference_yz / circumference_xy));

const fp m_irreducible = sqrt(area / (16*PI));;
const fp areal_radius = sqrt(area / (4*PI));
const fp areal_radius_gradient = sqrt(1 / (16*PI*area)) * area_gradient;
fprintf(fileptr,
     //  area    m_irre- areal   expan-  inner   prod.   mean    areal   expan-  inner   prod.   mean    mean    mean    mean
     //          ducible radius  sion    expan-  of the  curva-  radius  sion    expan-  of the  curva-	 curva-  curva-  curva-
     //                                  sion    expan-  ture    grad.   grad.   sion    exp.s   ture    ture    ture    ture
     //                                          sions                           grad.   grad.   grad.   min.    max.    integ.
     //  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======
     //                                                                                                
     //  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======  ======
	"%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\n",
	double(area), double(m_irreducible), double(areal_radius),
        double(expansion), double(inner_expansion), double(product_expansion), double(mean_curvature),
        double(areal_radius_gradient),
        double(expansion_gradient), double(inner_expansion_gradient), double(product_expansion_gradient), double(mean_curvature_gradient),
	double(mean_curvature_minimum), double(mean_curvature_maximum), double(mean_curvature_integral));

fflush(fileptr);
}


//******************************************************************************
//******************************************************************************
//******************************************************************************

	  }	// namespace AHFinderDirect
