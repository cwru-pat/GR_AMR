/* structs.h -- definitions of data/coeffsicient structures */
/* $Header$ */

/*
 * prerequisite headers:
 *	"cctk.h"
 *	"../InterpLocalUniform.h"
 */

/******************************************************************************/

/*
 * Each of these structures holds a molecule-sized piece of a single
 * real data array, or of a single real/imaginary component of a single
 * complex data array.
 */

struct	data_struct_1d_cube_size2
	{
	#include "1d.cube.size2/data-dcl.h"
	};
struct	data_struct_1d_cube_size3
	{
	#include "1d.cube.size3/data-dcl.h"
	};
struct	data_struct_1d_cube_size4
	{
	#include "1d.cube.size4/data-dcl.h"
	};
struct	data_struct_1d_cube_size5
	{
	#include "1d.cube.size5/data-dcl.h"
	};
struct	data_struct_1d_cube_size6
	{
	#include "1d.cube.size6/data-dcl.h"
	};
struct	data_struct_1d_cube_size7
	{
	#include "1d.cube.size7/data-dcl.h"
	};

struct	data_struct_2d_cube_size2
	{
	#include "2d.cube.size2/data-dcl.h"
	};
struct	data_struct_2d_cube_size3
	{
	#include "2d.cube.size3/data-dcl.h"
	};
struct	data_struct_2d_cube_size4
	{
	#include "2d.cube.size4/data-dcl.h"
	};
struct	data_struct_2d_cube_size5
	{
	#include "2d.cube.size5/data-dcl.h"
	};
struct	data_struct_2d_cube_size6
	{
	#include "2d.cube.size6/data-dcl.h"
	};

struct	data_struct_3d_cube_size2
	{
	#include "3d.cube.size2/data-dcl.h"
	};
struct	data_struct_3d_cube_size3
	{
	#include "3d.cube.size3/data-dcl.h"
	};
struct	data_struct_3d_cube_size4
	{
	#include "3d.cube.size4/data-dcl.h"
	};
struct	data_struct_3d_cube_size5
	{
	#include "3d.cube.size5/data-dcl.h"
	};
struct	data_struct_3d_cube_size6
	{
	#include "3d.cube.size6/data-dcl.h"
	};

/******************************************************************************/

/*
 * Each of these structures holds all the coeffsicients for a single
 * molecule (written as a linear combination of the input data).
 */

struct	coeffs_struct_1d_cube_size2
	{
	#include "1d.cube.size2/coeffs-dcl.h"
	};
struct	coeffs_struct_1d_cube_size3
	{
	#include "1d.cube.size3/coeffs-dcl.h"
	};
struct	coeffs_struct_1d_cube_size4
	{
	#include "1d.cube.size4/coeffs-dcl.h"
	};
struct	coeffs_struct_1d_cube_size5
	{
	#include "1d.cube.size5/coeffs-dcl.h"
	};
struct	coeffs_struct_1d_cube_size6
	{
	#include "1d.cube.size6/coeffs-dcl.h"
	};
struct	coeffs_struct_1d_cube_size7
	{
	#include "1d.cube.size7/coeffs-dcl.h"
	};

struct	coeffs_struct_2d_cube_size2
	{
	#include "2d.cube.size2/coeffs-dcl.h"
	};
struct	coeffs_struct_2d_cube_size3
	{
	#include "2d.cube.size3/coeffs-dcl.h"
	};
struct	coeffs_struct_2d_cube_size4
	{
	#include "2d.cube.size4/coeffs-dcl.h"
	};
struct	coeffs_struct_2d_cube_size5
	{
	#include "2d.cube.size5/coeffs-dcl.h"
	};
struct	coeffs_struct_2d_cube_size6
	{
	#include "2d.cube.size6/coeffs-dcl.h"
	};

struct	coeffs_struct_3d_cube_size2
	{
	#include "3d.cube.size2/coeffs-dcl.h"
	};
struct	coeffs_struct_3d_cube_size3
	{
	#include "3d.cube.size3/coeffs-dcl.h"
	};
struct	coeffs_struct_3d_cube_size4
	{
	#include "3d.cube.size4/coeffs-dcl.h"
	};
struct	coeffs_struct_3d_cube_size5
	{
	#include "3d.cube.size5/coeffs-dcl.h"
	};
struct	coeffs_struct_3d_cube_size6
	{
	#include "3d.cube.size6/coeffs-dcl.h"
	};
