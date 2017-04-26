#ifndef COSMO_TYPES
#define COSMO_TYPES

//#include "cosmo_includes.h"
//#include "cosmo_macros.h"
#include "SAMRAI/pdat/MDA_Access.h"

namespace cosmo
{

#define DIM 3
#define EPS (1e-11)
#define INF 1e100
// changing this affects FFTs:
typedef double real_t; /**< real type; changing this may require changes to HDF5 and FFTW functionality */
// see http://www.fftw.org/doc/Precision.html

typedef int idx_t; /**< indexing type, must be long enough to support large arrays */

typedef MDA_Access<double, DIM, MDA_OrderColMajor<DIM>> arr_t; /**< base array type */


} /* namespace cosmo */

#endif
