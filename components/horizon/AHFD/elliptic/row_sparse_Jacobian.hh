// Jacobian.hh -- data structures for the Jacobian matrix
// $Header$
//
#ifdef HAVE_ROW_SPARSE_JACOBIAN
// row_sparse_Jacobian -- Jacobian stored as row-oriented sparse matrix
#endif
#ifdef HAVE_ROW_SPARSE_JACOBIAN__ILUCG
// row_sparse_Jacobian__ILUCG -- ... with ILUCG linear solver
#endif
#ifdef HAVE_ROW_SPARSE_JACOBIAN__UMFPACK
// row_sparse_Jacobian__UMFPACK -- ... with UMFPACK linear solver
#endif
//

#ifndef AHFINDERDIRECT__ROW_SPARSE_JACOBIAN_HH
#define AHFINDERDIRECT__ROW_SPARSE_JACOBIAN_HH

//
// prerequisites:
//	"../patch/patch_system.hh"
//	"Jacobian.hh"
//

#undef  DEBUG_ROW_SPARSE_JACOBIAN	// define this for detailed
					// debugging of data structures
					// (quite slow, not for production use)
#ifdef DEBUG_ROW_SPARSE_JACOBIAN
  // define this to sanity-check the matrix data structures
  // at the start of each new row
  #define ROW_SPARSE_JACOBIAN__CHECK_AT_ROW_START

  // define this to sanity-check the matrix data structures
  // after each new element is inserted
  #undef  ROW_SPARSE_JACOBIAN__CHECK_AFTER_INSERT

  // define this to print a line or two of information
  // for each of the main data-structure operations,
  // e.g. zero the matrix, start a new row, insert a new matrix element, etc
  // ... note this produces a lot of output
  //     (tens of megabytes for a reasonable-sized angular grid)
  #define ROW_SPARSE_JACOBIAN__LOG_MAIN_OPS

  // define this to print and sanity-check the matrix data structures
  // at the start of each new row; note this produces a *huge* amount of
  // output (many hundreds of megabytes for a reasonable-sized angular grid),
  // and is only suitable for tiny examples
  #undef  ROW_SPARSE_JACOBIAN__PRINT_AT_ROW_START

  // define this to print and sanity-check the matrix data structures
  // after each new element is inserted; note this produces even more
  // output than "just" printing at the start of each row!
  #undef  ROW_SPARSE_JACOBIAN__PRINT_AFTER_INSERT
#endif

// everything in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN
//
// This class stores the Jacobian as a row-oriented sparse matrix.
// This sparse matrix format is widely used by sparse matrix libraries.
//
// There are two variants of the format, depending on whether the array
// indices are C-style 0-origin or Fortran-style 1-origin values.  Note
// that these variants actually store different values in memory -- the
// difference is *not* just an adjustment of the subscripting.  To handle
// both variants with common code, we store an integer "index origin" IO,
// which is 0 for C, 1 for Fortran.
//
// The matrix is represented by 3 arrays:
//	IA[N_rows+1] = IO-relative indices in JA and A of the start of
//		       each row's matrix elements;
//			IA(IO       ) = IO
//			IA(IO+N_rows) = IO + N_nonzeros
//	JA[N_nonzeros] = IO-relative column indices of the corresponding
//			 entries in A
//	A[N_nonzeros] = nonzero matrix elements, ordered by rows;
//			within a row the matrix elements may either be
//			in arbitrary order, or may be sorted in increasing
//			order of their column indices
// A[] may include explicitly stored zeros if necessary.
// IA[] and JA[] are arrays of  integer  (a C typedef which matches the
// Fortran "integer" datatype) so they can be passed by reference to/from
// Fortran code if necessary.
//
// For example, the matrix
//	[ 11.   12.   13.             ]
//	[ 21.   22.         24.   25. ]
//	[       32.         34.       ]
//	[                   44.   45. ]
//	[             53.         55. ]
// could be represented by the IO=0 arrays
//	IA[] = {  0,           3,               7,       9,      11,     13 }
//	JA[] = {  0,  1,  2,   0,  1,  3,  4,   1,  3,   3,  4,   2,  4 }
//	A [] = { 11.,12.,13., 21.,22.,24.,25., 32.,34., 44.,45., 53.,55.}
// or the IO=1 arrays
//	IA[] = {  1,           4,               8,      10,      12,     14 }
//	JA[] = {  1,  2,  3,   1,  2,  4,  5,   2,  4,   4,  5,   3,  5 }
//	A [] = { 11.,12.,13., 21.,22.,24.,25., 32.,34., 44.,45., 53.,55.}
// Notice that the A[] array is identical in both cases, but the IA[] and
// JA[] arrays vary.
//
class	row_sparse_Jacobian
	: public Jacobian
	{
public:
	//
	// index origin can be either C or Fortran
	//

	int IO() const { return IO_; }
	enum {C_index_origin = 0, Fortran_index_origin = 1};


	//
	// routines to access the matrix
	//

	// get a matrix element
	fp element(int II, int JJ) const;

	// is the matrix element (II,JJ) stored explicitly?
	bool is_explicitly_stored(int II, int JJ) const
		{ return find_element(II,JJ) > 0; }


	//
	// routines for setting values in the matrix
	// ... with the current implementation, clients *must* set up
	//     the matrix in row order.
	//

	// zero the entire matrix
	void zero_matrix();

	// set a matrix element to a specified value
	void set_element(int II, int JJ, fp value);

	// sum a value into a matrix element
	void sum_into_element(int II, int JJ, fp value);


	//
	// internal routines
	//
protected:
	// return 0-origin position of (II,JJ) in A_ and JA_,
	// or -1 if not found
	int find_element(int II, int JJ) const;

	// insert new array element (II,JJ) in matrix,
	// starting new row and/or growing arrays if necessary
	// return 0-origin position (1-origin) in A_ and JA_
	int insert_element(int II, int JJ, fp value);

	// grow arrays to make room for more elements
	// ... growth is geometric ==> amortized cost remains O(current size)
	void grow_arrays();		// this fn does the actual growing

	// parameter for growth policy
	// ... this should be > 0
	// ... for debugging it may be useful to set this to a fairly
	//     small value (even 1), to force some  grow_array()  calls
	//     when the arrays are still small enough for easy examination
	enum {base_growth_amount = 1000};	// grow by this much
						// plus half the current size

	// sort each row's JA_[] and A_[] array elements
	// so the columns (JA_[] values) are in increasing order
	// ... this doesn't change the abstract value of the matrix
	void sort_each_row_into_column_order();

#ifdef DEBUG_ROW_SPARSE_JACOBIAN
	// sanity-check and optionally print matrix data structure
	// (i.e. the IA_ and JA_ arrays)
	void check_and_print_data_structure(bool print_flag) const;
#endif


	//
	// constructor, destructor
	//
protected:
	// the constructor only uses ps to get the size of the matrix
	// ... print_msg_flag controls messages about data structure setup
	//     during construction and element-setting process
	row_sparse_Jacobian(patch_system& ps, int IO_in,
			    bool print_msg_flag = false);
	~row_sparse_Jacobian();


protected:
	int IO_;		// index origin (0=C, 1=Fortran)
#ifdef DEBUG_ROW_SPARSE_JACOBIAN
	const char*const IOstr_;// --> "c" or "f" as appropriate
#endif

	int N_nonzeros_;	// current number of nonzeros in matrix
				// = number of elements valid in JA_[], A_[]
	int current_N_rows_;	// current number of rows in matrix
				// i.e. current rows have 0-origin indices
				//      0 <= II < current_N_rows_

	int N_nonzeros_allocated_;	// allocated size of JA_[], A_[]

	// --> new[]-allocated arrays
	// ***** due to constructor initialization list ordering,
	// ***** these must be declared *AFTER* N_nonzeros_allocated_

	// this is allocated in the constructor
	integer* IA_;		// --> array of size N_rows_+1
				// valid for 0-origin indices
				//       0 <= II <= current_N_rows_

	// this is allocated in the constructor,
	// but may be reallocated by  grow_arrays()
	// ... using STL vector would be much cleaner here, but alas
	//     there are too many broken compilers/linkers in the world
	//     world which don't fully grok vector :( :(
	integer* JA_;		// --> array of size N_nonzeros_allocated_
				// valid for 0-origin indices
				//       0 <= posn < II[current_N_rows]
	fp* A_;			// ditto
	};
#endif	/* HAVE_ROW_SPARSE_JACOBIAN */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN__ILUCG
//
// This class defines the linear solver routine using ILUCG.
//
// This class stores the matrix with IO=1 (Fortran-style indices),
// and does not sort the elements in a row into increasing-column order.
//
class row_sparse_Jacobian__ILUCG
	: public row_sparse_Jacobian
	{
public:
	// solve linear system J.x = rhs via ILUCG
	// ... rhs and x are nominal-grid gridfns
	// ... overwrites Jacobian matrix with its LU decomposition
	// ... returns -1.0 to signal that condition number is unknown
	fp solve_linear_system(int rhs_gfn, int x_gfn,
			       const struct linear_solver_pars& pars,
			       bool print_msg_flag);

	// constructor, destructor
public:
	// the constructor only uses ps to get the size of the matrix
	row_sparse_Jacobian__ILUCG(patch_system& ps,
				   bool print_msg_flag = false);
	~row_sparse_Jacobian__ILUCG();

private:
	// work vectors for ILUCG subroutine
	// ... allocated by  solve_linear_system()  on first call
	integer* itemp_;
	fp*      rtemp_;
	};
#endif	/* HAVE_ROW_SPARSE_JACOBIAN__ILUCG */

//******************************************************************************

#ifdef HAVE_ROW_SPARSE_JACOBIAN__UMFPACK
//
// This class defines the linear solver routine using UMFPACK.
//
// This class stores the matrix with IO=0 (C-style indices),
// and sorts the elements in a row into increasing-column order.
//
// Notice that due to the pImpl-style design of the UMFPACK classes,
// and our using Fortran "integer" instead of UMFPACK "Int" for the
// integer arrays, "umfpack.h" is *not* a prerequisite for including
// this header file.
//
class row_sparse_Jacobian__UMFPACK
	: public row_sparse_Jacobian
	{
public:
	// solve linear system J.x = rhs via UMFPACK
	// ... rhs and x are nominal-grid gridfns
	// ... does symbolic factorization on first call,
	//     reuses this on all following calls
	// ... returns estimated reciprocal condition number
	//	       ... rcond = 0.0         ==> exactly singular
	//		   rcond < DBL_EPSILON ==> numerically singular
	//		   rcond = 1.0         ==> orthogonal matrix
	fp solve_linear_system(int rhs_gfn, int x_gfn,
			       const struct linear_solver_pars& pars,
			       bool print_msg_flag);

	// constructor, destructor
public:
	// the constructor only uses ps to get the size of the matrix
	row_sparse_Jacobian__UMFPACK(patch_system& ps,
				     bool print_msg_flag = false);
	~row_sparse_Jacobian__UMFPACK();

private:
	// pointers to UMFPACK control arrays
	double *Control_;
	double *Info_;

	// pointers to UMFPACK workspace arrays
	// ... allocated by  solve_linear_system()  on first call
	integer *solve_workspace_integer_;
	double  *solve_workspace_double_;

	// pointers to UMFPACK data objects
	// (UMFPACK takes these as  void*  and casts them internally)
	void* Symbolic_;
	void* Numeric_;

	// 
	};
#endif	/* HAVE_ROW_SPARSE_JACOBIAN__UMFPACK */

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif		/* AHFINDERDIRECT__ROW_SPARSE_JACOBIAN_HH */
