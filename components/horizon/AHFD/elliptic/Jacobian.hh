// Jacobian.hh -- generic data structures for the Jacobian matrix
// $Header$
//
// Jacobian -- ABC to describe Jacobian matrix
// decode_Jacobian_store_solve_method - decode string into internal enum
// new_Jacobian - factory method
//

#ifndef AHFINDERDIRECT__JACOBIAN_HH
#define AHFINDERDIRECT__JACOBIAN_HH

//
// prerequisites:
//	"../patch/patch_system.hh"
//

// everything in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

//
// These classes are used to store and manipulate Jacobian matrices
// for a patch system.
//
// Jacobian is an abstract base class (ABC) that defines the generic
// Jacobian API.  We derive classes from it for each type of sparsity
// pattern, then we derive classes from them for each linear solver.
//
// At present the inheritance graph looks like this:
//	Jacobian
//	    dense_Jacobian
//		dense_Jacobian__LAPACK
//	    row_sparse_Jacobian
//		row_sparse_Jacobian__ILUCG
//		row_sparse_Jacobian__UMFPACK
// each derived class is inside a corresponding #ifdef (set or unset
// as appropriate, in "../include/config.h").
//

//******************************************************************************

//
// We assume a "traversal" interface for computing Jacobians, defined
// by the  Jacobian  API.  The typical code to construct and use a
// Jacobian matrix looks like this:
//	// prototype
//	// ... calls Jac.zero_matrix()
//	//           Jac.set_element()
//	//           Jac.sum_into_element()
//	void traverse_Jacobian(const patch_system& ps, Jacobian& Jac);
//
//	Jacobian& Jac = new_Jacobian(Jac_type, ps);
//	traverse_Jacobian(ps, Jac);
//	Jac.solve_linear_system(...);
//
//	Jac.zero_matrix()
//	traverse_Jacobian(ps, Jac);	// must specify same sparsity pattern
//					// as before
//	Jac.solve_linear_system(...);
//	
//

//
// A row/column index of the Jacobian (denoted II/JJ) is a 0-origin grid
// point number within the patch system.
//
// Since many of the derived classes use Fortran routines, we also use
// 1-origin indices; these have a leading "f", eg fII/fJJ.  Finally, we
// use generic indices that might be either 0-origin or 1-origin; these
// have a leading "g", eg gII/gJJ.
//

//
// Note that the APIs here implicitly assume there is only a *single* gridfn
// in the Jacobian computation.  (If we had N gridfns for this, then the
// Jacobian would really be a block matrix with N*N blocks.)  This isn't
// very general, but matches our present use in this thorn.
//

//******************************************************************************

// forward declarations
class linear_solver_pars;

// ABC to define Jacobian matrix
class	Jacobian
	{
public:
	// basic meta-info
	patch_system& my_patch_system() const { return ps_; }
	int N_rows() const { return N_rows_; }

	// convert (patch,irho,isigma) <--> row/column index
	int II_of_patch_irho_isigma(const patch& p, int irho, int isigma)
		const
		{ return ps_.gpn_of_patch_irho_isigma(p, irho,isigma); }
	const patch& patch_irho_isigma_of_II(int II, int& irho, int& isigma)
		const
		{ return ps_.patch_irho_isigma_of_gpn(II, irho,isigma); }


#ifdef NOT_USED
	//
	// convert C <--> Fortran indices
	//
	int csub(int f) const { return f-1; }
	int fsub(int c) const { return c+1; }
#endif


	//
	// routines to access the matrix
	//

	// get a matrix element
	virtual fp element(int II, int JJ) const = 0;

	// is a given element explicitly stored, or implicitly 0 via sparsity
	virtual bool is_explicitly_stored(int II, int JJ) const = 0;


	//
	// routines for setting values in the matrix
	//

	// zero the entire matrix
	virtual void zero_matrix() = 0;

	// set a matrix element to a specified value
	virtual void set_element(int II, int JJ, fp value) = 0;

	// sum a value into a matrix element
	virtual void sum_into_element(int II, int JJ, fp value) = 0;


	//
	// solve linear system J.x = rhs
	// ... rhs and x are nominal-grid gridfns
	// ... may modify Jacobian matrix (eg for LU decomposition)
	// ... returns estimated reciprocal condition number if known,
	//	       -1.0 if condition number is unknown
	//	       ... rcond = 0.0         ==> exactly singular
	//		   rcond < DBL_EPSILON ==> numerically singular
	//		   rcond = 1.0         ==> orthogonal matrix
	// ... once this has been called, the sparsity pattern should
	//     not be changed, i.e. no new nonzeros should be introduced
	//     into the matrix
	//
	virtual fp solve_linear_system(int rhs_gfn, int x_gfn,
				       const struct linear_solver_pars& pars,
				       bool print_msg_flag) = 0;


	//
	// constructor, destructor
	//
protected:
	// the constructor only uses  ps  to get the size of the matrix
	Jacobian(patch_system& ps)
		: ps_(ps),
		  N_rows_(ps.N_grid_points() + ps.N_additional_points())
		{ }
public:
	virtual ~Jacobian() { }

private:
	// we forbid copying and passing by value
	// by declaring the copy constructor and assignment operator
	// private, but never defining them
	Jacobian(const Jacobian &rhs);
	Jacobian& operator=(const Jacobian &rhs);

protected:
	patch_system& ps_;
	int N_rows_;
	};

//**************************************

//
// this class defines parameters for  solve_linear_system()
// for all of our derived classes
//
struct	linear_solver_pars
	{
	struct	LAPACK_pars
		{
		} LAPACK_pars;
	struct	ILUCG_pars
		{
		fp  error_tolerance;
		// should we limit to N_rows_ conjugate gradient iterations?
		bool limit_CG_iterations;
		} ILUCG_pars;
	struct	UMFPACK_pars
		{
		// number of iterative-improvement iterations to do
		// inside UMFPACK after the sparse LU decompose/solve,
		// each time we solve a linear system
		int N_II_iterations;
		} UMFPACK_pars;
	};

//******************************************************************************

enum	Jacobian_store_solve_method
	{
#ifdef HAVE_DENSE_JACOBIAN__LAPACK
	Jacobian__dense_matrix__LAPACK,
#endif
#ifdef HAVE_ROW_SPARSE_JACOBIAN__ILUCG
	Jacobian__row_sparse_matrix__ILUCG,
#endif
#ifdef HAVE_ROW_SPARSE_JACOBIAN__UMFPACK
	Jacobian__row_sparse_matrix__UMFPACK // no comma on last entry in enum
#endif
	};

//
// prototypes of various Jacobian-related functions
//

// decode string into internal enum
enum Jacobian_store_solve_method
  decode_Jacobian_store_solve_method
	(const char Jacobian_store_solve_method_string[]);

// "object factory" routines to construct and return
// pointers to new-allocated objects of specified derived types
Jacobian* new_Jacobian(enum Jacobian_store_solve_method Jac_method,
		       patch_system& ps,
		       bool print_msg_flag = false);

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif		/* AHFINDERDIRECT__JACOBIAN_HH */
