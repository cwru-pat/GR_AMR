// dense_Jacobian.hh -- dense-matrix Jacobian 
// $Header$
//
#ifdef HAVE_DENSE_JACOBIAN
// dense_Jacobian -- Jacobian stored as a dense matrix
#endif
#ifdef HAVE_DENSE_JACOBIAN__LAPACK
// dense_Jacobian__LAPACK -- dense_Jacobian with LAPACK linear solver
#endif
//

#ifndef AHFINDERDIRECT__DENSE_JACOBIAN_HH
#define AHFINDERDIRECT__DENSE_JACOBIAN_HH

//
// prerequisites:
//	"../patch/patch_system.hh"
//	"Jacobian.hh"
//

// everything in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

#ifdef HAVE_DENSE_JACOBIAN
//
// This class stores the Jacobian as a dense matrix in Fortran (column)
// order.
//
class	dense_Jacobian
	: public Jacobian
	{
public:
	//
	// routines to access the matrix
	//

	// get a matrix element
	fp element(int II, int JJ) const
		{ return matrix_(JJ,II); }

	// dense matrix ==> all elements are explicitly stored
	bool is_explicitly_stored(int II, int JJ) const
		{ return true; }


	//
	// routines for setting values in the matrix
	//

	// zero the entire matrix
	void zero_matrix();

	// set a matrix element to a specified value
	void set_element(int II, int JJ, fp value)
		{ matrix_(JJ,II) = value; }

	// sum a value into a matrix element
	void sum_into_element(int II, int JJ, fp value)
		{ matrix_(JJ,II) += value; }

	//
	// constructor, destructor
	//
protected:
	// the constructor only uses ps to get the size of the matrix
	// ... print_msg_flag controls messages about data structure setup
	//     during construction and element-setting process
	dense_Jacobian(patch_system& ps,
		       bool print_msg_flag = false);
	~dense_Jacobian() { }

protected:
	// Fortran storage order ==> subscripts are (JJ,II)
	jtutil::array2d<fp> matrix_;
	};
#endif	/* HAVE_DENSE_JACOBIAN */

//******************************************************************************

#ifdef HAVE_DENSE_JACOBIAN__LAPACK
//
// This class defines the linear solver routine using LAPACK.
//
class	dense_Jacobian__LAPACK
	: public dense_Jacobian
	{
public:
	// solve linear system J.x = rhs via LAPACK
	// ... rhs and x are nominal-grid gridfns
	// ... overwrites Jacobian matrix with its LU decomposition
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
	dense_Jacobian__LAPACK(patch_system& ps,
			       bool print_msg_flag = false);
	~dense_Jacobian__LAPACK();

private:
	// pivot vector for LAPACK routines
	integer *pivot_;	// size N_rows_

	// work vectors for LAPACK condition number computation
	integer *iwork_;	// size N_rows_
	fp *rwork_;		// size 4*N_rows_
	};
#endif	/* HAVE_DENSE_JACOBIAN__LAPACK */

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif		/* AHFINDERDIRECT__DENSE_JACOBIAN_HH */
