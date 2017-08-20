// horizon_sequence.cc -- describes sequence of horizons a processor works on
// $Header$

//
// horizon_sequence::horizon_sequence
// horizon_sequence::~horizon_sequence
// horizon_sequence::sequence_string
// horizon_sequence::append_hn
// horizon_sequence::next_posn
// horizon_sequence::is_hn_genuine
//

#include <stdio.h>
#include <assert.h>
#include <cmath>

#include "../jtutil/util_String.h"
#include "../jtutil/util.hh"
#include "../AHFD_macros.h"

#include "horizon_sequence.hh"

// all the code in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

//
// This function constructs an horizon sequence, which can hold a maximum
// of N_horizons horizon numbers.
//
horizon_sequence::horizon_sequence(int N_horizons_in)
	: N_horizons_(N_horizons_in),
	  my_N_horizons_(0),		// sequence starts out empty
	  posn_(-1),
	  my_hn_(new int[N_horizons_in])
{ }

//******************************************************************************

//
// This function destroys a horizon sequence.
//
horizon_sequence::~horizon_sequence()
{
delete[] my_hn_;
}

//******************************************************************************

//
// This function returns a (pointer to a) C-style string showing all
// the horizon numbers in the sequence, separated by the C-style string
//  sep , eg. "2,3,5".
//
// Note the result points into a private static buffer.
//
char* horizon_sequence::sequence_string(const char sep[]) const
{
const int N_hn_buffer = 10;
char hn_buffer[N_hn_buffer];

const int N_buffer = 100;
static char buffer[N_buffer];

buffer[0] = '\0';
	for (int pos = 0 ; pos < my_N_horizons_ ; ++pos)
	{
	if (pos > 0)
	   then Util_Strlcat(buffer, sep, N_buffer);
	snprintf(hn_buffer, N_hn_buffer, "%d", my_hn_[pos]);
	Util_Strlcat(buffer, hn_buffer, N_buffer);
	}

return buffer;
}

//******************************************************************************

//
// This function appends  hn  to the sequence.  It returns the new value
// of my_N_horizons().
//
int horizon_sequence::append_hn(int hn)
{
assert( hn > 0 );			// can only append genuine horizons
assert( my_N_horizons_ < N_horizons_ );	// make sure there's space for it
my_hn_[my_N_horizons_++] = hn;
posn_ = 0;
return my_N_horizons_;
}

//******************************************************************************

//
// This function computes the internal position immediately following
// a given internal position in the sequence.
//
// Arguments:
// p = (in) The current internal position, with posn_ semantics
//
// Results:
// This function returns the next internal position after p.
//
int horizon_sequence::next_posn(int pos)
	const
{
return   (pos < 0) ? pos-1
       : (pos+1 < my_N_horizons_) ? pos+1
       : -1;
}

//******************************************************************************

//
// This function determines whether or not a given  hn  is genuine.
//
bool horizon_sequence::is_hn_genuine(int hn)
	const
{
	for (int pos = 0 ; pos < my_N_horizons_ ; ++pos)
	{
	if (my_hn_[pos] == hn)
	   then return true;
	}

return false;
}

//******************************************************************************

	  }	// namespace AHFinderDirect
