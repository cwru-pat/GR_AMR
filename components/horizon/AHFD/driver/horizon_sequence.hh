#ifndef AHFD_DRIVER_HORIZON_SEQUENCE_H
#define AHFD_DRIVER_HORIZON_SEQUENCE_H
// horizon_sequence.hh -- describes sequence of horizons a processor works on
// $Header$

// everything in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

//
// A  horizon_sequence  object describes the sequence of horizons a
// (the current) procesor works on.  This is some sequence of genuine
// horizons, followed by a dummy horizon, the latter repeating indefinitely.
//
// A horizon is specified by its (globally unique) "horizon number",
// which is 0 for a dummy horizon, or [1,N_horizons] for a genuine horizon.
//
// Thus a typical sequence of horizon numbers might be
//	1, 3, 4, 0, 0, 0, 0, ...
//
class	horizon_sequence
	{
public:
	//
	// ***** query functions *****
	//

	// how many (genuine) horizons are there in total?
	int N_horizons() const { return N_horizons_; }

	// how many genuine horizons are in the sequence (for this processor)?
	int my_N_horizons() const { return my_N_horizons_; }

	// are there any genuine horizons in the sequence (for this processor)?
	bool has_genuine_horizons() const { return my_N_horizons_ > 0; }

	// C-style string showing all horizon numbers,
	// separated by the C-style string sep, eg. "2,3,5"
	// ... result points into a private static buffer
	char* sequence_string(const char sep[]) const;

	// is the current horizon in the sequence dummy/genuine?
	bool is_dummy  () const { return posn_is_dummy  (posn_); }
	bool is_genuine() const { return posn_is_genuine(posn_); }

	// what will  is_genuine()  return after  next_hn()  is called?
	// i.e. is this *not* the final genuine  hn  in the sequence?
	bool is_next_genuine() const
		{ return posn_is_genuine( next_posn(posn_) ); }

	// return 0 if current hn is genuine,
	// or 1-origin ordinal number of dummy if it's dummy,
	//    i.e. 1 for 1st dummy, 2 for 2nd dummy, 3 for 3rd dummy, ...
	int dummy_number() const { return is_genuine() ? 0 : -posn_; }

	// get current hn in sequence
	int get_hn() const
		{ return posn_is_genuine(posn_) ? my_hn_[posn_] : 0; }

	// is a given hn genuine?
	bool is_hn_genuine(int hn) const;

	//
	// ***** traverse the sequence *****
	//
	// the idiom to traverse the genuine horizons is
	//    for (int hn = hs.init_hn() ; hs.is_genuine() ; hn = hs.next_hn())
	//    {
	//    }
	//
	// the idiom to traverse the genuine horizons, followed by
	// infinite repetition of the dummy horizon, is
	//    for (int hn = hs.init_hn() ; ; hn = hs.next_hn())
	//    {
	//    }
	//

	// reset sequence, return starting genuine hn
	int init_hn()
		{
		posn_ = (my_N_horizons_ == 0) ? -1 : 0;
		return get_hn();
		}

	// get next hn in sequence
	int next_hn() { posn_ = next_posn(posn_); return get_hn(); }


	//
	// ***** set up the sequence *****
	//

	// construct an empty horizon sequence
	// which can hold <= N_horizons horizon numbers
	horizon_sequence(int N_horizons);
	~horizon_sequence();

	// append hn to the sequence
	// ... returns new my_N_horizons()
	int append_hn(int hn);

private:
	// is a specified posn dummy/genuine
	bool posn_is_genuine(int pos) const
		{ return (pos >= 0) && (pos < my_N_horizons_); }
	bool posn_is_dummy(int pos) const
		{ return !posn_is_genuine(pos); }

	// what is the next posn in the sequence
	int next_posn(int pos) const;

private:
	// we forbid copying and passing by value
	// by declaring the copy constructor and assignment operator
	// private, but never defining them
	horizon_sequence(const horizon_sequence& rhs);
	horizon_sequence& operator=(const horizon_sequence& rhs);

private:
	const int N_horizons_;
	int my_N_horizons_;

	// "internal position" in sequence
	// this is a [0,my_N_horizons_) index into my_hn_[]
	//           for genuine horizons,
	//         or < 0 for the dummy horizon, with the absolute value
	//            counting the number of times we've returned hn=0
	int posn_;

	int* my_hn_;	// --> new[]-allocated array of genuine horizon numbers
	};

//******************************************************************************

	  }	// namespace AHFinderDirect
#endif
