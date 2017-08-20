#ifndef AHFD_GR_CG_H
#define AHFD_GR_CG_H

// cg.hh -- interface between machine-generated and other code
// $Header$

//
// This header file defines the "virtual machine" used by machine-generated
// code.  It is "dangerous" in that it #defines macros and local variables
// for all the gridfns, which will break lots of other code.  Thus this file
// can only be #included within a function, and should only be #included
// just before #including machine-generated code, and this should be as
// close as possible to the end of the source file (to minimize the amount
// of other code affected by this file).
//
// FIXME: we should have an "uncg.hh" header file to #undef all these macros
//
// prerequisites
//	gfns.hh
//

//******************************************************************************

//
// The machine-generated code uses the following variables:
//	patch& p		// current patch
//	int irho,isigma		// current generic integer coords within patch
//	fp xx, yy, zz;		// local (x,y,z) coords of this grid point
//	fp X_ud_[12][123]			// 1st derivative coefficients
//	fp X_udd_[12]{11,12,13,22,23,33}	// 2nd derivative coefficients
//

//
// rational numbers
//
#define RATIONAL(num,den)	(num/den)

//
// partial derivatives
//
#define PARTIAL_RHO(ghosted_gridfn_name)	\
	p.partial_rho(gfns::gfn__ ## ghosted_gridfn_name, irho,isigma)
#define PARTIAL_SIGMA(ghosted_gridfn_name)	\
	p.partial_sigma(gfns::gfn__ ## ghosted_gridfn_name, irho,isigma)
#define PARTIAL_RHO_RHO(ghosted_gridfn_name)	\
	p.partial_rho_rho(gfns::gfn__ ## ghosted_gridfn_name, irho,isigma)
#define PARTIAL_RHO_SIGMA(ghosted_gridfn_name)	\
	p.partial_rho_sigma(gfns::gfn__ ## ghosted_gridfn_name, irho,isigma)
#define PARTIAL_SIGMA_SIGMA(ghosted_gridfn_name)\
	p.partial_sigma_sigma(gfns::gfn__ ## ghosted_gridfn_name, irho,isigma)

//******************************************************************************

//
// ghosted-grid gridfns
// n.b. since we're always evaluating on the surface, r == h
//
#define h	p.ghosted_gridfn(gfns::gfn__h, irho,isigma)
#define r	h

//******************************************************************************

//
// nominal-grid gridfns
//
#define g_dd_11	p.gridfn(gfns::gfn__g_dd_11, irho,isigma)
#define g_dd_12	p.gridfn(gfns::gfn__g_dd_12, irho,isigma)
#define g_dd_13	p.gridfn(gfns::gfn__g_dd_13, irho,isigma)
#define g_dd_22	p.gridfn(gfns::gfn__g_dd_22, irho,isigma)
#define g_dd_23	p.gridfn(gfns::gfn__g_dd_23, irho,isigma)
#define g_dd_33	p.gridfn(gfns::gfn__g_dd_33, irho,isigma)
#define K_dd_11	p.gridfn(gfns::gfn__K_dd_11, irho,isigma)
#define K_dd_12	p.gridfn(gfns::gfn__K_dd_12, irho,isigma)
#define K_dd_13	p.gridfn(gfns::gfn__K_dd_13, irho,isigma)
#define K_dd_22	p.gridfn(gfns::gfn__K_dd_22, irho,isigma)
#define K_dd_23	p.gridfn(gfns::gfn__K_dd_23, irho,isigma)
#define K_dd_33	p.gridfn(gfns::gfn__K_dd_33, irho,isigma)

#define partial_d_g_dd_111	p.gridfn(gfns::gfn__partial_d_g_dd_111, irho,isigma)
#define partial_d_g_dd_112	p.gridfn(gfns::gfn__partial_d_g_dd_112, irho,isigma)
#define partial_d_g_dd_113	p.gridfn(gfns::gfn__partial_d_g_dd_113, irho,isigma)
#define partial_d_g_dd_122	p.gridfn(gfns::gfn__partial_d_g_dd_122, irho,isigma)
#define partial_d_g_dd_123	p.gridfn(gfns::gfn__partial_d_g_dd_123, irho,isigma)
#define partial_d_g_dd_133	p.gridfn(gfns::gfn__partial_d_g_dd_133, irho,isigma)
#define partial_d_g_dd_211	p.gridfn(gfns::gfn__partial_d_g_dd_211, irho,isigma)
#define partial_d_g_dd_212	p.gridfn(gfns::gfn__partial_d_g_dd_212, irho,isigma)
#define partial_d_g_dd_213	p.gridfn(gfns::gfn__partial_d_g_dd_213, irho,isigma)
#define partial_d_g_dd_222	p.gridfn(gfns::gfn__partial_d_g_dd_222, irho,isigma)
#define partial_d_g_dd_223	p.gridfn(gfns::gfn__partial_d_g_dd_223, irho,isigma)
#define partial_d_g_dd_233	p.gridfn(gfns::gfn__partial_d_g_dd_233, irho,isigma)
#define partial_d_g_dd_311	p.gridfn(gfns::gfn__partial_d_g_dd_311, irho,isigma)
#define partial_d_g_dd_312	p.gridfn(gfns::gfn__partial_d_g_dd_312, irho,isigma)
#define partial_d_g_dd_313	p.gridfn(gfns::gfn__partial_d_g_dd_313, irho,isigma)
#define partial_d_g_dd_322	p.gridfn(gfns::gfn__partial_d_g_dd_322, irho,isigma)
#define partial_d_g_dd_323	p.gridfn(gfns::gfn__partial_d_g_dd_323, irho,isigma)
#define partial_d_g_dd_333	p.gridfn(gfns::gfn__partial_d_g_dd_333, irho,isigma)

#define Theta	p.gridfn(gfns::gfn__Theta, irho,isigma)

#define partial_Theta_wrt_partial_d_h_1	\
	p.gridfn(gfns::gfn__partial_Theta_wrt_partial_d_h_1, irho,isigma)
#define partial_Theta_wrt_partial_d_h_2	\
	p.gridfn(gfns::gfn__partial_Theta_wrt_partial_d_h_2, irho,isigma)
#define partial_Theta_wrt_partial_dd_h_11	\
	p.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_11, irho,isigma)
#define partial_Theta_wrt_partial_dd_h_12	\
	p.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_12, irho,isigma)
#define partial_Theta_wrt_partial_dd_h_22	\
	p.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_22, irho,isigma)

#define save_Theta	p.gridfn(gfns::gfn__save_Theta, irho,isigma)
#define old_Theta	p.gridfn(gfns::gfn__old_Theta, irho,isigma)
#define Delta_h		p.gridfn(gfns::gfn__Delta_h, irho,isigma)

//******************************************************************************

//
// pseudo-gridfns, i.e. temporaries used only at a single grid point
// (these are not actually stored as gridfns)
//
fp g_uu_11;
fp g_uu_12;
fp g_uu_13;
fp g_uu_22;
fp g_uu_23;
fp g_uu_33;
fp K;
fp K_uu_11;
fp K_uu_12;
fp K_uu_13;
fp K_uu_22;
fp K_uu_23;
fp K_uu_33;

fp partial_d_ln_sqrt_g_1;
fp partial_d_ln_sqrt_g_2;
fp partial_d_ln_sqrt_g_3;

fp partial_d_g_uu_111;
fp partial_d_g_uu_112;
fp partial_d_g_uu_113;
fp partial_d_g_uu_122;
fp partial_d_g_uu_123;
fp partial_d_g_uu_133;
fp partial_d_g_uu_211;
fp partial_d_g_uu_212;
fp partial_d_g_uu_213;
fp partial_d_g_uu_222;
fp partial_d_g_uu_223;
fp partial_d_g_uu_233;
fp partial_d_g_uu_311;
fp partial_d_g_uu_312;
fp partial_d_g_uu_313;
fp partial_d_g_uu_322;
fp partial_d_g_uu_323;
fp partial_d_g_uu_333;

fp Theta_A;
fp Theta_B;
fp Theta_C;
fp Theta_D;

#endif
