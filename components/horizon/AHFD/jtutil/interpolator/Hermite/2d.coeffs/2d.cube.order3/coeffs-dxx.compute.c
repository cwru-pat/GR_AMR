fp t875;
fp t840;
fp t874;
fp t851;
fp t852;
fp t873;
fp t836;
fp t872;
fp t849;
fp t871;
fp t870;
fp t869;
fp t848;
fp t830;
fp t795;
fp t868;
fp t829;
fp t785;
fp t825;
fp t867;
fp t810;
fp t866;
fp t819;
fp t796;
fp t865;
fp t828;
fp t799;
fp t864;
fp t788;
fp t863;
fp t862;
fp t861;
fp t860;
fp t821;
fp t859;
fp t858;
fp t850;
fp t827;
fp t789;
fp t857;
fp t845;
fp t856;
fp t855;
fp t820;
fp t854;
fp t847;
fp t846;
fp t844;
fp t843;
fp t842;
fp t841;
fp t839;
fp t838;
fp t837;
fp t835;
fp t834;
fp t832;
fp t831;
fp t826;
fp t824;
fp t823;
fp t822;
fp t803;
fp t802;
fp t801;
fp t800;
fp t798;
fp t797;
fp t794;
fp t793;
fp t791;
fp t787;
fp t786;
fp t784;
fp t783;
fp t782;
fp t780;
      t875 = x*y;
      t840 = RATIONAL(2.0,3.0);
      t874 = t840*y;
      t851 = y*y;
      t852 = t851*y;
      t873 = y+t852;
      t836 = RATIONAL(-1.0,24.0);
      t872 = t836*x;
      t849 = RATIONAL(-2.0,3.0);
      t871 = t849*y;
      t870 = x*t852;
      t869 = x*t851;
      t848 = RATIONAL(-5.0,3.0);
      t830 = RATIONAL(-14.0,3.0);
      t795 = t830*t870;
      t868 = t848*t851+t795;
      t829 = RATIONAL(-7.0,24.0);
      t785 = t829*t870;
      t825 = RATIONAL(-1.0,12.0);
      t867 = t785+t825*t851;
      t810 = RATIONAL(-35.0,6.0)*t851;
      t866 = t810+t795;
      t819 = RATIONAL(5.0,18.0);
      t796 = t840*t870;
      t865 = t819*t851+t796;
      t828 = RATIONAL(-7.0,18.0);
      t799 = t849*t870;
      t864 = t828*t851+t799;
      t788 = RATIONAL(14.0,3.0)*t870;
      t863 = t788+RATIONAL(25.0,6.0)*t851;
      t862 = RATIONAL(-5.0,12.0)*t851+t785;
      t861 = RATIONAL(-32.0,3.0)*t870+RATIONAL(-70.0,9.0)*t851;
      t860 = RATIONAL(-5.0,9.0)*t851+t799;
      t821 = RATIONAL(-1.0,36.0);
      t859 = t836*t870+t821*t851;
      t858 = RATIONAL(7.0,9.0)*t851+t796;
      t850 = RATIONAL(1.0,6.0);
      t827 = RATIONAL(7.0,24.0);
      t789 = t827*t870;
      t857 = t850*t851+t789;
      t845 = RATIONAL(7.0,3.0);
      t856 = t788+t845*t851;
      t855 = RATIONAL(-49.0,24.0)*t870+RATIONAL(-5.0,4.0)*t851;
      t820 = RATIONAL(5.0,24.0);
      t854 = t789+t820*t851;
      t847 = RATIONAL(-2.0,9.0);
      t846 = RATIONAL(1.0,2.0);
      t844 = RATIONAL(-1.0,3.0);
      t843 = RATIONAL(2.0,9.0);
      t842 = RATIONAL(4.0,3.0);
      t841 = RATIONAL(-4.0,3.0);
      t839 = RATIONAL(1.0,72.0);
      t838 = RATIONAL(-5.0,24.0);
      t837 = RATIONAL(-7.0,12.0);
      t835 = RATIONAL(1.0,12.0);
      t834 = RATIONAL(7.0,18.0);
      t832 = RATIONAL(1.0,24.0);
      t831 = RATIONAL(10.0,3.0);
      t826 = RATIONAL(-5.0,18.0);
      t824 = RATIONAL(7.0,12.0);
      t823 = RATIONAL(-1.0,72.0);
      t822 = RATIONAL(1.0,36.0);
      t803 = RATIONAL(-7.0,3.0)*t875;
      t802 = RATIONAL(1.0,3.0)*t875;
      t801 = t845*t875;
      t800 = t844*t875;
      t798 = x*t871;
      t797 = x*t874;
      t794 = RATIONAL(32.0,3.0)*t870;
      t793 = RATIONAL(-16.0,3.0)*t875;
      t791 = y*t872;
      t787 = t827*t875;
      t786 = RATIONAL(16.0,3.0)*t875;
      t784 = t829*t875;
      t783 = t832*t870;
      t782 = t832*t875;
      t780 = RATIONAL(49.0,24.0)*t870;
      coeffs_dxx->coeff_m2_m2 = t782+t783+t873*t821+(t825*x+RATIONAL(1.0,18.0))
*t851;
      coeffs_dxx->coeff_m1_m2 = t784+t824*t869+t873*t820+t862;
      coeffs_dxx->coeff_0_m2 = t797+t841*t869+t873*t828+t858;
      coeffs_dxx->coeff_p1_m2 = t798+t842*t869+t873*t819+t860;
      coeffs_dxx->coeff_p2_m2 = t837*t869+t787+t873*t825+t857;
      coeffs_dxx->coeff_p3_m2 = t835*t869+t791+t873*t839+t859;
      coeffs_dxx->coeff_m2_m1 = t843*y+RATIONAL(7.0,36.0)*t852+RATIONAL(5.0,8.0
)*t869+t800+t862;
      coeffs_dxx->coeff_m1_m1 = RATIONAL(-35.0,24.0)*t852+t780+t801+t848*y+(
RATIONAL(25.0,8.0)+RATIONAL(-35.0,8.0)*x)*t851;
      coeffs_dxx->coeff_0_m1 = RATIONAL(28.0,9.0)*y+t793+RATIONAL(49.0,18.0)*
t852+RATIONAL(10.0,1.0)*t869+t866;
      coeffs_dxx->coeff_p1_m1 = t786+RATIONAL(-10.0,1.0)*t869+RATIONAL(-20.0,
9.0)*y+RATIONAL(-35.0,18.0)*t852+t863;
      coeffs_dxx->coeff_p2_m1 = RATIONAL(35.0,8.0)*t869+t874+t803+t824*t852+
t855;
      coeffs_dxx->coeff_p3_m1 = RATIONAL(-7.0,72.0)*t852+RATIONAL(-1.0,9.0)*y+
t802+RATIONAL(-5.0,8.0)*t869+t854;
      coeffs_dxx->coeff_m2_0 = RATIONAL(-4.0,9.0)*t852+t844+(t846+RATIONAL(-7.0
,6.0)*t851)*x+t858;
      coeffs_dxx->coeff_m1_0 = t831*t852+RATIONAL(5.0,2.0)+(RATIONAL(-7.0,2.0)+
RATIONAL(49.0,6.0)*t851)*x+t866;
      coeffs_dxx->coeff_0_0 = RATIONAL(-56.0,9.0)*t852+t830+RATIONAL(98.0,9.0)*
t851+t794+(RATIONAL(8.0,1.0)+RATIONAL(-56.0,3.0)*t851)*x;
      coeffs_dxx->coeff_p1_0 = t831+RATIONAL(40.0,9.0)*t852+(RATIONAL(56.0,3.0)
*t851+RATIONAL(-8.0,1.0))*x+t861;
      coeffs_dxx->coeff_p2_0 = t841*t852+RATIONAL(-1.0,1.0)+(RATIONAL(7.0,2.0)+
RATIONAL(-49.0,6.0)*t851)*x+t856;
      coeffs_dxx->coeff_p3_0 = t843*t852+t850+(RATIONAL(-1.0,2.0)+RATIONAL(7.0,
6.0)*t851)*x+t864;
      coeffs_dxx->coeff_m2_p1 = RATIONAL(5.0,6.0)*t869+RATIONAL(4.0,9.0)*t852+
t802+t847*y+t860;
      coeffs_dxx->coeff_m1_p1 = RATIONAL(5.0,3.0)*y+RATIONAL(-10.0,3.0)*t852+
t803+x*t810+t863;
      coeffs_dxx->coeff_0_p1 = RATIONAL(56.0,9.0)*t852+RATIONAL(-28.0,9.0)*y+
RATIONAL(40.0,3.0)*t869+t786+t861;
      coeffs_dxx->coeff_p1_p1 = t793+t794+RATIONAL(20.0,9.0)*y+RATIONAL(-40.0,
9.0)*t852+(RATIONAL(-40.0,3.0)*x+RATIONAL(50.0,9.0))*t851;
      coeffs_dxx->coeff_p2_p1 = t801+RATIONAL(35.0,6.0)*t869+t842*t852+t871+
t868;
      coeffs_dxx->coeff_p3_p1 = t800+t847*t852+RATIONAL(1.0,9.0)*y+RATIONAL(
-5.0,6.0)*t869+t865;
      coeffs_dxx->coeff_m2_p2 = RATIONAL(-7.0,36.0)*t852+t791+RATIONAL(-1.0,4.0
)*t869+t822*y+t857;
      coeffs_dxx->coeff_m1_p2 = t838*y+RATIONAL(35.0,24.0)*t852+RATIONAL(7.0,
4.0)*t869+t787+t855;
      coeffs_dxx->coeff_0_p2 = RATIONAL(-49.0,18.0)*t852+t834*y+t798+RATIONAL(
-4.0,1.0)*t869+t856;
      coeffs_dxx->coeff_p1_p2 = t797+RATIONAL(4.0,1.0)*t869+t826*y+RATIONAL(
35.0,18.0)*t852+t868;
      coeffs_dxx->coeff_p2_p2 = t837*t852+t835*y+t784+t780+(t846+RATIONAL(-7.0,
4.0)*x)*t851;
      coeffs_dxx->coeff_p3_p2 = RATIONAL(7.0,72.0)*t852+t823*y+t782+RATIONAL(
1.0,4.0)*t869+t867;
      coeffs_dxx->coeff_m2_p3 = t832*t869+t822*t852+t859;
      coeffs_dxx->coeff_m1_p3 = t838*t852+t829*t869+t854;
      coeffs_dxx->coeff_0_p3 = t834*t852+t840*t869+t864;
      coeffs_dxx->coeff_p1_p3 = t849*t869+t826*t852+t865;
      coeffs_dxx->coeff_p2_p3 = t827*t869+t835*t852+t867;
      coeffs_dxx->coeff_p3_p3 = t823*t852+t783+(t839+t872)*t851;
