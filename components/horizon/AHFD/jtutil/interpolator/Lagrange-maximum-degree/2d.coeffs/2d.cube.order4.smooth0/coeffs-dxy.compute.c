fp t774;
fp t775;
fp t793;
fp t763;
fp t747;
fp t765;
fp t732;
fp t792;
fp t740;
fp t753;
fp t791;
fp t767;
fp t771;
fp t790;
fp t789;
fp t788;
fp t787;
fp t786;
fp t731;
fp t785;
fp t784;
fp t728;
fp t783;
fp t730;
fp t782;
fp t768;
fp t773;
fp t781;
fp t780;
fp t779;
fp t778;
fp t777;
fp t776;
fp t770;
fp t766;
fp t764;
fp t762;
fp t755;
fp t754;
fp t751;
fp t746;
fp t745;
fp t742;
fp t736;
fp t733;
fp t729;
fp t727;
      t774 = x*x;
      t775 = y*y;
      t793 = t774+t775;
      t763 = RATIONAL(1.0,35.0);
      t747 = t763*x;
      t765 = RATIONAL(-1.0,35.0);
      t732 = t765*y;
      t792 = t747+t732;
      t740 = t765*x;
      t753 = t763*y;
      t791 = t740+t753;
      t767 = RATIONAL(1.0,10.0);
      t771 = RATIONAL(-1.0,40.0);
      t790 = t767*t775+t771*t774;
      t789 = t767*t774+t771*t775;
      t788 = x*y;
      t787 = t793*RATIONAL(-1.0,20.0);
      t786 = t793*RATIONAL(1.0,20.0);
      t731 = RATIONAL(-2.0,49.0)*t788;
      t785 = t731+RATIONAL(21.0,200.0);
      t784 = t731+RATIONAL(-21.0,200.0);
      t728 = RATIONAL(1.0,49.0)*t788;
      t783 = t728+RATIONAL(-37.0,300.0)+t786;
      t730 = RATIONAL(4.0,49.0)*t788;
      t782 = t730+RATIONAL(-11.0,150.0)+t786;
      t768 = RATIONAL(-1.0,10.0);
      t773 = RATIONAL(1.0,40.0);
      t781 = t773*t775+t768*t774+t785;
      t780 = t768*t775+t773*t774+t785;
      t779 = t740+t732+t784;
      t778 = t730+RATIONAL(11.0,150.0)+t787;
      t777 = t747+t753+t784;
      t776 = t728+RATIONAL(37.0,300.0)+t787;
      t770 = RATIONAL(-2.0,35.0);
      t766 = RATIONAL(2.0,35.0);
      t764 = RATIONAL(-1.0,70.0);
      t762 = RATIONAL(1.0,70.0);
      t755 = t762*y;
      t754 = t770*y;
      t751 = t770*x;
      t746 = t766*y;
      t745 = t766*x;
      t742 = t764*y;
      t736 = t762*x;
      t733 = t764*x;
      t729 = RATIONAL(-4.0,49.0)*t788;
      t727 = RATIONAL(2.0,49.0)*t788;
      coeffs_dxy->coeff_m2_m2 = t754+t751+t782;
      coeffs_dxy->coeff_m1_m2 = t781+t792;
      coeffs_dxy->coeff_0_m2 = t729+t745;
      coeffs_dxy->coeff_p1_m2 = t777+t789;
      coeffs_dxy->coeff_p2_m2 = t746+t751+t778;
      coeffs_dxy->coeff_m2_m1 = t780+t791;
      coeffs_dxy->coeff_m1_m1 = t755+t736+t776;
      coeffs_dxy->coeff_0_m1 = t727+t747;
      coeffs_dxy->coeff_p1_m1 = t742+t736+t783;
      coeffs_dxy->coeff_p2_m1 = t779+t790;
      coeffs_dxy->coeff_m2_0 = t746+t729;
      coeffs_dxy->coeff_m1_0 = t727+t753;
      coeffs_dxy->coeff_0_0 = t730;
      coeffs_dxy->coeff_p1_0 = t732+t727;
      coeffs_dxy->coeff_p2_0 = t754+t729;
      coeffs_dxy->coeff_m2_p1 = t777+t790;
      coeffs_dxy->coeff_m1_p1 = t733+t755+t783;
      coeffs_dxy->coeff_0_p1 = t727+t740;
      coeffs_dxy->coeff_p1_p1 = t733+t742+t776;
      coeffs_dxy->coeff_p2_p1 = t780+t792;
      coeffs_dxy->coeff_m2_p2 = t745+t754+t778;
      coeffs_dxy->coeff_m1_p2 = t779+t789;
      coeffs_dxy->coeff_0_p2 = t729+t751;
      coeffs_dxy->coeff_p1_p2 = t781+t791;
      coeffs_dxy->coeff_p2_p2 = t745+t746+t782;
