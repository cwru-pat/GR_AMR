fp t10;
fp t9;
fp t8;
fp t7;
      t10 = x*x;
      t9 = x*t10;
      t8 = RATIONAL(1.0,12.0);
      t7 = RATIONAL(-1.0,12.0);
      coeffs_I->coeff_m2 = RATIONAL(-1.0,6.0)*t10+(x+t9)*t8;
      coeffs_I->coeff_m1 = RATIONAL(-2.0,3.0)*x+RATIONAL(5.0,4.0)*t10+RATIONAL(
-7.0,12.0)*t9;
      coeffs_I->coeff_0 = RATIONAL(-7.0,3.0)*t10+RATIONAL(1.0,1.0)+RATIONAL(4.0
,3.0)*t9;
      coeffs_I->coeff_p1 = RATIONAL(-4.0,3.0)*t9+RATIONAL(2.0,3.0)*x+RATIONAL(
5.0,3.0)*t10;
      coeffs_I->coeff_p2 = RATIONAL(7.0,12.0)*t9+t7*x+RATIONAL(-1.0,2.0)*t10;
      coeffs_I->coeff_p3 = t7*t9+t8*t10;
