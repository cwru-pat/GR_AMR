fp t6;
      t6 = x*x;
      coeffs_dx->coeff_m1 = RATIONAL(2.0,1.0)*x+RATIONAL(-3.0,2.0)*t6+RATIONAL(
-1.0,2.0);
      coeffs_dx->coeff_0 = RATIONAL(-5.0,1.0)*x+RATIONAL(9.0,2.0)*t6;
      coeffs_dx->coeff_p1 = RATIONAL(1.0,2.0)+RATIONAL(4.0,1.0)*x+RATIONAL(-9.0
,2.0)*t6;
      coeffs_dx->coeff_p2 = -x+RATIONAL(3.0,2.0)*t6;
