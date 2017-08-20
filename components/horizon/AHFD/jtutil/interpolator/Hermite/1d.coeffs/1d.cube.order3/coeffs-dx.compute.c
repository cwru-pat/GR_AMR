fp t12;
      t12 = x*x;
      coeffs_dx->coeff_m2 = RATIONAL(-1.0,3.0)*x+RATIONAL(1.0,12.0)+RATIONAL(
1.0,4.0)*t12;
      coeffs_dx->coeff_m1 = RATIONAL(-2.0,3.0)+RATIONAL(5.0,2.0)*x+RATIONAL(
-7.0,4.0)*t12;
      coeffs_dx->coeff_0 = RATIONAL(-14.0,3.0)*x+RATIONAL(4.0,1.0)*t12;
      coeffs_dx->coeff_p1 = RATIONAL(-4.0,1.0)*t12+RATIONAL(10.0,3.0)*x+
RATIONAL(2.0,3.0);
      coeffs_dx->coeff_p2 = RATIONAL(7.0,4.0)*t12+RATIONAL(-1.0,12.0)-x;
      coeffs_dx->coeff_p3 = RATIONAL(1.0,6.0)*x+RATIONAL(-1.0,4.0)*t12;
