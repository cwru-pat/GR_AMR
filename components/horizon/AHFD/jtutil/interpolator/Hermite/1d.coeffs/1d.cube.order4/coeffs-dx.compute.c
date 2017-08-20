fp t30;
fp t29;
fp t28;
fp t27;
      t30 = x*x;
      t29 = x*t30;
      t28 = t30*t30;
      t27 = RATIONAL(-1.0,8.0)*t30;
      coeffs_dx->coeff_m2 = RATIONAL(-5.0,24.0)*t28+RATIONAL(1.0,2.0)*t29+
RATIONAL(1.0,12.0)+t27+RATIONAL(-1.0,4.0)*x;
      coeffs_dx->coeff_m1 = RATIONAL(13.0,6.0)*x+RATIONAL(25.0,24.0)*t28+t27+
RATIONAL(-2.0,3.0)+RATIONAL(-7.0,3.0)*t29;
      coeffs_dx->coeff_0 = RATIONAL(5.0,4.0)*t30+RATIONAL(13.0,3.0)*t29+
RATIONAL(-25.0,6.0)*x+RATIONAL(-25.0,12.0)*t28;
      coeffs_dx->coeff_p1 = RATIONAL(2.0,3.0)+RATIONAL(25.0,12.0)*t28+RATIONAL(
-7.0,4.0)*t30+RATIONAL(-4.0,1.0)*t29+RATIONAL(3.0,1.0)*x;
      coeffs_dx->coeff_p2 = RATIONAL(-25.0,24.0)*t28+RATIONAL(-1.0,12.0)+
RATIONAL(11.0,6.0)*t29+RATIONAL(-11.0,12.0)*x+RATIONAL(7.0,8.0)*t30;
      coeffs_dx->coeff_p3 = RATIONAL(-1.0,3.0)*t29+RATIONAL(1.0,6.0)*x+RATIONAL
(5.0,24.0)*t28+t27;
