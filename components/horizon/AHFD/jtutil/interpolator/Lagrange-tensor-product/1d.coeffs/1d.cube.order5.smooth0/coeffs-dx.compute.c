fp t71;
fp t70;
fp t76;
fp t75;
fp t72;
fp t67;
      t71 = x*x;
      t70 = t71*x;
      t76 = RATIONAL(-2.0,3.0)*t70+RATIONAL(4.0,3.0)*x;
      t75 = RATIONAL(-1.0,12.0)*x+RATIONAL(1.0,6.0)*t70;
      t72 = t71*t71;
      t67 = RATIONAL(-1.0,8.0)*t71;
      coeffs_dx->coeff_m2 = RATIONAL(-1.0,24.0)*t72+RATIONAL(1.0,20.0)+t67+t75;
      coeffs_dx->coeff_m1 = t67+RATIONAL(-1.0,2.0)+RATIONAL(5.0,24.0)*t72+t76;
      coeffs_dx->coeff_0 = RATIONAL(-5.0,2.0)*x+t70+RATIONAL(-5.0,12.0)*t72+
RATIONAL(-1.0,3.0)+RATIONAL(5.0,4.0)*t71;
      coeffs_dx->coeff_p1 = RATIONAL(5.0,12.0)*t72+RATIONAL(-7.0,4.0)*t71+
RATIONAL(1.0,1.0)+t76;
      coeffs_dx->coeff_p2 = RATIONAL(7.0,8.0)*t71+RATIONAL(-1.0,4.0)+RATIONAL(
-5.0,24.0)*t72+t75;
      coeffs_dx->coeff_p3 = RATIONAL(1.0,24.0)*t72+t67+RATIONAL(1.0,30.0);
