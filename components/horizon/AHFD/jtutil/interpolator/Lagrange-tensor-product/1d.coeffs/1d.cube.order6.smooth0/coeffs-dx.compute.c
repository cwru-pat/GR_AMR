fp t123;
fp t125;
fp t128;
fp t126;
fp t133;
fp t132;
fp t131;
fp t127;
      t123 = RATIONAL(-3.0,20.0);
      t125 = x*x;
      t128 = x*t125;
      t126 = t128*t125;
      t133 = t123*x+RATIONAL(-1.0,20.0)*t126+RATIONAL(1.0,3.0)*t128;
      t132 = RATIONAL(1.0,120.0)*t126+RATIONAL(1.0,90.0)*x+RATIONAL(-1.0,36.0)*
t128;
      t131 = RATIONAL(3.0,2.0)*x+RATIONAL(-13.0,12.0)*t128+RATIONAL(1.0,8.0)*
t126;
      t127 = t125*t125;
      coeffs_dx->coeff_m3 = RATIONAL(-1.0,48.0)*t127+RATIONAL(-1.0,60.0)+
RATIONAL(1.0,16.0)*t125+t132;
      coeffs_dx->coeff_m2 = RATIONAL(-1.0,2.0)*t125+RATIONAL(3.0,20.0)+RATIONAL
(1.0,12.0)*t127+t133;
      coeffs_dx->coeff_m1 = RATIONAL(-5.0,48.0)*t127+RATIONAL(13.0,16.0)*t125+
RATIONAL(-3.0,4.0)+t131;
      coeffs_dx->coeff_0 = RATIONAL(-1.0,6.0)*t126+RATIONAL(14.0,9.0)*t128+
RATIONAL(-49.0,18.0)*x;
      coeffs_dx->coeff_p1 = RATIONAL(5.0,48.0)*t127+RATIONAL(3.0,4.0)+RATIONAL(
-13.0,16.0)*t125+t131;
      coeffs_dx->coeff_p2 = RATIONAL(1.0,2.0)*t125+RATIONAL(-1.0,12.0)*t127+
t123+t133;
      coeffs_dx->coeff_p3 = RATIONAL(-1.0,16.0)*t125+RATIONAL(1.0,60.0)+
RATIONAL(1.0,48.0)*t127+t132;
