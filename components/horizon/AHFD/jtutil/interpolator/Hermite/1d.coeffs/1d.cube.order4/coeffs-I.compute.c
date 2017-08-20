fp t23;
fp t25;
fp t22;
fp t20;
fp t19;
fp t18;
fp t17;
fp t16;
fp t15;
fp t14;
fp t13;
      t23 = x*x;
      t25 = t23*t23;
      t22 = x*t23;
      t20 = x*t25;
      t19 = RATIONAL(-1.0,24.0);
      t18 = RATIONAL(1.0,12.0);
      t17 = RATIONAL(-7.0,12.0);
      t16 = RATIONAL(5.0,12.0);
      t15 = RATIONAL(-1.0,12.0);
      t14 = RATIONAL(13.0,12.0);
      t13 = t19*t22;
      coeffs_I->coeff_m2 = t18*x+RATIONAL(-1.0,8.0)*t23+t13+RATIONAL(1.0,8.0)*
t25+t19*t20;
      coeffs_I->coeff_m1 = RATIONAL(-2.0,3.0)*x+t14*t23+t13+t17*t25+RATIONAL(
5.0,24.0)*t20;
      coeffs_I->coeff_0 = RATIONAL(-25.0,12.0)*t23+t16*t22+t14*t25+RATIONAL(1.0
,1.0)+RATIONAL(-5.0,12.0)*t20;
      coeffs_I->coeff_p1 = t17*t22+RATIONAL(2.0,3.0)*x+RATIONAL(3.0,2.0)*t23+
t16*t20-t25;
      coeffs_I->coeff_p2 = t15*x+RATIONAL(-11.0,24.0)*t23+RATIONAL(7.0,24.0)*
t22+RATIONAL(-5.0,24.0)*t20+RATIONAL(11.0,24.0)*t25;
      coeffs_I->coeff_p3 = t13+RATIONAL(1.0,24.0)*t20+t18*t23+t15*t25;
