fp t25;
fp t23;
fp t29;
fp t21;
fp t22;
fp t28;
fp t24;
fp t20;
fp t19;
      t25 = x*x;
      t23 = t25*t25;
      t29 = RATIONAL(1.0,24.0)*t23+RATIONAL(-1.0,24.0)*t25;
      t21 = RATIONAL(-1.0,6.0);
      t22 = RATIONAL(2.0,3.0);
      t28 = t21*t23+t22*t25;
      t24 = x*t25;
      t20 = RATIONAL(1.0,12.0);
      t19 = RATIONAL(-1.0,12.0);
      coeffs_I->coeff_m2 = t20*x+t19*t24+t29;
      coeffs_I->coeff_m1 = RATIONAL(-2.0,3.0)*x+RATIONAL(1.0,6.0)*t24+t28;
      coeffs_I->coeff_0 = RATIONAL(-5.0,4.0)*t25+RATIONAL(1.0,1.0)+RATIONAL(1.0
,4.0)*t23;
      coeffs_I->coeff_p1 = t21*t24+t22*x+t28;
      coeffs_I->coeff_p2 = t20*t24+t19*x+t29;
