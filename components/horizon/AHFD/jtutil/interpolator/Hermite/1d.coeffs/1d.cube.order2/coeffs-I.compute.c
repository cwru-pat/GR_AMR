fp t4;
fp t3;
fp t2;
fp t1;
      t4 = x*x;
      t3 = x*t4;
      t2 = RATIONAL(1.0,2.0);
      t1 = RATIONAL(-1.0,2.0);
      coeffs_I->coeff_m1 = t4+(x+t3)*t1;
      coeffs_I->coeff_0 = RATIONAL(-5.0,2.0)*t4+RATIONAL(1.0,1.0)+RATIONAL(3.0,
2.0)*t3;
      coeffs_I->coeff_p1 = RATIONAL(-3.0,2.0)*t3+t2*x+RATIONAL(2.0,1.0)*t4;
      coeffs_I->coeff_p2 = t2*t3+t1*t4;
