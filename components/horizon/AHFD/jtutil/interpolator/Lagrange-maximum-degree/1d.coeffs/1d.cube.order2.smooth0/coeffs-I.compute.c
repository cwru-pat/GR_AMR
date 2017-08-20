fp t3;
fp t2;
fp t1;
      t3 = x*x;
      t2 = RATIONAL(1.0,2.0);
      t1 = t2*t3;
      coeffs_I->coeff_m1 = RATIONAL(-1.0,2.0)*x+t1;
      coeffs_I->coeff_0 = RATIONAL(1.0,1.0)-t3;
      coeffs_I->coeff_p1 = t2*x+t1;
