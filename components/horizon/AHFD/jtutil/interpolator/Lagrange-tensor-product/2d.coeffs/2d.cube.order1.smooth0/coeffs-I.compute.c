fp t1;
fp t2;
      t1 = x*y;
      t2 = y-t1;
      coeffs_I->coeff_0_0 = RATIONAL(1.0,1.0)-x-t2;
      coeffs_I->coeff_p1_0 = x-t1;
      coeffs_I->coeff_0_p1 = t2;
      coeffs_I->coeff_p1_p1 = t1;
