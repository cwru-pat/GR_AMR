fp t14;
      t14 = RATIONAL(1.0,1.0);
      coeffs_dxx->coeff_m1 = t14-x;
      coeffs_dxx->coeff_0 = RATIONAL(-2.0,1.0)+RATIONAL(3.0,1.0)*x;
      coeffs_dxx->coeff_p1 = t14+RATIONAL(-3.0,1.0)*x;
      coeffs_dxx->coeff_p2 = x;
