fp t36;
fp t42;
fp t41;
fp t40;
fp t39;
fp t37;
      t36 = RATIONAL(1.0,3.0)*x;
      t42 = RATIONAL(1.0,4.0)*y+t36;
      t41 = t36+RATIONAL(-1.0,4.0)*y;
      t40 = RATIONAL(1.0,6.0);
      t39 = RATIONAL(-1.0,6.0);
      t37 = RATIONAL(-2.0,3.0)*x;
      coeffs_dx->coeff_m1_m1 = t39+t42;
      coeffs_dx->coeff_0_m1 = t37;
      coeffs_dx->coeff_p1_m1 = t40+t41;
      coeffs_dx->coeff_m1_0 = t36+t39;
      coeffs_dx->coeff_0_0 = t37;
      coeffs_dx->coeff_p1_0 = t36+t40;
      coeffs_dx->coeff_m1_p1 = t39+t41;
      coeffs_dx->coeff_0_p1 = t37;
      coeffs_dx->coeff_p1_p1 = t40+t42;