fp t13;
fp t12;
      t13 = x*x;
      t12 = RATIONAL(-1.0,2.0);
      coeffs_dx->coeff_m1 = x+t12*t13+RATIONAL(-1.0,3.0);
      coeffs_dx->coeff_0 = t12+RATIONAL(-2.0,1.0)*x+RATIONAL(3.0,2.0)*t13;
      coeffs_dx->coeff_p1 = x+RATIONAL(1.0,1.0)+RATIONAL(-3.0,2.0)*t13;
      coeffs_dx->coeff_p2 = RATIONAL(1.0,2.0)*t13+RATIONAL(-1.0,6.0);
