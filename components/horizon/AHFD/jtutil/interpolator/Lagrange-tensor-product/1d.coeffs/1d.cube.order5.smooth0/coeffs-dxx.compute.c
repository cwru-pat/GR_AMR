fp t83;
fp t86;
fp t85;
fp t82;
fp t78;
      t83 = x*x;
      t86 = RATIONAL(1.0,2.0)*t83+RATIONAL(-1.0,12.0);
      t85 = RATIONAL(-2.0,1.0)*t83+RATIONAL(4.0,3.0);
      t82 = x*t83;
      t78 = RATIONAL(-1.0,4.0)*x;
      coeffs_dxx->coeff_m2 = RATIONAL(-1.0,6.0)*t82+t78+t86;
      coeffs_dxx->coeff_m1 = RATIONAL(5.0,6.0)*t82+t78+t85;
      coeffs_dxx->coeff_0 = RATIONAL(-5.0,3.0)*t82+RATIONAL(5.0,2.0)*x+RATIONAL
(3.0,1.0)*t83+RATIONAL(-5.0,2.0);
      coeffs_dxx->coeff_p1 = RATIONAL(5.0,3.0)*t82+RATIONAL(-7.0,2.0)*x+t85;
      coeffs_dxx->coeff_p2 = RATIONAL(-5.0,6.0)*t82+RATIONAL(7.0,4.0)*x+t86;
      coeffs_dxx->coeff_p3 = t78+RATIONAL(1.0,6.0)*t82;
