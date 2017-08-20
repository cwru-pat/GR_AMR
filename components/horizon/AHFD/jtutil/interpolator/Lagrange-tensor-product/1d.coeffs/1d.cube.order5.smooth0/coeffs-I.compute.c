fp t54;
fp t55;
fp t57;
fp t58;
fp t64;
fp t63;
fp t60;
fp t59;
fp t51;
      t54 = RATIONAL(1.0,24.0);
      t55 = RATIONAL(-1.0,24.0);
      t57 = x*x;
      t58 = t57*t57;
      t64 = t54*t58+t55*t57;
      t63 = RATIONAL(-1.0,6.0)*t58+RATIONAL(2.0,3.0)*t57;
      t60 = x*t57;
      t59 = t60*t57;
      t51 = t55*t60;
      coeffs_I->coeff_m2 = RATIONAL(1.0,20.0)*x+t51+RATIONAL(-1.0,120.0)*t59+
t64;
      coeffs_I->coeff_m1 = RATIONAL(-1.0,2.0)*x+t51+t54*t59+t63;
      coeffs_I->coeff_0 = RATIONAL(1.0,1.0)+RATIONAL(-1.0,3.0)*x+RATIONAL(-5.0,
4.0)*t57+RATIONAL(5.0,12.0)*t60+RATIONAL(1.0,4.0)*t58+RATIONAL(-1.0,12.0)*t59;
      coeffs_I->coeff_p1 = x+RATIONAL(-7.0,12.0)*t60+RATIONAL(1.0,12.0)*t59+t63
;
      coeffs_I->coeff_p2 = RATIONAL(-1.0,4.0)*x+RATIONAL(7.0,24.0)*t60+t55*t59+
t64;
      coeffs_I->coeff_p3 = RATIONAL(1.0,30.0)*x+RATIONAL(1.0,120.0)*t59+t51;
