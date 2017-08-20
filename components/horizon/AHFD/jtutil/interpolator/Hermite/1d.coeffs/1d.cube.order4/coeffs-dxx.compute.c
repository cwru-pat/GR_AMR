fp t37;
fp t36;
fp t35;
fp t34;
fp t33;
      t37 = x*x;
      t36 = x*t37;
      t35 = RATIONAL(-1.0,4.0);
      t34 = RATIONAL(-25.0,6.0);
      t33 = t35*x;
      coeffs_dxx->coeff_m2 = t35+t33+RATIONAL(3.0,2.0)*t37+RATIONAL(-5.0,6.0)*
t36;
      coeffs_dxx->coeff_m1 = RATIONAL(13.0,6.0)+t33+RATIONAL(-7.0,1.0)*t37+
RATIONAL(25.0,6.0)*t36;
      coeffs_dxx->coeff_0 = RATIONAL(-25.0,3.0)*t36+RATIONAL(13.0,1.0)*t37+t34+
RATIONAL(5.0,2.0)*x;
      coeffs_dxx->coeff_p1 = RATIONAL(3.0,1.0)+RATIONAL(-7.0,2.0)*x+RATIONAL(
-12.0,1.0)*t37+RATIONAL(25.0,3.0)*t36;
      coeffs_dxx->coeff_p2 = RATIONAL(11.0,2.0)*t37+t34*t36+RATIONAL(7.0,4.0)*x
+RATIONAL(-11.0,12.0);
      coeffs_dxx->coeff_p3 = RATIONAL(5.0,6.0)*t36+RATIONAL(1.0,6.0)+t33-t37;
