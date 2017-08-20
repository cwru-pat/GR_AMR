fp t144;
fp t145;
fp t150;
fp t149;
fp t141;
fp t148;
fp t143;
      t144 = x*x;
      t145 = t144*t144;
      t150 = RATIONAL(-13.0,4.0)*t144+RATIONAL(3.0,2.0)+RATIONAL(5.0,8.0)*t145;
      t149 = t144+RATIONAL(-3.0,20.0)+RATIONAL(-1.0,4.0)*t145;
      t141 = RATIONAL(-1.0,12.0);
      t148 = RATIONAL(1.0,90.0)+RATIONAL(1.0,24.0)*t145+t141*t144;
      t143 = x*t144;
      coeffs_dxx->coeff_m3 = RATIONAL(1.0,8.0)*x+t141*t143+t148;
      coeffs_dxx->coeff_m2 = -x+RATIONAL(1.0,3.0)*t143+t149;
      coeffs_dxx->coeff_m1 = RATIONAL(-5.0,12.0)*t143+RATIONAL(13.0,8.0)*x+t150
;
      coeffs_dxx->coeff_0 = RATIONAL(14.0,3.0)*t144+RATIONAL(-5.0,6.0)*t145+
RATIONAL(-49.0,18.0);
      coeffs_dxx->coeff_p1 = RATIONAL(5.0,12.0)*t143+RATIONAL(-13.0,8.0)*x+t150
;
      coeffs_dxx->coeff_p2 = x+RATIONAL(-1.0,3.0)*t143+t149;
      coeffs_dxx->coeff_p3 = RATIONAL(-1.0,8.0)*x+RATIONAL(1.0,12.0)*t143+t148;
