fp t103;
fp t106;
fp t107;
fp t105;
fp t113;
fp t112;
fp t100;
fp t96;
fp t101;
fp t111;
fp t104;
fp t99;
fp t98;
fp t97;
      t103 = x*x;
      t106 = t103*t103;
      t107 = x*t103;
      t105 = t107*t107;
      t113 = RATIONAL(1.0,720.0)*t105+RATIONAL(1.0,180.0)*t103+RATIONAL(-1.0,
144.0)*t106;
      t112 = RATIONAL(-1.0,120.0)*t105+RATIONAL(-3.0,40.0)*t103+RATIONAL(1.0,
12.0)*t106;
      t100 = RATIONAL(1.0,48.0);
      t96 = RATIONAL(-13.0,48.0);
      t101 = RATIONAL(3.0,4.0);
      t111 = t100*t105+t96*t106+t101*t103;
      t104 = t107*t103;
      t99 = RATIONAL(-1.0,60.0);
      t98 = RATIONAL(1.0,60.0);
      t97 = RATIONAL(-1.0,48.0);
      coeffs_I->coeff_m3 = t99*x+t100*t107+RATIONAL(-1.0,240.0)*t104+t113;
      coeffs_I->coeff_m2 = RATIONAL(3.0,20.0)*x+RATIONAL(-1.0,6.0)*t107+t98*
t104+t112;
      coeffs_I->coeff_m1 = RATIONAL(-3.0,4.0)*x+RATIONAL(13.0,48.0)*t107+t97*
t104+t111;
      coeffs_I->coeff_0 = RATIONAL(7.0,18.0)*t106+RATIONAL(-49.0,36.0)*t103+
RATIONAL(1.0,1.0)+RATIONAL(-1.0,36.0)*t105;
      coeffs_I->coeff_p1 = t96*t107+t101*x+t100*t104+t111;
      coeffs_I->coeff_p2 = RATIONAL(1.0,6.0)*t107+RATIONAL(-3.0,20.0)*x+t99*
t104+t112;
      coeffs_I->coeff_p3 = t97*t107+t98*x+RATIONAL(1.0,240.0)*t104+t113;
