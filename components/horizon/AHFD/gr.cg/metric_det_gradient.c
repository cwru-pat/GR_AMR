/*
 * inputs = {g_uu, partial_d_g_dd}
 * outputs = {partial_d_ln_sqrt_g}
 * cost = 9*assignments+15*additions+27*multiplications
 */
fp t1, t5, t8, t11, t15, t18;
      t1 = g_uu_11;
      t5 = g_uu_12;
      t8 = g_uu_13;
      t11 = g_uu_22;
      t15 = g_uu_23;
      t18 = g_uu_33;
      partial_d_ln_sqrt_g_1 = RATIONAL(1.0,2.0)*t1*partial_d_g_dd_111+t5*
partial_d_g_dd_112+t8*partial_d_g_dd_113+RATIONAL(1.0,2.0)*t11*
partial_d_g_dd_122+t15*partial_d_g_dd_123+RATIONAL(1.0,2.0)*t18*
partial_d_g_dd_133;
      partial_d_ln_sqrt_g_2 = RATIONAL(1.0,2.0)*t1*partial_d_g_dd_211+t5*
partial_d_g_dd_212+t8*partial_d_g_dd_213+RATIONAL(1.0,2.0)*t11*
partial_d_g_dd_222+t15*partial_d_g_dd_223+RATIONAL(1.0,2.0)*t18*
partial_d_g_dd_233;
      partial_d_ln_sqrt_g_3 = RATIONAL(1.0,2.0)*t1*partial_d_g_dd_311+t5*
partial_d_g_dd_312+t8*partial_d_g_dd_313+RATIONAL(1.0,2.0)*t11*
partial_d_g_dd_322+t15*partial_d_g_dd_323+RATIONAL(1.0,2.0)*t18*
partial_d_g_dd_333;
