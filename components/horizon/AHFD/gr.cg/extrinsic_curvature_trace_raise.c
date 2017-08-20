/*
 * inputs = {K_dd, g_uu}
 * outputs = {K_uu, K}
 * cost = 37*assignments+44*additions+84*multiplications
 */
fp t1, t2, t4, t5, t8, t9, t12, t13, t15, t16;
fp t19, t20, t22, t24, t27, t30, t32, t35, t42, t44;
fp t46, t48, t50, t60, t62, t69, t71, t74, t85, t95;
      t1 = g_uu_11;
      t2 = K_dd_11;
      t4 = g_uu_12;
      t5 = K_dd_12;
      t8 = g_uu_13;
      t9 = K_dd_13;
      t12 = g_uu_22;
      t13 = K_dd_22;
      t15 = g_uu_23;
      t16 = K_dd_23;
      t19 = g_uu_33;
      t20 = K_dd_33;
      K = t1*t2+2.0*t4*t5+2.0*t8*t9+t12*t13+2.0*t15*t16+t19*t20;
      t22 = t1*t1;
      t24 = t4*t1;
      t27 = t8*t1;
      t30 = t4*t4;
      t32 = t8*t4;
      t35 = t8*t8;
      K_uu_11 = t22*t2+2.0*t24*t5+2.0*t27*t9+t30*t13+2.0*t32*t16+t35*t20;
      t42 = t4*t12;
      t44 = t8*t12;
      t46 = t1*t15;
      t48 = t15*t4;
      t50 = t8*t15;
      K_uu_12 = t24*t2+t30*t5+t32*t9+t1*t12*t5+t42*t13+t44*t16+t46*t9+t48*t16+
t50*t20;
      t60 = t4*t19;
      t62 = t8*t19;
      K_uu_13 = t27*t2+t32*t5+t35*t9+t46*t5+t48*t13+t50*t16+t1*t19*t9+t60*t16+
t62*t20;
      t69 = t12*t12;
      t71 = t15*t12;
      t74 = t15*t15;
      K_uu_22 = t30*t2+2.0*t42*t5+2.0*t48*t9+t69*t13+2.0*t71*t16+t74*t20;
      t85 = t15*t19;
      K_uu_23 = t32*t2+t44*t5+t50*t9+t48*t5+t71*t13+t74*t16+t60*t9+t12*t19*t16+
t85*t20;
      t95 = t19*t19;
      K_uu_33 = t35*t2+2.0*t50*t5+2.0*t62*t9+t74*t13+2.0*t85*t16+t95*t20;
