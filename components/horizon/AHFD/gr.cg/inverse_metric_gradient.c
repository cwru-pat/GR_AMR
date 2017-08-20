/*
 * inputs = {g_uu, partial_d_g_dd}
 * outputs = {partial_d_g_uu}
 * cost = 63*assignments+183*multiplications+117*additions
 */
fp t1, t2, t3, t5, t6, t7, t10, t11, t12, t15;
fp t16, t18, t19, t22, t23, t28, t29, t31, t33, t35;
fp t36, t38, t40, t48, t49, t51, t53, t60, t62, t65;
fp t74, t76, t86, t88, t90, t93, t96, t98, t101, t148;
fp t150, t153, t156, t158, t161;
      t1 = g_uu_11;
      t2 = t1*t1;
      t3 = partial_d_g_dd_111;
      t5 = g_uu_12;
      t6 = t5*t1;
      t7 = partial_d_g_dd_112;
      t10 = g_uu_13;
      t11 = t10*t1;
      t12 = partial_d_g_dd_113;
      t15 = t5*t5;
      t16 = partial_d_g_dd_122;
      t18 = t10*t5;
      t19 = partial_d_g_dd_123;
      t22 = t10*t10;
      t23 = partial_d_g_dd_133;
      partial_d_g_uu_111 = -t2*t3-2.0*t6*t7-2.0*t11*t12-t15*t16-2.0*t18*t19-t22
*t23;
      t28 = g_uu_22;
      t29 = t1*t28;
      t31 = t5*t28;
      t33 = t10*t28;
      t35 = g_uu_23;
      t36 = t1*t35;
      t38 = t5*t35;
      t40 = t10*t35;
      partial_d_g_uu_112 = -t6*t3-t15*t7-t18*t12-t29*t7-t31*t16-t33*t19-t36*t12
-t38*t19-t40*t23;
      t48 = g_uu_33;
      t49 = t1*t48;
      t51 = t48*t5;
      t53 = t10*t48;
      partial_d_g_uu_113 = -t11*t3-t18*t7-t22*t12-t36*t7-t38*t16-t40*t19-t49*
t12-t51*t19-t53*t23;
      t60 = t28*t28;
      t62 = t35*t28;
      t65 = t35*t35;
      partial_d_g_uu_122 = -t15*t3-2.0*t31*t7-2.0*t38*t12-t60*t16-2.0*t62*t19-
t65*t23;
      t74 = t28*t48;
      t76 = t35*t48;
      partial_d_g_uu_123 = -t18*t3-t33*t7-t40*t12-t38*t7-t62*t16-t65*t19-t51*
t12-t74*t19-t76*t23;
      t86 = t48*t48;
      partial_d_g_uu_133 = -t22*t3-2.0*t40*t7-2.0*t53*t12-t65*t16-2.0*t76*t19-
t86*t23;
      t88 = partial_d_g_dd_211;
      t90 = partial_d_g_dd_212;
      t93 = partial_d_g_dd_213;
      t96 = partial_d_g_dd_222;
      t98 = partial_d_g_dd_223;
      t101 = partial_d_g_dd_233;
      partial_d_g_uu_211 = -t2*t88-2.0*t6*t90-2.0*t11*t93-t15*t96-2.0*t18*t98-
t22*t101;
      partial_d_g_uu_212 = -t6*t88-t15*t90-t18*t93-t29*t90-t31*t96-t33*t98-t36*
t93-t38*t98-t40*t101;
      partial_d_g_uu_213 = -t11*t88-t18*t90-t22*t93-t36*t90-t38*t96-t40*t98-t49
*t93-t51*t98-t53*t101;
      partial_d_g_uu_222 = -t15*t88-2.0*t31*t90-2.0*t38*t93-t60*t96-2.0*t62*t98
-t65*t101;
      partial_d_g_uu_223 = -t18*t88-t33*t90-t40*t93-t38*t90-t62*t96-t65*t98-t51
*t93-t74*t98-t76*t101;
      partial_d_g_uu_233 = -t22*t88-2.0*t40*t90-2.0*t53*t93-t65*t96-2.0*t76*t98
-t86*t101;
      t148 = partial_d_g_dd_311;
      t150 = partial_d_g_dd_312;
      t153 = partial_d_g_dd_313;
      t156 = partial_d_g_dd_322;
      t158 = partial_d_g_dd_323;
      t161 = partial_d_g_dd_333;
      partial_d_g_uu_311 = -t2*t148-2.0*t6*t150-2.0*t11*t153-t15*t156-2.0*t18*
t158-t22*t161;
      partial_d_g_uu_312 = -t6*t148-t15*t150-t18*t153-t29*t150-t31*t156-t33*
t158-t36*t153-t38*t158-t40*t161;
      partial_d_g_uu_313 = -t11*t148-t18*t150-t22*t153-t36*t150-t38*t156-t40*
t158-t49*t153-t51*t158-t53*t161;
      partial_d_g_uu_322 = -t15*t148-2.0*t31*t150-2.0*t38*t153-t60*t156-2.0*t62
*t158-t65*t161;
      partial_d_g_uu_323 = -t18*t148-t33*t150-t40*t153-t38*t150-t62*t156-t65*
t158-t51*t153-t74*t158-t76*t161;
      partial_d_g_uu_333 = -t22*t148-2.0*t40*t150-2.0*t53*t153-t65*t156-2.0*t76
*t158-t86*t161;
