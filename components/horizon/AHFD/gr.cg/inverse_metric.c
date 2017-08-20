/*
 * inputs = {g_dd}
 * outputs = {g_uu}
 * cost = 18*assignments+24*multiplications+divisions+13*additions
 */
fp t1, t2, t4, t5, t7, t8, t11, t12, t14, t15;
fp t18, t21;
      t1 = g_dd_22;
      t2 = g_dd_33;
      t4 = g_dd_23;
      t5 = t4*t4;
      t7 = g_dd_11;
      t8 = t7*t1;
      t11 = g_dd_12;
      t12 = t11*t11;
      t14 = g_dd_13;
      t15 = t11*t14;
      t18 = t14*t14;
      t21 = 1/(t8*t2-t7*t5-t12*t2+2.0*t15*t4-t18*t1);
      g_uu_11 = (t1*t2-t5)*t21;
      g_uu_12 = -(t11*t2-t14*t4)*t21;
      g_uu_13 = -(-t11*t4+t14*t1)*t21;
      g_uu_22 = (t7*t2-t18)*t21;
      g_uu_23 = -(t7*t4-t15)*t21;
      g_uu_33 = (t8-t12)*t21;
