/*
 * inputs = {partial_d_ln_sqrt_g, partial_d_g_uu, r, X_ud, X_udd, g_uu, K_uu, h}
 * outputs = {Theta_A, Theta_B, Theta_C, Theta_D}
 * cost = 134*assignments+3*divisions+5*functions+401*multiplications+173*additions
 */
fp t1, t2, t4, t5, t7, t8, t10, t11, t12, t14;
fp t16, t18, t19, t21, t23, t29, t31, t33, t35, t40;
fp t42, t46, t47, t48, t49, t56, t61, t63, t67, t83;
fp t87, t102, t106, t107, t108, t109, t113, t116, t123, t128;
fp t130, t134, t135, t138, t139, t140, t141, t144, t145, t148;
fp t149, t155, t156, t159, t160, t161, t162, t163, t166, t167;
fp t168, t171, t172, t173, t176, t177, t180, t181, t182, t185;
fp t186, t197, t198, t205, t207, t210, t213, t221, t226, t228;
fp t230, t241, t242, t248, t255, t258, t259, t260, t263, t270;
fp t271, t274, t277, t278, t281, t282, t285, t291, t293, t294;
fp t295, t297, t298, t301, t304, t309, t318, t323, t324, t330;
fp t333, t340, t341, t344, t363, t369, t376, t378, t380, t384;
fp t393, t396, t398, t411, t431, t440, t444, t447, t450, t465;
      t1 = g_uu_13;
      t2 = 1/r;
      t4 = X_ud_13;
      t5 = PARTIAL_RHO(h);
      t7 = X_ud_23;
      t8 = PARTIAL_SIGMA(h);
      t10 = zz*t2-t4*t5-t7*t8;
      t11 = t1*t10;
      t12 = g_uu_23;
      t14 = X_ud_12;
      t16 = X_ud_22;
      t18 = yy*t2-t14*t5-t16*t8;
      t19 = t12*t18;
      t21 = r*r;
      t23 = 1/t21/r;
      t29 = X_ud_11;
      t31 = PARTIAL_RHO_RHO(h);
      t33 = X_ud_21;
      t35 = PARTIAL_RHO_SIGMA(h);
      t40 = PARTIAL_SIGMA_SIGMA(h);
      t42 = -xx*zz*t23-X_udd_113*t5-X_udd_213*t8-t29*t4*t31-t33*t4*t35-t29*t7*
t35-t33*t7*t40;
      t46 = g_uu_12;
      t47 = t46*t18;
      t48 = yy*yy;
      t49 = zz*zz;
      t56 = t29*t29;
      t61 = t33*t33;
      t63 = (t48+t49)*t23-X_udd_111*t5-X_udd_211*t8-t56*t31-2.0*t33*t29*t35-t61
*t40;
      t67 = t12*t12;
      t83 = -yy*zz*t23-X_udd_123*t5-X_udd_223*t8-t14*t4*t31-t16*t4*t35-t14*t7*
t35-t16*t7*t40;
      t87 = t12*t10;
      t102 = -xx*yy*t23-X_udd_112*t5-X_udd_212*t8-t29*t14*t31-t33*t14*t35-t29*
t16*t35-t33*t16*t40;
      t106 = g_uu_33;
      t107 = t106*t10;
      t108 = partial_d_g_uu_312;
      t109 = t108*t18;
      t113 = xx*t2-t29*t5-t33*t8;
      t116 = xx*xx;
      t123 = t4*t4;
      t128 = t7*t7;
      t130 = (t116+t48)*t23-X_udd_133*t5-X_udd_233*t8-t123*t31-2.0*t7*t4*t35-
t128*t40;
      t134 = partial_d_g_uu_212;
      t135 = t134*t18;
      t138 = g_uu_22;
      t139 = t138*t18;
      t140 = partial_d_g_uu_213;
      t141 = t140*t10;
      t144 = partial_d_g_uu_113;
      t145 = t144*t10;
      t148 = partial_d_g_uu_313;
      t149 = t148*t10;
      t155 = t10*t10;
      t156 = t12*t155;
      t159 = -2.0*t11*t19*t42-2.0*t11*t47*t63-2.0*t67*t10*t18*t83-2.0*t87*t47*
t102-t107*t109*t113-2.0*t107*t19*t130-t87*t135*t113-t139*t141*t113-t47*t145*
t113-t19*t149*t113-2.0*t107*t139*t83-t156*t140*t113;
      t160 = g_uu_11;
      t161 = t160*t113;
      t162 = partial_d_g_uu_133;
      t163 = t162*t155;
      t166 = t46*t113;
      t167 = partial_d_g_uu_233;
      t168 = t167*t155;
      t171 = partial_d_g_uu_222;
      t172 = t18*t18;
      t173 = t171*t172;
      t176 = partial_d_g_uu_122;
      t177 = t176*t172;
      t180 = t1*t113;
      t181 = partial_d_g_uu_333;
      t182 = t181*t155;
      t185 = partial_d_g_uu_322;
      t186 = t185*t172;
      t197 = t113*t113;
      t198 = t46*t197;
      t205 = RATIONAL(-1.0,2.0)*t161*t163+RATIONAL(-1.0,2.0)*t166*t168+RATIONAL
(-1.0,2.0)*t166*t173+RATIONAL(-1.0,2.0)*t161*t177+RATIONAL(-1.0,2.0)*t180*t182+
RATIONAL(-1.0,2.0)*t180*t186+RATIONAL(-1.0,2.0)*t47*t163+RATIONAL(-1.0,2.0)*
t139*t168+RATIONAL(-1.0,2.0)*t19*t182+RATIONAL(-1.0,2.0)*t11*t177-2.0*t198*t160
*t102-2.0*t11*t139*t102;
      t207 = t160*t160;
      t210 = t106*t106;
      t213 = t46*t46;
      t221 = t14*t14;
      t226 = t16*t16;
      t228 = (t116+t49)*t23-X_udd_122*t5-X_udd_222*t8-t221*t31-2.0*t16*t14*t35-
t226*t40;
      t230 = t1*t197;
      t241 = partial_d_g_uu_112;
      t242 = t241*t18;
      t248 = t172*t18;
      t255 = t180*t130;
      t258 = -t207*t197*t63-t210*t155*t130-t213*t197*t228-2.0*t230*t160*t42-2.0
*t230*t46*t83+RATIONAL(-1.0,2.0)*t107*t186+RATIONAL(-1.0,2.0)*t87*t173-t11*t242
*t113-2.0*t47*t180*t42+RATIONAL(-1.0,2.0)*t12*t248*t185-2.0*t87*t139*t228-2.0*
t19*t255;
      t259 = t106*t155;
      t260 = t12*t83;
      t263 = t180*t83;
      t270 = partial_d_g_uu_123;
      t271 = t270*t10;
      t274 = t166*t228;
      t277 = partial_d_g_uu_211;
      t278 = t277*t197;
      t281 = partial_d_g_uu_111;
      t282 = t281*t197;
      t285 = t161*t102;
      t291 = t46*t172;
      t293 = t138*t172;
      t294 = partial_d_g_uu_223;
      t295 = t294*t10;
      t297 = partial_d_g_uu_311;
      t298 = t297*t197;
      t301 = t161*t42;
      t304 = -2.0*t259*t260-2.0*t139*t263-2.0*t213*t18*t113*t102-t161*t271*t18
-2.0*t139*t274+RATIONAL(-1.0,2.0)*t139*t278+RATIONAL(-1.0,2.0)*t47*t282-2.0*
t139*t285+RATIONAL(-1.0,2.0)*t138*t248*t171-t291*t271-t293*t295+RATIONAL(-1.0,
2.0)*t107*t298-2.0*t19*t301;
      t309 = t1*t1;
      t318 = t166*t83;
      t323 = partial_d_g_uu_323;
      t324 = t323*t10;
      t330 = t161*t63;
      t333 = t155*t10;
      t340 = RATIONAL(-1.0,2.0)*t87*t278-t309*t197*t130-t213*t172*t63+RATIONAL(
-1.0,2.0)*t11*t282+RATIONAL(-1.0,2.0)*t19*t298-2.0*t19*t318-t166*t295*t18-t180*
t324*t18+RATIONAL(-1.0,2.0)*t46*t248*t176-2.0*t47*t330+RATIONAL(-1.0,2.0)*t12*
t333*t167+RATIONAL(-1.0,2.0)*t1*t333*t162;
      t341 = t46*t102;
      t344 = t12*t172;
      t363 = t138*t138;
      t369 = t1*t42;
      t376 = -2.0*t293*t341-2.0*t344*t46*t42-2.0*t344*t138*t83+RATIONAL(-1.0,
2.0)*t106*t333*t181-t344*t324-t67*t172*t130-t309*t155*t63-t67*t155*t228-t156*
t294*t18-t363*t172*t228-2.0*t156*t1*t102-2.0*t259*t369-2.0*t309*t10*t113*t42;
      t378 = t323*t18;
      t380 = t160*t197;
      t384 = t134*t113;
      t393 = t1*t155;
      t396 = t148*t113;
      t398 = -t259*t378-t380*t242-t198*t135-t230*t109-t293*t384-t291*t241*t113-
t344*t108*t113-t198*t141-t380*t145-t230*t149-t393*t144*t113-t259*t396;
      t411 = t197*t113;
      t431 = -t393*t270*t18-2.0*t107*t47*t42-2.0*t11*t166*t102-2.0*t87*t274-2.0
*t107*t255+RATIONAL(-1.0,2.0)*t160*t411*t281+RATIONAL(-1.0,2.0)*t46*t411*t277+
RATIONAL(-1.0,2.0)*t1*t411*t297-2.0*t87*t263-2.0*t87*t285-2.0*t11*t330-2.0*t107
*t318-2.0*t107*t301;
      Theta_A = t159+t205+t258+t304+t340+t376+t398+t431;
      t440 = t281*t113+t384+t396+t242+t171*t18+t378+t145+t295+t181*t10+t160*t63
+2.0*t341+2.0*t369;
      t444 = partial_d_ln_sqrt_g_1;
      t447 = partial_d_ln_sqrt_g_2;
      t450 = partial_d_ln_sqrt_g_3;
      t465 = t138*t228+2.0*t260+t106*t130+t160*t444*t113+t46*t447*t113+t1*t450*
t113+t46*t444*t18+t138*t447*t18+t12*t450*t18+t1*t444*t10+t12*t447*t10+t106*t450
*t10;
      Theta_B = t440+t465;
      Theta_C = K_uu_11*t197+2.0*K_uu_12*t18*t113+2.0*K_uu_13*t10*t113+K_uu_22*
t172+2.0*K_uu_23*t10*t18+K_uu_33*t155;
      Theta_D = t380+2.0*t47*t113+2.0*t11*t113+t293+2.0*t87*t18+t259;
