fp t1, t2, t3, t5, t6, t7, t9, t10, t12, t28;
fp t30, t33, t35, t36, t40, t42, t43, t48, t49, t52;
fp t55;
      t1 = a*a;
      t2 = b*b;
      t3 = t1*t2;
      t5 = t3*zcos*CW;
      t6 = c*c;
      t7 = t1*t6;
      t9 = t7*ycos*BV;
      t10 = t2*t6;
      t12 = t10*xcos*AU;
      t28 = xcos*xcos;
      t30 = CW*CW;
      t33 = BV*BV;
      t35 = t10*t28;
      t36 = ycos*ycos;
      t40 = AU*AU;
      t42 = t7*t36;
      t43 = zcos*zcos;
      t48 = t3*t43;
      t49 = -2.0*t1*zcos*CW*ycos*BV-2.0*t2*zcos*CW*xcos*AU-2.0*t6*ycos*BV*xcos*
AU+t2*t28*t30+t6*t28*t33-t35+t1*t36*t30+t6*t36*t40-t42+t1*t43*t33+t2*t43*t40-
t48;
      t52 = sqrt(-t3*t6*t49);
      t55 = 1/(t35+t42+t48);
      r_plus = (t5+t9+t12+t52)*t55;
      r_minus = (t5+t9+t12-t52)*t55;
