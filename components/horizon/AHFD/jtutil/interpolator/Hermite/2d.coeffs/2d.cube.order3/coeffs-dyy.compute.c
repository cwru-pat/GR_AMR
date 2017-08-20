fp t1116;
fp t1091;
fp t1115;
fp t1071;
fp t1114;
fp t1070;
fp t1113;
fp t1080;
fp t1112;
fp t1088;
fp t1111;
fp t1075;
fp t1110;
fp t1090;
fp t1109;
fp t1108;
fp t1032;
fp t1107;
fp t1050;
fp t1064;
fp t1028;
fp t1106;
fp t1066;
fp t1105;
fp t1104;
fp t1040;
fp t1058;
fp t1103;
fp t1059;
fp t1068;
fp t1034;
fp t1102;
fp t1101;
fp t1087;
fp t1100;
fp t1027;
fp t1099;
fp t1089;
fp t1098;
fp t1097;
fp t1042;
fp t1096;
fp t1060;
fp t1095;
fp t1082;
fp t1094;
fp t1069;
fp t1093;
fp t1086;
fp t1085;
fp t1084;
fp t1083;
fp t1081;
fp t1079;
fp t1078;
fp t1077;
fp t1076;
fp t1074;
fp t1073;
fp t1067;
fp t1065;
fp t1063;
fp t1062;
fp t1061;
fp t1041;
fp t1039;
fp t1038;
fp t1037;
fp t1036;
fp t1035;
fp t1033;
fp t1031;
fp t1030;
fp t1029;
fp t1025;
fp t1024;
fp t1023;
fp t1022;
fp t1020;
      t1116 = x*y;
      t1091 = x*x;
      t1115 = y*t1091;
      t1071 = RATIONAL(1.0,24.0);
      t1114 = y*t1071;
      t1070 = RATIONAL(-7.0,24.0);
      t1113 = y*t1070;
      t1080 = RATIONAL(2.0,3.0);
      t1112 = y*t1080;
      t1088 = RATIONAL(-2.0,3.0);
      t1111 = y*t1088;
      t1075 = RATIONAL(-1.0,24.0);
      t1110 = t1075*y;
      t1090 = t1091*x;
      t1109 = y*t1090;
      t1108 = x+t1090;
      t1032 = t1070*t1109;
      t1107 = RATIONAL(-5.0,12.0)*t1091+t1032;
      t1050 = RATIONAL(-35.0,6.0)*t1091;
      t1064 = RATIONAL(-14.0,3.0);
      t1028 = t1064*t1109;
      t1106 = t1050+t1028;
      t1066 = RATIONAL(-1.0,12.0);
      t1105 = t1066*t1091+t1032;
      t1104 = RATIONAL(-70.0,9.0)*t1091+RATIONAL(-32.0,3.0)*t1109;
      t1040 = t1080*t1109;
      t1058 = RATIONAL(5.0,18.0);
      t1103 = t1040+t1058*t1091;
      t1059 = RATIONAL(5.0,24.0);
      t1068 = RATIONAL(7.0,24.0);
      t1034 = t1068*t1109;
      t1102 = t1059*t1091+t1034;
      t1101 = RATIONAL(-5.0,4.0)*t1091+RATIONAL(-49.0,24.0)*t1109;
      t1087 = RATIONAL(-5.0,3.0);
      t1100 = t1087*t1091+t1028;
      t1027 = RATIONAL(14.0,3.0)*t1109;
      t1099 = RATIONAL(25.0,6.0)*t1091+t1027;
      t1089 = RATIONAL(1.0,6.0);
      t1098 = t1089*t1091+t1034;
      t1097 = t1040+RATIONAL(7.0,9.0)*t1091;
      t1042 = t1088*t1109;
      t1096 = t1042+RATIONAL(-5.0,9.0)*t1091;
      t1060 = RATIONAL(-1.0,36.0);
      t1095 = t1060*t1091+t1075*t1109;
      t1082 = RATIONAL(7.0,3.0);
      t1094 = t1027+t1082*t1091;
      t1069 = RATIONAL(-7.0,18.0);
      t1093 = t1042+t1069*t1091;
      t1086 = RATIONAL(-2.0,9.0);
      t1085 = RATIONAL(-1.0,3.0);
      t1084 = RATIONAL(2.0,9.0);
      t1083 = RATIONAL(4.0,3.0);
      t1081 = RATIONAL(-4.0,3.0);
      t1079 = RATIONAL(1.0,2.0);
      t1078 = RATIONAL(1.0,72.0);
      t1077 = RATIONAL(-5.0,24.0);
      t1076 = RATIONAL(-7.0,12.0);
      t1074 = RATIONAL(1.0,12.0);
      t1073 = RATIONAL(7.0,18.0);
      t1067 = RATIONAL(-5.0,18.0);
      t1065 = RATIONAL(7.0,12.0);
      t1063 = RATIONAL(10.0,3.0);
      t1062 = RATIONAL(-1.0,72.0);
      t1061 = RATIONAL(1.0,36.0);
      t1041 = RATIONAL(1.0,3.0)*t1116;
      t1039 = t1085*t1116;
      t1038 = t1082*t1116;
      t1037 = RATIONAL(-7.0,3.0)*t1116;
      t1036 = x*t1112;
      t1035 = x*t1111;
      t1033 = x*t1110;
      t1031 = t1071*t1109;
      t1030 = x*t1114;
      t1029 = x*t1113;
      t1025 = t1068*t1116;
      t1024 = RATIONAL(-16.0,3.0)*t1116;
      t1023 = RATIONAL(32.0,3.0)*t1109;
      t1022 = RATIONAL(16.0,3.0)*t1116;
      t1020 = RATIONAL(49.0,24.0)*t1109;
      coeffs_dyy->coeff_m2_m2 = t1030+t1031+(t1066*y+RATIONAL(1.0,18.0))*t1091+
t1108*t1060;
      coeffs_dyy->coeff_m1_m2 = RATIONAL(5.0,8.0)*t1115+RATIONAL(7.0,36.0)*
t1090+t1039+t1084*x+t1107;
      coeffs_dyy->coeff_0_m2 = t1085+RATIONAL(-4.0,9.0)*t1090+(t1079+RATIONAL(
-7.0,6.0)*t1091)*y+t1097;
      coeffs_dyy->coeff_p1_m2 = RATIONAL(5.0,6.0)*t1115+t1041+RATIONAL(4.0,9.0)
*t1090+t1086*x+t1096;
      coeffs_dyy->coeff_p2_m2 = t1033+RATIONAL(-1.0,4.0)*t1115+t1061*x+RATIONAL
(-7.0,36.0)*t1090+t1098;
      coeffs_dyy->coeff_p3_m2 = t1091*t1114+t1061*t1090+t1095;
      coeffs_dyy->coeff_m2_m1 = t1065*t1115+t1029+t1108*t1059+t1107;
      coeffs_dyy->coeff_m1_m1 = t1020+t1038+t1087*x+RATIONAL(-35.0,24.0)*t1090+
(RATIONAL(-35.0,8.0)*y+RATIONAL(25.0,8.0))*t1091;
      coeffs_dyy->coeff_0_m1 = RATIONAL(5.0,2.0)+t1063*t1090+(RATIONAL(49.0,6.0
)*t1091+RATIONAL(-7.0,2.0))*y+t1106;
      coeffs_dyy->coeff_p1_m1 = RATIONAL(5.0,3.0)*x+t1037+RATIONAL(-10.0,3.0)*
t1090+y*t1050+t1099;
      coeffs_dyy->coeff_p2_m1 = RATIONAL(7.0,4.0)*t1115+t1077*x+RATIONAL(35.0,
24.0)*t1090+t1025+t1101;
      coeffs_dyy->coeff_p3_m1 = t1091*t1113+t1077*t1090+t1102;
      coeffs_dyy->coeff_m2_0 = t1036+t1081*t1115+t1108*t1069+t1097;
      coeffs_dyy->coeff_m1_0 = RATIONAL(10.0,1.0)*t1115+RATIONAL(28.0,9.0)*x+
t1024+RATIONAL(49.0,18.0)*t1090+t1106;
      coeffs_dyy->coeff_0_0 = t1064+RATIONAL(-56.0,9.0)*t1090+RATIONAL(98.0,9.0
)*t1091+t1023+(RATIONAL(-56.0,3.0)*t1091+RATIONAL(8.0,1.0))*y;
      coeffs_dyy->coeff_p1_0 = t1022+RATIONAL(-28.0,9.0)*x+RATIONAL(56.0,9.0)*
t1090+RATIONAL(40.0,3.0)*t1115+t1104;
      coeffs_dyy->coeff_p2_0 = t1035+RATIONAL(-49.0,18.0)*t1090+t1073*x+
RATIONAL(-4.0,1.0)*t1115+t1094;
      coeffs_dyy->coeff_p3_0 = t1073*t1090+t1091*t1112+t1093;
      coeffs_dyy->coeff_m2_p1 = t1083*t1115+t1035+t1108*t1058+t1096;
      coeffs_dyy->coeff_m1_p1 = RATIONAL(-35.0,18.0)*t1090+RATIONAL(-20.0,9.0)*
x+t1022+RATIONAL(-10.0,1.0)*t1115+t1099;
      coeffs_dyy->coeff_0_p1 = t1063+RATIONAL(40.0,9.0)*t1090+(RATIONAL(56.0,
3.0)*t1091+RATIONAL(-8.0,1.0))*y+t1104;
      coeffs_dyy->coeff_p1_p1 = t1023+t1024+RATIONAL(20.0,9.0)*x+RATIONAL(-40.0
,9.0)*t1090+(RATIONAL(50.0,9.0)+RATIONAL(-40.0,3.0)*y)*t1091;
      coeffs_dyy->coeff_p2_p1 = t1036+RATIONAL(35.0,18.0)*t1090+RATIONAL(4.0,
1.0)*t1115+t1067*x+t1100;
      coeffs_dyy->coeff_p3_p1 = t1091*t1111+t1067*t1090+t1103;
      coeffs_dyy->coeff_m2_p2 = t1076*t1115+t1025+t1108*t1066+t1098;
      coeffs_dyy->coeff_m1_p2 = RATIONAL(35.0,8.0)*t1115+t1037+t1080*x+t1065*
t1090+t1101;
      coeffs_dyy->coeff_0_p2 = t1081*t1090+RATIONAL(-1.0,1.0)+(RATIONAL(-49.0,
6.0)*t1091+RATIONAL(7.0,2.0))*y+t1094;
      coeffs_dyy->coeff_p1_p2 = RATIONAL(35.0,6.0)*t1115+t1083*t1090+t1088*x+
t1038+t1100;
      coeffs_dyy->coeff_p2_p2 = t1020+t1029+t1076*t1090+t1074*x+(RATIONAL(-7.0,
4.0)*y+t1079)*t1091;
      coeffs_dyy->coeff_p3_p2 = t1074*t1090+t1068*t1115+t1105;
      coeffs_dyy->coeff_m2_p3 = t1074*t1115+t1033+t1108*t1078+t1095;
      coeffs_dyy->coeff_m1_p3 = RATIONAL(-7.0,72.0)*t1090+t1041+RATIONAL(-1.0,
9.0)*x+RATIONAL(-5.0,8.0)*t1115+t1102;
      coeffs_dyy->coeff_0_p3 = t1084*t1090+t1089+(RATIONAL(7.0,6.0)*t1091+
RATIONAL(-1.0,2.0))*y+t1093;
      coeffs_dyy->coeff_p1_p3 = RATIONAL(1.0,9.0)*x+t1086*t1090+RATIONAL(-5.0,
6.0)*t1115+t1039+t1103;
      coeffs_dyy->coeff_p2_p3 = t1030+RATIONAL(1.0,4.0)*t1115+t1062*x+RATIONAL(
7.0,72.0)*t1090+t1105;
      coeffs_dyy->coeff_p3_p3 = t1062*t1090+t1031+(t1110+t1078)*t1091;
