fp t1543;
fp t1525;
fp t1539;
fp t1578;
fp t1531;
fp t1519;
fp t1529;
fp t1577;
fp t1518;
fp t1538;
fp t1527;
fp t1576;
fp t1524;
fp t1528;
fp t1575;
fp t1574;
fp t1573;
fp t1521;
fp t1572;
fp t1571;
fp t1540;
fp t1520;
fp t1570;
fp t1569;
fp t1568;
fp t1567;
fp t1526;
fp t1566;
fp t1565;
fp t1564;
fp t1563;
fp t1562;
fp t1533;
fp t1561;
fp t1544;
fp t1560;
fp t1541;
fp t1559;
fp t1558;
fp t1557;
fp t1556;
fp t1555;
fp t1554;
fp t1553;
fp t1552;
fp t1551;
fp t1550;
fp t1535;
fp t1549;
fp t1548;
fp t1547;
fp t1546;
fp t1545;
fp t1542;
fp t1537;
fp t1534;
fp t1530;
fp t1523;
fp t1522;
      t1543 = RATIONAL(1.0,16.0);
      t1525 = t1543*x;
      t1539 = RATIONAL(1.0,80.0);
      t1578 = t1525+t1539;
      t1531 = RATIONAL(-3.0,80.0);
      t1519 = t1531*y;
      t1529 = t1531*z;
      t1577 = t1519+t1529;
      t1518 = t1539*y;
      t1538 = RATIONAL(-1.0,80.0);
      t1527 = t1538*z;
      t1576 = t1518+t1527;
      t1524 = t1538*y;
      t1528 = t1539*z;
      t1575 = t1524+t1528;
      t1574 = t1518+t1528;
      t1573 = t1518+t1529;
      t1521 = RATIONAL(3.0,16.0)*x;
      t1572 = t1518+t1521;
      t1571 = t1521+t1528;
      t1540 = RATIONAL(3.0,80.0);
      t1520 = t1540*y;
      t1570 = t1520+t1525;
      t1569 = t1524+t1529;
      t1568 = t1521+RATIONAL(-1.0,10.0);
      t1567 = t1519+t1527;
      t1526 = t1540*z;
      t1566 = t1525+t1526;
      t1565 = t1520+t1521;
      t1564 = t1524+t1527;
      t1563 = t1521+t1526;
      t1562 = t1519+t1525;
      t1533 = RATIONAL(3.0,40.0);
      t1561 = t1519+t1528+t1533;
      t1544 = RATIONAL(1.0,20.0);
      t1560 = t1544+t1574;
      t1541 = RATIONAL(1.0,40.0);
      t1559 = t1520+t1526+t1541;
      t1558 = t1520+t1527+t1544;
      t1557 = RATIONAL(1.0,10.0)+t1577;
      t1556 = t1524+t1526+t1544;
      t1555 = t1519+t1526+t1543;
      t1554 = t1533+t1573;
      t1553 = t1518+t1526+t1540;
      t1552 = t1520+t1528+t1540;
      t1551 = t1520+t1529+t1543;
      t1550 = t1543+t1575;
      t1535 = RATIONAL(7.0,80.0);
      t1549 = t1535+t1569;
      t1548 = t1533+t1564;
      t1547 = t1543+t1576;
      t1546 = t1535+t1567;
      t1545 = RATIONAL(-1.0,8.0);
      t1542 = RATIONAL(-1.0,40.0);
      t1537 = RATIONAL(-9.0,80.0);
      t1534 = RATIONAL(-3.0,20.0);
      t1530 = RATIONAL(-11.0,80.0);
      t1523 = RATIONAL(-3.0,16.0)*x;
      t1522 = RATIONAL(-1.0,16.0)*x;
      coeffs_dxx->coeff_m1_m1_m1 = t1522+t1557;
      coeffs_dxx->coeff_0_m1_m1 = RATIONAL(-13.0,80.0)+t1520+t1563;
      coeffs_dxx->coeff_p1_m1_m1 = t1523+t1559;
      coeffs_dxx->coeff_p2_m1_m1 = t1540+t1529+t1562;
      coeffs_dxx->coeff_m1_0_m1 = t1522+t1549;
      coeffs_dxx->coeff_0_0_m1 = t1534+t1518+t1563;
      coeffs_dxx->coeff_p1_0_m1 = t1523+t1553;
      coeffs_dxx->coeff_p2_0_m1 = t1525+t1541+t1569;
      coeffs_dxx->coeff_m1_p1_m1 = t1522+t1554;
      coeffs_dxx->coeff_0_p1_m1 = t1524+t1530+t1563;
      coeffs_dxx->coeff_p1_p1_m1 = t1523+t1556;
      coeffs_dxx->coeff_p2_p1_m1 = t1573+t1578;
      coeffs_dxx->coeff_m1_p2_m1 = t1522+t1551;
      coeffs_dxx->coeff_0_p2_m1 = t1545+t1519+t1563;
      coeffs_dxx->coeff_p1_p2_m1 = t1523+t1555;
      coeffs_dxx->coeff_p2_p2_m1 = t1529+t1570;
      coeffs_dxx->coeff_m1_m1_0 = t1522+t1546;
      coeffs_dxx->coeff_0_m1_0 = t1528+t1534+t1565;
      coeffs_dxx->coeff_p1_m1_0 = t1523+t1552;
      coeffs_dxx->coeff_p2_m1_0 = t1541+t1527+t1562;
      coeffs_dxx->coeff_m1_0_0 = t1522+t1548;
      coeffs_dxx->coeff_0_0_0 = t1530+t1518+t1571;
      coeffs_dxx->coeff_p1_0_0 = t1523+t1560;
      coeffs_dxx->coeff_p2_0_0 = t1564+t1578;
      coeffs_dxx->coeff_m1_p1_0 = t1522+t1547;
      coeffs_dxx->coeff_0_p1_0 = t1524+t1545+t1571;
      coeffs_dxx->coeff_p1_p1_0 = t1523+t1550;
      coeffs_dxx->coeff_p2_p1_0 = t1525+t1576;
      coeffs_dxx->coeff_m1_p2_0 = t1522+t1558;
      coeffs_dxx->coeff_0_p2_0 = t1537+t1519+t1571;
      coeffs_dxx->coeff_p1_p2_0 = t1523+t1561;
      coeffs_dxx->coeff_p2_p2_0 = t1527+t1538+t1570;
      coeffs_dxx->coeff_m1_m1_p1 = t1522+t1561;
      coeffs_dxx->coeff_0_m1_p1 = t1527+t1530+t1565;
      coeffs_dxx->coeff_p1_m1_p1 = t1523+t1558;
      coeffs_dxx->coeff_p2_m1_p1 = t1528+t1539+t1562;
      coeffs_dxx->coeff_m1_0_p1 = t1522+t1550;
      coeffs_dxx->coeff_0_0_p1 = t1545+t1527+t1572;
      coeffs_dxx->coeff_p1_0_p1 = t1523+t1547;
      coeffs_dxx->coeff_p2_0_p1 = t1525+t1575;
      coeffs_dxx->coeff_m1_p1_p1 = t1522+t1560;
      coeffs_dxx->coeff_0_p1_p1 = t1537+t1521+t1564;
      coeffs_dxx->coeff_p1_p1_p1 = t1523+t1548;
      coeffs_dxx->coeff_p2_p1_p1 = t1538+t1525+t1574;
      coeffs_dxx->coeff_m1_p2_p1 = t1522+t1552;
      coeffs_dxx->coeff_0_p2_p1 = t1567+t1568;
      coeffs_dxx->coeff_p1_p2_p1 = t1523+t1546;
      coeffs_dxx->coeff_p2_p2_p1 = t1528+t1542+t1570;
      coeffs_dxx->coeff_m1_m1_p2 = t1522+t1555;
      coeffs_dxx->coeff_0_m1_p2 = t1545+t1529+t1565;
      coeffs_dxx->coeff_p1_m1_p2 = t1523+t1551;
      coeffs_dxx->coeff_p2_m1_p2 = t1526+t1562;
      coeffs_dxx->coeff_m1_0_p2 = t1522+t1556;
      coeffs_dxx->coeff_0_0_p2 = t1529+t1537+t1572;
      coeffs_dxx->coeff_p1_0_p2 = t1523+t1554;
      coeffs_dxx->coeff_p2_0_p2 = t1538+t1524+t1566;
      coeffs_dxx->coeff_m1_p1_p2 = t1522+t1553;
      coeffs_dxx->coeff_0_p1_p2 = t1568+t1569;
      coeffs_dxx->coeff_p1_p1_p2 = t1523+t1549;
      coeffs_dxx->coeff_p2_p1_p2 = t1542+t1518+t1566;
      coeffs_dxx->coeff_m1_p2_p2 = t1522+t1559;
      coeffs_dxx->coeff_0_p2_p2 = RATIONAL(-7.0,80.0)+t1521+t1577;
      coeffs_dxx->coeff_p1_p2_p2 = t1523+t1557;
      coeffs_dxx->coeff_p2_p2_p2 = t1531+t1520+t1566;
