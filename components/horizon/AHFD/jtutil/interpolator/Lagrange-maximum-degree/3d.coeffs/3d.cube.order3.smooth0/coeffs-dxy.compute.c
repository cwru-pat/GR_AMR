fp t1622;
fp t1588;
fp t1590;
fp t1660;
fp t1611;
fp t1584;
fp t1659;
fp t1621;
fp t1587;
fp t1658;
fp t1583;
fp t1657;
fp t1613;
fp t1585;
fp t1656;
fp t1586;
fp t1599;
fp t1655;
fp t1581;
fp t1593;
fp t1654;
fp t1653;
fp t1623;
fp t1594;
fp t1652;
fp t1651;
fp t1650;
fp t1649;
fp t1648;
fp t1647;
fp t1603;
fp t1580;
fp t1646;
fp t1645;
fp t1589;
fp t1644;
fp t1643;
fp t1642;
fp t1620;
fp t1591;
fp t1592;
fp t1641;
fp t1640;
fp t1639;
fp t1638;
fp t1637;
fp t1636;
fp t1635;
fp t1634;
fp t1633;
fp t1632;
fp t1579;
fp t1631;
fp t1630;
fp t1629;
fp t1628;
fp t1627;
fp t1607;
fp t1582;
fp t1626;
fp t1625;
fp t1624;
fp t1618;
fp t1616;
fp t1612;
fp t1609;
fp t1604;
fp t1602;
fp t1601;
fp t1598;
fp t1597;
fp t1595;
      t1622 = RATIONAL(1.0,80.0);
      t1588 = t1622*y;
      t1590 = t1622*x;
      t1660 = t1588+t1590;
      t1611 = RATIONAL(-3.0,1000.0);
      t1584 = t1611*z;
      t1659 = t1584+RATIONAL(-31.0,1000.0);
      t1621 = RATIONAL(-1.0,80.0);
      t1587 = t1621*x;
      t1658 = t1587+t1588;
      t1583 = RATIONAL(3.0,1000.0)*z;
      t1657 = t1583+RATIONAL(37.0,2000.0);
      t1613 = RATIONAL(-9.0,1000.0);
      t1585 = t1613*z;
      t1656 = t1585+RATIONAL(49.0,2000.0);
      t1586 = RATIONAL(9.0,1000.0)*z;
      t1599 = RATIONAL(-19.0,2000.0);
      t1655 = t1586+t1599;
      t1581 = RATIONAL(1.0,1000.0)*z;
      t1593 = t1621*y;
      t1654 = t1581+t1593;
      t1653 = t1583+RATIONAL(-17.0,500.0);
      t1623 = RATIONAL(3.0,80.0);
      t1594 = t1623*x;
      t1652 = t1594+t1593;
      t1651 = t1583+RATIONAL(-1.0,250.0);
      t1650 = t1584+t1593;
      t1649 = t1584+t1587;
      t1648 = t1583+RATIONAL(2.0,125.0);
      t1647 = t1581+t1588;
      t1603 = RATIONAL(-27.0,1000.0);
      t1580 = t1603*z;
      t1646 = t1580+t1613;
      t1645 = t1585+RATIONAL(-9.0,500.0);
      t1589 = t1623*y;
      t1644 = t1589+t1594;
      t1643 = t1585+RATIONAL(-7.0,250.0);
      t1642 = t1585+t1594;
      t1620 = RATIONAL(-3.0,80.0);
      t1591 = t1620*x;
      t1592 = t1620*y;
      t1641 = t1591+t1592;
      t1640 = t1589+t1590;
      t1639 = t1586+t1589;
      t1638 = t1588+t1594;
      t1637 = t1592+t1594;
      t1636 = t1590+t1592;
      t1635 = t1588+t1591;
      t1634 = t1591+t1593;
      t1633 = t1584+t1591;
      t1632 = t1589+t1591;
      t1579 = RATIONAL(27.0,1000.0)*z;
      t1631 = t1592+t1579;
      t1630 = t1587+t1592;
      t1629 = t1586+t1592;
      t1628 = t1585+RATIONAL(11.0,500.0);
      t1627 = t1586+t1591;
      t1607 = RATIONAL(-1.0,1000.0);
      t1582 = t1607*z;
      t1626 = t1582+t1593;
      t1625 = t1583+RATIONAL(-13.0,2000.0);
      t1624 = t1587+t1589;
      t1618 = RATIONAL(-1.0,500.0);
      t1616 = RATIONAL(-9.0,250.0);
      t1612 = RATIONAL(-1.0,2000.0);
      t1609 = RATIONAL(-7.0,2000.0);
      t1604 = RATIONAL(19.0,1000.0);
      t1602 = RATIONAL(-37.0,1000.0);
      t1601 = RATIONAL(13.0,1000.0);
      t1598 = RATIONAL(43.0,2000.0);
      t1597 = RATIONAL(31.0,2000.0);
      t1595 = RATIONAL(-21.0,2000.0);
      coeffs_dxy->coeff_m1_m1_m1 = t1580+RATIONAL(147.0,2000.0)+t1641;
      coeffs_dxy->coeff_0_m1_m1 = t1593+t1612+t1642;
      coeffs_dxy->coeff_p1_m1_m1 = t1586+t1602+t1638;
      coeffs_dxy->coeff_p2_m1_m1 = t1579+t1616+t1632;
      coeffs_dxy->coeff_m1_0_m1 = t1585+t1612+t1624;
      coeffs_dxy->coeff_0_0_m1 = t1584+RATIONAL(-17.0,2000.0)+t1660;
      coeffs_dxy->coeff_p1_0_m1 = t1593+t1590+t1651;
      coeffs_dxy->coeff_p2_0_m1 = t1587+t1601+t1629;
      coeffs_dxy->coeff_m1_p1_m1 = t1602+t1590+t1639;
      coeffs_dxy->coeff_0_p1_m1 = t1651+t1658;
      coeffs_dxy->coeff_p1_p1_m1 = RATIONAL(33.0,2000.0)+t1593+t1649;
      coeffs_dxy->coeff_p2_p1_m1 = t1636+t1656;
      coeffs_dxy->coeff_m1_p2_m1 = t1594+t1616+t1631;
      coeffs_dxy->coeff_0_p2_m1 = t1601+t1593+t1627;
      coeffs_dxy->coeff_p1_p2_m1 = t1635+t1656;
      coeffs_dxy->coeff_p2_p2_m1 = RATIONAL(-3.0,2000.0)+t1580+t1644;
      coeffs_dxy->coeff_m1_m1_0 = RATIONAL(129.0,2000.0)+t1585+t1641;
      coeffs_dxy->coeff_0_m1_0 = t1609+t1594+t1650;
      coeffs_dxy->coeff_p1_m1_0 = t1638+t1653;
      coeffs_dxy->coeff_p2_m1_0 = t1589+t1603+t1627;
      coeffs_dxy->coeff_m1_0_0 = t1609+t1584+t1624;
      coeffs_dxy->coeff_0_0_0 = t1599+t1582+t1660;
      coeffs_dxy->coeff_p1_0_0 = t1611+t1590+t1654;
      coeffs_dxy->coeff_p2_0_0 = t1630+t1648;
      coeffs_dxy->coeff_m1_p1_0 = t1640+t1653;
      coeffs_dxy->coeff_0_p1_0 = t1587+t1611+t1647;
      coeffs_dxy->coeff_p1_p1_0 = t1587+t1597+t1626;
      coeffs_dxy->coeff_p2_p1_0 = t1598+t1584+t1636;
      coeffs_dxy->coeff_m1_p2_0 = t1603+t1594+t1629;
      coeffs_dxy->coeff_0_p2_0 = t1634+t1648;
      coeffs_dxy->coeff_p1_p2_0 = t1588+t1598+t1633;
      coeffs_dxy->coeff_p2_p2_0 = t1595+t1589+t1642;
      coeffs_dxy->coeff_m1_m1_p1 = t1592+RATIONAL(111.0,2000.0)+t1627;
      coeffs_dxy->coeff_0_m1_p1 = t1625+t1652;
      coeffs_dxy->coeff_p1_m1_p1 = t1638+t1659;
      coeffs_dxy->coeff_p2_m1_p1 = t1632+t1645;
      coeffs_dxy->coeff_m1_0_p1 = t1624+t1625;
      coeffs_dxy->coeff_0_0_p1 = t1590+t1595+t1647;
      coeffs_dxy->coeff_p1_0_p1 = t1618+t1590+t1626;
      coeffs_dxy->coeff_p2_0_p1 = t1584+t1604+t1630;
      coeffs_dxy->coeff_m1_p1_p1 = t1640+t1659;
      coeffs_dxy->coeff_0_p1_p1 = t1582+t1618+t1658;
      coeffs_dxy->coeff_p1_p1_p1 = t1587+RATIONAL(29.0,2000.0)+t1654;
      coeffs_dxy->coeff_p2_p1_p1 = t1636+t1657;
      coeffs_dxy->coeff_m1_p2_p1 = t1637+t1645;
      coeffs_dxy->coeff_0_p2_p1 = t1604+t1593+t1633;
      coeffs_dxy->coeff_p1_p2_p1 = t1635+t1657;
      coeffs_dxy->coeff_p2_p2_p1 = t1594+RATIONAL(-39.0,2000.0)+t1639;
      coeffs_dxy->coeff_m1_m1_p2 = t1591+RATIONAL(93.0,2000.0)+t1631;
      coeffs_dxy->coeff_0_m1_p2 = t1652+t1655;
      coeffs_dxy->coeff_p1_m1_p2 = t1638+t1643;
      coeffs_dxy->coeff_p2_m1_p2 = t1632+t1646;
      coeffs_dxy->coeff_m1_0_p2 = t1624+t1655;
      coeffs_dxy->coeff_0_0_p2 = RATIONAL(-23.0,2000.0)+t1583+t1660;
      coeffs_dxy->coeff_p1_0_p2 = t1590+t1607+t1650;
      coeffs_dxy->coeff_p2_0_p2 = t1628+t1630;
      coeffs_dxy->coeff_m1_p1_p2 = t1640+t1643;
      coeffs_dxy->coeff_0_p1_p2 = t1607+t1588+t1649;
      coeffs_dxy->coeff_p1_p1_p2 = t1593+RATIONAL(27.0,2000.0)+t1583+t1587;
      coeffs_dxy->coeff_p2_p1_p2 = t1590+t1597+t1629;
      coeffs_dxy->coeff_m1_p2_p2 = t1637+t1646;
      coeffs_dxy->coeff_0_p2_p2 = t1628+t1634;
      coeffs_dxy->coeff_p1_p2_p2 = t1597+t1588+t1627;
      coeffs_dxy->coeff_p2_p2_p2 = t1579+RATIONAL(-57.0,2000.0)+t1644;
