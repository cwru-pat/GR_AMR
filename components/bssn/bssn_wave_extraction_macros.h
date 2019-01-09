#ifndef BSSN_WAVE_MACROS
#define BSSN_WAVE_MACROS

#define BSSN_WAVE_CALCULATE_CHRISTOFFEL(I, J, K) GammaUDD[I-1][J-1][K-1] = bd.G##I##J##K \
    - (double)(I==J)*bd.d##K##chi/ bd.chi - (double)(I==K)*bd.d##J##chi/bd.chi \
    + bd.gamma##J##K *                                                  \
    (bd.gammai##I##1 * bd.d1chi + bd.gammai##I##2 * bd.d2chi + bd.gammai##I##3 * bd.d3chi) / bd.chi

#define BSSN_WAVE_CALCULATE_GAMMADD_dDD(I, J, M, N) gammaDD_dDD[M-1][N-1][I-1][J-1] =   ( \
    + 6.0 * bd.d##I##chi * bd.d##J##chi * bd.gamma##M##N / pw2(bd.chi)    \
    - 2.0 * bd.d##I##d##J##chi * bd.gamma##M##N / bd.chi                   \
    - 2.0 * bd.d##J##chi * bd.d##I##g##M##N / bd.chi                    \
    - 2.0 * bd.d##I##chi * bd.d##J##g##M##N / bd.chi                    \
    + bd.d##I##d##J##g##M##N                                               \
  )/pw2(bd.chi)


#define BSSN_WAVE_CALCULATE_KDD_dD(M, I, J) KDD_dD[I-1][J-1][M-1] = (         \
    -2.0 * bd.d##M##chi * K[I-1][J-1] * bd.chi                              \
    + derivative(bd.i, bd.j, bd.k, M, A##I##J##_a, dx)                  \
    + derivative(bd.i, bd.j, bd.k, M, DIFFK_a, dx) * bd.gamma##I##J /3.0 \
    + bd.DIFFK * bd.d##M##g##I##J /3.0                                  \
  ) / pw2(bd.chi)


#endif
