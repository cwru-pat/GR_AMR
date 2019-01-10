#include "bssn.h"
#include "../../cosmo_includes.h"
#include "../../utils/math.h"
#include "bssn_wave_extraction_macros.h"

using namespace SAMRAI;

namespace cosmo
{
// calculating Weyl scalars \Psi0 ~ \Psi4
// following the way in
// https://github.com/zachetienne/nrpytutorial/blob/a7f1f5778be2a228fcf5bd9c29d602e6f384dd93/Tutorial-WeylScalarsInvariants-Cartesian.ipynb
void BSSN::cal_Weyl_scalars(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  int weight_idx)
{
  #if CAL_WEYL_SCALS

  // defining constant arrays
  double leviCvt[DIM][DIM][DIM];
  double tot_energy = 0;
  for(int i = 0; i < DIM; i ++)
    for(int j = 0; j < DIM; j ++)
      for(int k = 0; k < DIM; k++)
        leviCvt[i][j][k] = (double)(i-j)*(double)(j-k)*(double)(k-i)/2.0;
  
  for(int ln = 0; ln < hierarchy->getNumberOfLevels(); ln ++)
  {
    std::shared_ptr <hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    
    for( hier::PatchLevel::iterator pit(level->begin());
         pit != level->end(); ++pit)
    {
      const std::shared_ptr<hier::Patch> & patch = *pit;

      const hier::Box& box = patch->getBox();

      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));
      
      std::shared_ptr<pdat::CellData<double> > weight(
        SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(weight_idx)));
      

      arr_t weight_array =
        pdat::ArrayDataAccess::access<DIM, double>(
          weight->getArrayData());


      initPData(patch);

      initMDA(patch);
      
      const int * lower = &box.lower()[0];
      const int * upper = &box.upper()[0];
      
      const double *dx = &patch_geom->getDx()[0];

      idx_t NX = round(L[0] / dx[0]);
      idx_t NY = round(L[1] / dx[1]); 
      idx_t NZ = round(L[2] / dx[2]);

      for(int k = lower[2]; k <= upper[2]; k++)
      {
        for(int j = lower[1]; j <= upper[1]; j++)
        {
          for(int i = lower[0]; i <= upper[0]; i++)
          {
            double x = (dx[0] * ((real_t)i + 0.5)) - L[0] / 2.0 ;
            double y = (dx[1] * ((real_t)j + 0.5)) - L[1] / 2.0 ;
            double z = (dx[2] * ((real_t)k + 0.5)) - L[2] / 2.0 ;
            double r = sqrt(pw2(x) + pw2(y) + pw2(z));

            BSSNData bd = {0};

            // 3 orthogonal vectors
            double v1[3] = {0}, v2[3] = {0}, v3[3] = {0};
            double w1[3] = {0}, w2[3] = {0}, w3[3] = {0};
            double e1[3] = {0}, e2[3] = {0}, e3[3] = {0};

            double gamma[3][3] = {0}, gammai[3][3] = {0};
            double K[3][3] = {0}, KDD_dD[DIM][DIM][DIM];
            
            set_bd_values(i,j,k,&bd,dx);
            
            v1[0] = -y, v1[1] = x, v1[2] = 0;
            v2[0] = x, v2[1] = y, v2[2] = z;

            // initializing gamma and gammai
            // putting them in to arrays makes it easier
            // to sum over
            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
              {
                gamma[0][0] = bd.gamma11 / pw2(bd.chi);
                gamma[0][1] = gamma[1][0] = bd.gamma12 / pw2(bd.chi);
                gamma[0][2] = gamma[2][0] = bd.gamma13 / pw2(bd.chi);
                gamma[1][1] = bd.gamma22 / pw2(bd.chi);
                gamma[1][2] = gamma[2][1] = bd.gamma23 / pw2(bd.chi);
                gamma[2][2] = bd.gamma33 / pw2(bd.chi);

                gammai[0][0] = bd.gammai11 * pw2(bd.chi);
                gammai[0][1] = gammai[1][0] = bd.gammai12 * pw2(bd.chi);
                gammai[0][2] = gammai[2][0] = bd.gammai13 * pw2(bd.chi);
                gammai[1][1] = bd.gammai22 * pw2(bd.chi);
                gammai[1][2] = gammai[2][1] = bd.gammai23 * pw2(bd.chi);
                gammai[2][2] = bd.gammai33 * pw2(bd.chi);

                K[0][0] = (bd.A11 + bd.gamma11 * bd.K / 3.0) / pw2(bd.chi);
                K[1][1] = (bd.A22 + bd.gamma22 * bd.K / 3.0) / pw2(bd.chi);
                K[2][2] = (bd.A33 + bd.gamma33 * bd.K / 3.0) / pw2(bd.chi);
                K[0][1] = K[1][0] = (bd.A12 + bd.gamma12 * bd.K / 3.0) / pw2(bd.chi);
                K[0][2] = K[2][0] = (bd.A13 + bd.gamma13 * bd.K / 3.0) / pw2(bd.chi);
                K[1][2] = K[2][1] = (bd.A23 + bd.gamma23 * bd.K / 3.0) / pw2(bd.chi);
              }
            
            COSMO_APPLY_TO_IJK_PERMS(BSSN_WAVE_CALCULATE_KDD_dD);

            for(int a = 0; a < DIM; a++)
              for(int b = a+1; b < DIM; b++)
                for(int c = 0; c < DIM; c++)
                  KDD_dD[b][a][c] = KDD_dD[a][b][c];
            
            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
                for(int c = 0; c < DIM; c++)
                  for(int d = 0; d < DIM; d++)
                  {
                    v3[a] += 1.0 / pw3(bd.chi) * gammai[a][d] * leviCvt[d][b][c]
                      * v1[b] * v2[c];
                  }
            
            // Gram-Schmidt orthonormalization of vectors
            double omega11=0, omega12=0, omega13=0, omega22=0,omega23=0,omega33=0;
            
            for(int a = 0; a < DIM; a++)
              w1[a] = v1[a];
            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
                omega11 += w1[a] * w1[b] * gamma[a][b];
            for(int a = 0; a < DIM; a++)
              e1[a] = w1[a] / sqrt(omega11);

            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
                omega12 += e1[a] * v2[b] * gamma[a][b];
            for(int a = 0; a < DIM; a++)
              w2[a] = v2[a] - omega12 * e1[a];
            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
                omega22 += w2[a] * w2[b] * gamma[a][b];
            for(int a = 0; a < DIM; a++)
              e2[a] = w2[a] / sqrt(omega22);

            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
                omega13 += e1[a] * v3[b] * gamma[a][b];
            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
                omega23 += e2[a] * v3[b] * gamma[a][b];
            for(int a = 0; a < DIM; a++)
              w3[a] = v3[a] - omega13 * e1[a] - omega23 * e2[a];
            
            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
                omega33 += w3[a] * w3[b] * gamma[a][b];
            for(int a = 0; a < DIM; a++)
              e3[a] = w3[a] / sqrt(omega33);

            double isqrt2 = 1.0 / sqrt(2);
            double ltet[DIM] = {0}, ntet[DIM] = {0}, remtet[DIM] = {0}, immtet[DIM] = {0};

            for(int a = 0; a < DIM; a ++)
            {
              ltet[a] = isqrt2 * e2[a];
              ntet[a] = -isqrt2 * e2[a];
              remtet[a] = isqrt2 * e3[a];
              immtet[a] = isqrt2 * e1[a];
            }

            
            // calculating Christoffle from the conformal one!!!!
            double GammaUDD[DIM][DIM][DIM] = {0};

            COSMO_APPLY_TO_IJK_PERMS(BSSN_WAVE_CALCULATE_CHRISTOFFEL);
            
            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
                for(int c = b+1; c < DIM; c++)
                {
                  GammaUDD[a][c][b] = GammaUDD[a][b][c]; 
                }

            // calculating R_{ijkl}
            double gammaDD_dDD[DIM][DIM][DIM][DIM] = {0};
            double Riemann[DIM][DIM][DIM][DIM] = {0};

            COSMO_APPLY_TO_IJMN_PERMS(BSSN_WAVE_CALCULATE_GAMMADD_dDD);

            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b ++)
              {
                for(int c = 0; c < DIM; c++)
                  for(int d = 0; d < DIM; d++)
                  {
                    int aa = a, bb = b, cc = c, dd = d;                   
                    if(b < a)
                      std::swap(aa, bb);
                    if(d < c)
                      std::swap(dd, cc);
                    gammaDD_dDD[a][b][c][d] = gammaDD_dDD[aa][bb][cc][dd];
                  }
                   
              }

            
            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
                for(int c = 0; c < DIM; c++)
                  for(int d = 0; d < DIM; d++)
                  {
                    Riemann[a][b][c][d] = (  gammaDD_dDD[a][d][c][b]
                                           + gammaDD_dDD[b][c][d][a]
                                           - gammaDD_dDD[a][c][b][d]
                                           - gammaDD_dDD[b][d][a][c]
                    ) / 2.0;
                    for(int m = 0; m < DIM; m++)
                      for(int n = 0; n < DIM; n++)
                        Riemann[a][b][c][d] +=
                          gamma[n][m] * GammaUDD[n][b][c] * GammaUDD[m][a][d]
                          - gamma[n][m] * GammaUDD[n][b][d] * GammaUDD[m][a][c];
                  }
            double CodazziDDD[DIM][DIM][DIM] = {0};
            double GaussDDDD[DIM][DIM][DIM][DIM] = {0};
            double RojoDD[DIM][DIM] = {0};

            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
                for(int c = 0; c < DIM; c++)
                  for(int d = 0; d < DIM; d++)
                    GaussDDDD[a][b][c][d] =
                      Riemann[a][b][c][d]
                      + K[a][c] * K[d][b] - K[a][d] * K [c][b];

            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
                for(int c = 0; c < DIM; c++)
                {
                  CodazziDDD[a][b][c] = KDD_dD[a][c][b] - KDD_dD[a][b][c];
                  for(int d = 0; d < DIM; d++)
                    CodazziDDD[a][b][c] +=
                      GammaUDD[d][a][c] * K[b][d] - GammaUDD[d][a][b] * K[c][d];
                }
            
            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
              {
                RojoDD[a][b] = bd.K * K[a][b];
                for(int c = 0; c < DIM; c++)
                  for(int d = 0; d < DIM; d++)
                    RojoDD[a][b] += gammai[c][d] * Riemann[a][c][b][d]
                      - K[a][c] * gammai[c][d] * K[d][b];
              }

            // Psi0r_a(i, j, k) = 0, Psi0i_a(i, j, k) = 0;
            // Psi1r_a(i, j, k) = 0, Psi1i_a(i, j, k) = 0;
            // Psi2r_a(i, j, k) = 0, Psi2i_a(i, j, k) = 0;
            // Psi3r_a(i, j, k) = 0, Psi3i_a(i, j, k) = 0;
            // Psi4r_a(i, j, k) = 0, Psi4i_a(i, j, k) = 0;
            
            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
              {
                Psi4r_a(i, j, k) =
                  0.5 * RojoDD[b][a] * (remtet[b] * remtet[a] - immtet[b] * immtet[a]);
                Psi4i_a(i, j, k) =
                  0.5 * RojoDD[b][a] * (-remtet[b] * immtet[a] - immtet[b] * remtet[a]);

                Psi3r_a(i, j, k) =
                  - 0.5 * RojoDD[b][a] * (ntet[b] - ltet[b]) * remtet[a];
                Psi3i_a(i, j, k) =
                  0.5 * RojoDD[b][a] * (ntet[b] - ltet[b]) * immtet[a];

                Psi2r_a(i, j, k) =
                  - 0.5 * RojoDD[b][a] * (remtet[a] * remtet[b] + immtet[b] * immtet[a]);
                Psi2i_a(i, j, k) =
                  - 0.5 * RojoDD[b][a] * (immtet[a] * remtet[b] - remtet[b] * immtet[a]);

                Psi1r_a(i, j, k) =
                  0.5 * RojoDD[b][a] * (ntet[b] * remtet[a] - ltet[b] * remtet[a]);
                Psi1i_a(i, j, k) =
                  0.5 * RojoDD[b][a] * (ntet[b] * immtet[a] - ltet[b] * immtet[a]);

                Psi0r_a(i, j, k) =
                  0.5 * RojoDD[b][a] * (remtet[b] * remtet[a] - immtet[b] * immtet[a]);
                Psi0i_a(i, j, k) =
                  0.5 * RojoDD[b][a] * (remtet[b] * immtet[a] + immtet[b] * remtet[a]);
                
              }


            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
                for(int c = 0; c < DIM; c++)
                {
                  Psi4r_a(i, j, k) += sqrt(2) * CodazziDDD[b][c][a] * ntet[c]
                    * (remtet[b] * remtet[a] - immtet[b] * immtet[a]);
                  Psi4i_a(i, j, k) += sqrt(2) * CodazziDDD[b][c][a] * ntet[c]
                    * (-remtet[b] * immtet[a] - immtet[b] * remtet[a]);

                  Psi3r_a(i, j, k) += 1.0 * CodazziDDD[b][c][a]
                    * ((ntet[b] - ltet[b])
                       * remtet[c] * ntet[a] - remtet[b] * ltet[c] * ntet[a]);
                  Psi3i_a(i, j, k) += -1.0 * CodazziDDD[b][c][a] 
                    * ((ntet[b] - ltet[b])
                       * immtet[c] * ntet[a] - immtet[b] * ltet[c] * ntet[a]);

                  Psi2r_a(i, j, k) +=  CodazziDDD[b][c][a] * isqrt2
                    * (ntet[a] * (remtet[b] * remtet[c] + immtet[b] * immtet[c])
                       - ltet[c] * (remtet[b] * remtet[a] + immtet[b] * immtet[a]));
                  Psi2i_a(i, j, k) += isqrt2 * CodazziDDD[b][c][a]
                    *(ntet[a] * (immtet[b] * remtet[c] - remtet[b] * immtet[c])
                      - ltet[c] * (remtet[b] * immtet[a] - immtet[b] * remtet[a]));

                  Psi1r_a(i, j, k) += isqrt2 * CodazziDDD[b][c][a]
                    *(ltet[b] * remtet[c] * ltet[a] - remtet[b] * ntet[c] * ltet[a]
                      - ntet[b] * remtet[c] * ltet[a]);
                  Psi1i_a(i, j, k) += isqrt2 * CodazziDDD[b][c][a]
                    * (ltet[b] * immtet[c] * ltet[a] - immtet[b] * ntet[c] * ltet[a]
                       - ntet[b] * immtet[c] * ltet[a]);

                  Psi0r_a(i, j, k) += sqrt(2) * CodazziDDD[b][c][a] * ltet[c]
                    * (remtet[b] * remtet[a] - immtet[b] * immtet[a]);
                  Psi0i_a(i, j, k) += sqrt(2) * CodazziDDD[b][c][a] * ltet[c]
                    * (remtet[b] * immtet[a] + immtet[b] * remtet[a]);
                }

            for(int a = 0; a < DIM; a++)
              for(int b = 0; b < DIM; b++)
                for(int c = 0; c < DIM; c++)
                  for(int d = 0; d < DIM; d++)
                  {
                    Psi4r_a(i, j, k) += GaussDDDD[d][b][c][a] * ntet[d] * ntet[c]
                      * (remtet[b] * remtet[a] - immtet[b] * immtet[a]);
                    Psi4i_a(i, j, k) += GaussDDDD[d][b][c][a] * ntet[d] * ntet[c]
                      * (-remtet[b] * immtet[a] - immtet[b] * remtet[a]);

                    Psi3r_a(i, j, k) += GaussDDDD[d][b][c][a] * ltet[d] * ntet[b]
                      * remtet[c] * ntet[a];
                    Psi3i_a(i, j, k) += -GaussDDDD[d][b][c][a] * ltet[d] * ntet[b]
                      * immtet[c] * ntet[a];

                    Psi2r_a(i, j, k) += GaussDDDD[d][b][c][a] * ltet[d] * ntet[a]
                      * (remtet[b] * remtet[c] + immtet[b] * immtet[c]);
                    Psi2i_a(i, j, k) += GaussDDDD[d][b][c][a] * ltet[d] * ntet[a]
                      * (immtet[b] * remtet[c] - remtet[b] * immtet[c]);

                    Psi1r_a(i, j, k) += GaussDDDD[d][b][c][a] * ntet[d] * ltet[b]
                      * remtet[c] * ltet[a];
                    Psi1i_a(i, j, k) += GaussDDDD[d][b][c][a] * ntet[d] * ltet[b]
                      * immtet[c] * ltet[a];

                    Psi0r_a(i, j, k) += GaussDDDD[d][b][c][a] * ltet[d] * ltet[c]
                      * (remtet[b] * remtet[a] - immtet[b] * immtet[a]);
                    Psi0i_a(i, j, k) += GaussDDDD[d][b][c][a] * ltet[d] * ltet[c]
                      * (remtet[b] * immtet[a] + immtet[b] * remtet[a]);
                    
                  }

            tot_energy += 1.0 / pw3(bd.chi) * sqrt(pw2(Psi4r_a(i, j, k)) + pw2(Psi4i_a(i, j, k))) * weight_array(i, j, k);
          }
        }
      }

      
    }
  }

  const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
  mpi.Barrier();
  if (mpi.getSize() > 1) {
    mpi.AllReduce(&tot_energy, 1, MPI_SUM);
  }

  tbox::pout<<"Total energy of GWs is "<<tot_energy<<"\n";
  #endif
}




}
