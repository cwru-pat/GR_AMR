#include "ICs.h"
#include "../cosmo_includes.h"

using namespace SAMRAI;

namespace cosmo
{

ICsData cosmo_get_ICsData(
  std::shared_ptr<tbox::Database> cosmo_ICs_db, real_t domain_size)
{
  ICsData icd = {0};

  icd.rho_K_matter = 3.0/PI/8.0; // matter density term from FRW equation

  real_t rho_K_lambda_frac = cosmo_ICs_db->getDouble("rho_K_lambda_frac"); // DE density
  icd.rho_K_lambda = rho_K_lambda_frac*icd.rho_K_matter;

  // power spectrum amplitude as a fraction of the density
  real_t peak_amplitude_frac = cosmo_ICs_db->getDouble("peak_amplitude_frac"); // fluctuation amplitude
  real_t peak_amplitude = peak_amplitude_frac*(1.0e-15); // scaling in arb. units

  real_t ic_spec_cut = cosmo_ICs_db->getDouble("ic_spec_cut"); // power spectrum cutoff parameter

  /* (peak scale in hubble units) * (to pixel scale) */
  icd.peak_k = (1.0/0.07)*domain_size;
  icd.peak_amplitude = peak_amplitude;
  icd.ic_spec_cut = ic_spec_cut; // cut spectrum off around p ~ ic_spec_cut

  icd.viol_amp = cosmo_ICs_db->getDouble("IC_viol_amp");

  return icd;
}

/**
 * @brief A cosmologically-motivated power spectrum for the conformal factor
 * @details Return the power spectrum amplitude for a particular wavenumber.
 *  The power spectrum describing the density goes roughly as 
 *    P(k) ~ k^{ 1} for small k
 *    P(k) ~ k^{-3} for large k
 *  (see, eg, http://ned.ipac.caltech.edu/level5/Sept11/Norman/Norman2.html)
 *  An analytic function with this form is the function
 *    P(k) ~ k / ( 1 + k^4 )
 *  The conformal factor obeys k^2 \phi ~ \rho; so \phi ~ \rho / k^2.
 *  The power spectrum < \phi \phi > then goes as 1/k^4 < \rho \rho >.
 *  Thus, the power spectrum for the conformal factor scales as
 *    P(k) ~ 1 / k^3 / ( 1 + k^4 )
 *  Which after additional scalings, this function returns.
 * 
 * @param k wavenumber of interest
 * @param icd initial condition data structre
 * 
 * @return power spectrum amplitude
 */
real_t cosmo_power_spectrum(real_t k, ICsData *icd)
{
  real_t pre = icd->peak_amplitude;
  return pre/(1.0 + pow(fabs(k)/icd->peak_k, 4.0)/3.0)/pow(fabs(k)/icd->peak_k, 3.0);
}

// set a field to an arbitrary gaussian random field

void set_gaussian_random_field(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  idx_t ln, idx_t f_id, ICsData *icd)
{
  idx_t i, j, k;
  real_t px, py, pz, pmag;
  real_t scale;

  std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
    SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
      hierarchy->getGridGeometry()));
  TBOX_ASSERT(grid_geometry_);
  geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

  const double * domain_lower = &grid_geometry.getXLower()[0];
  const double * domain_upper = &grid_geometry.getXUpper()[0];

  real_t L[3];
  
  for(int i = 0 ; i < 3; i++)
    L[i] = domain_upper[i] - domain_lower[i];

  const double * dx = &grid_geometry.getDx()[0];

  idx_t NX = round(L[0] / dx[0]);
  idx_t NY = round(L[1] / dx[1]);
  idx_t NZ = round(L[2] / dx[2]);

  std::shared_ptr <hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
    
  
  // populate "field" with random values
  std::random_device rd;
  std::mt19937 gen(9.0 /*rd()*/);
  std::normal_distribution<real_t> gaussian_distribution;
  std::uniform_real_distribution<double> angular_distribution(0.0, 2.0*PI);
   // calling these here before looping suppresses a warning (bug)
  gaussian_distribution(gen);
  angular_distribution(gen);

  fftw_complex *f_field;

  
  double *r_field = new double[NX*NY*NZ];

  memset(r_field, 0, NX*NY*NZ*sizeof(double));

  
  f_field = (fftw_complex *) fftw_malloc(NX*NY*(NZ/2+1)
                                         *((long long) sizeof(fftw_complex)));
  // plans for taking FFTs
  fftw_plan p_c2r;

  p_c2r = fftw_plan_dft_c2r_3d(NX, NY, NZ,
                               f_field, r_field,
                               FFTW_MEASURE);
  // scale amplitudes in fourier space
  // don't expect to run at >512^3 anytime soon; loop over all momenta out to that resolution.
  // this won't work for a larger grid.
  idx_t NMAX = 512;
  real_t max_fft_index = NX*NY*(NZ/2+1) - 1;
  for(i=0; i<NMAX; i++)
  {
    px = (real_t) (i<=NMAX/2 ? i : i-NMAX);
    for(j=0; j<NMAX; j++)
    {
      py = (real_t) (j<=NMAX/2 ? j : j-NMAX);
      for(k=0; k<NMAX/2+1; k++)
      {
        pz = (real_t) k;

        // generate the same random modes for all resolutions (up to NMAX)
        real_t rand_mag = gaussian_distribution(gen);
        real_t rand_phase = angular_distribution(gen);

        // only store momentum values for relevant bins
        if( fabs(px) < (real_t) NX/2+1 + 0.01
            && fabs(py) < (real_t) NY/2+1 + 0.01
            && fabs(pz) < (real_t) NZ/2+1 + 0.01)
        {
          idx_t fft_index = FFT_NP_INDEX(
            px > -0.5 ? ROUND_2_IDXT(px) : (NX + ROUND_2_IDXT(px)),
            py > -0.5 ? ROUND_2_IDXT(py) : (NY + ROUND_2_IDXT(py)),
            pz > -0.5 ? ROUND_2_IDXT(pz) : (NZ + ROUND_2_IDXT(pz))
          );
 
          if(fft_index > max_fft_index)
          {
            tbox::pout<< "Warning: index " << fft_index << " is greater than max ("
                        << max_fft_index << ").\n";

            fft_index = max_fft_index;
          }

          pmag = sqrt(
            pw2(px * ( (real_t) NX / (real_t) NX ) )
             + pw2(py * ( (real_t) NX / (real_t) NY ))
             + pw2(pz * ( (real_t) NX / (real_t) NZ ))
            );

          // Scale by power spectrum
          // don't want much power on scales smaller than ~3 pixels
          // Or scales p > 1/(3*dx)
          real_t cutoff = 1.0 / (
              1.0 + exp(10.0*(pmag - icd->ic_spec_cut))
          );
          scale = cutoff*sqrt(cosmo_power_spectrum(pmag, icd));

          f_field[fft_index][0] = scale*rand_mag*cos(rand_phase);
          f_field[fft_index][1] = scale*rand_mag*sin(rand_phase);

        }
      }
    }
  }

  // zero-mode (mean density)... set this to something later
  f_field[FFT_NP_INDEX(0,0,0)][0] = 0;
  f_field[FFT_NP_INDEX(0,0,0)][1] = 0;

  // FFT back; 'field' array should now be populated with a gaussian random
  // field and power spectrum given by cosmo_power_spectrum.
  fftw_execute_dft_c2r(p_c2r, f_field, r_field);

  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {

    const std::shared_ptr<hier::Patch> & patch = *pit;



    const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch->getPatchGeometry()));


    std::shared_ptr<pdat::CellData<double> > f_pdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
        patch->getPatchData(f_id)));
    
    const hier::Box& box = f_pdata->getGhostBox();      

    MDA_Access<double, 3, MDA_OrderColMajor<3>> field =
      pdat::ArrayDataAccess::access<3,double>(
        f_pdata->getArrayData());


      
    const int * lower = &box.lower()[0];
    const int * upper = &box.upper()[0];
    
    for(k = lower[2]; k <= upper[2]; k++)
    {
      for(j = lower[1]; j <= upper[1]; j++)
      {
        for(i = lower[0]; i <= upper[0]; i++)
        {
          idx_t ii = (i + NX ) % NX, jj = (j + NY ) % NY, kk = (k + NZ ) % NZ; 
          idx_t fft_index = NP_INDEX(ii,jj,kk);
          field(i, j, k) = r_field[fft_index];
          if(tbox::MathUtilities< double >::isNaN(field(ii,jj,kk)))
          {
            TBOX_ERROR("Error: NaN field at "<<i<<" "<<j<<" "<<k<<"\n");
          }

        }
      }
    }
      
  }
  return;
}


} // namespace cosmo
