#include "../../cosmo_includes.h"
#include "geodesic_ic.h"
#include "SAMRAI/tbox/HDFDatabase.h"
#include "../elliptic_solver/full_multigrid.h"
#include "../elliptic_solver/multigrid_bd_handler.h"
#include "../../utils/Array.h"
#include <iostream>
#include <fstream>

#include "../horizon/AHFD/jtutil/util.hh"
#include "../horizon/AHFD/jtutil/array.hh"
#include "../horizon/AHFD/jtutil/util_String.h"
#include "../horizon/AHFD/jtutil/util_Table.h"
#include "../horizon/AHFD/jtutil/interpolator/InterpLocalUniform.h"
//#include "../horizon/AHFD/AHFD_types.h"
#include "../horizon/AHFD/AHFD_macros.h"

#include "geodesic.h"

using namespace SAMRAI;

namespace cosmo
{

  void Geodesic::geodesic_ic_transform_inertial_vectors(
    const std::shared_ptr<hier::Patch> & patch,
    BSSN *bssn,
    real_t ox, real_t oy, real_t oz,
    std::vector<real_t> &dirx, std::vector<real_t> &diry, std::vector<real_t> &dirz,
    const real_t dx[], double epsilon, int &p_id)
  {
    initPData(patch);
    bssn->initPData(patch);
    bssn->initMDA(patch);
    
    int i0 = floor((ox - domain_lower[0] ) / dx[0] );
    int j0 = floor((oy - domain_lower[1] ) / dx[1] );
    int k0 = floor((oz - domain_lower[2] ) / dx[2] );

    hier::Index *ic = new hier::Index(i0, j0, k0);
    
    real_t x0 = domain_lower[0] + (double)i0 * dx[0] + dx[0]/2.0;
    real_t y0 = domain_lower[1] + (double)j0 * dx[1] + dx[1]/2.0;
    real_t z0 = domain_lower[2] + (double)k0 * dx[2] + dx[2]/2.0;

    double xd = (ox - x0) / dx[0];
    double yd = (oy - y0) / dx[1];
    double zd = (oz - z0) / dx[2];

    GeodesicData gd_d = {0};
    GeodesicData *gd = &gd_d;

    GEODESIC_DEFINE_CRSPLINES_CHI;
    COSMO_APPLY_TO_IJ_PERMS(GEODESIC_DEFINE_CRSPLINES_M);
    BSSNData bd = {0};

    for(int i = 0; i < 4; i++)
      for(int j = 0; j < 4; j++)
        for(int k = 0; k < 4; k++)
        {
#if USE_COSMOTRACE
          bssn->set_bd_values_for_ray_tracing(i0-1+i, j0-1+j, k0-1+k, &bd, dx);
#endif
          GEODESIC_CRSPLINES_SET_F_CHI;
          COSMO_APPLY_TO_IJ_PERMS(GEODESIC_CRSPLINES_SET_F_M);
        }
    GEODESIC_CRSPLINES_CAL_COEF_CHI;
    COSMO_APPLY_TO_IJ_PERMS(GEODESIC_CRSPLINES_CAL_COEF_M);

    GEODESIC_CRSPLINES_EVAL_CHI;
    COSMO_APPLY_TO_IJ_PERMS(GEODESIC_CRSPLINES_EVAL_M);


    gd->m11 = gd->m11 / pw2(gd->chi);
    gd->m12 = gd->m12 / pw2(gd->chi);
    gd->m13 = gd->m13 / pw2(gd->chi);
    gd->m22 = gd->m22 / pw2(gd->chi);
    gd->m23 = gd->m23 / pw2(gd->chi);
    gd->m33 = gd->m33 / pw2(gd->chi);

    real_t det = gd->m11 * gd->m22 * gd->m33 + gd->m12 * gd->m23 * gd->m13
      + gd->m12 * gd->m23 * gd->m13 - gd->m13 * gd->m22 * gd->m13
      - gd->m12 * gd->m12 * gd->m33 - gd->m23 * gd->m23 * gd->m11;

    
    gd->mi11 = (gd->m22 * gd->m33 - pw2(gd->m23)) / det;
    gd->mi22 = (gd->m11 * gd->m33 - pw2(gd->m13)) / det;
    gd->mi33 = (gd->m11 * gd->m22 - pw2(gd->m12)) / det;
    gd->mi12 = (gd->m13*gd->m23 - gd->m12*(gd->m33)) / det;
    gd->mi13 = (gd->m12*gd->m23 - gd->m13*(gd->m22)) / det;
    gd->mi23 = (gd->m12*gd->m13 - gd->m23*(gd->m11)) / det;

    // real_t B = gd->mi11 * (gd->mi11 * gd->mi22 - pw2(gd->mi12));
    // real_t C = gd->mi11 * (gd->mi22 * gd->mi33 - pw2(gd->mi23))
    //   - pw2(gd->mi12) * gd->mi33 + 2 * gd->mi12 * gd->mi13 * gd->mi23 - pw2(gd->mi13) * gd->mi22;
    // real_t D = gd->mi11 * gd->mi22 - pw2(gd->mi23);

    // if(B < 0)
    //   TBOX_ERROR("B is smaller than 0 "<<gd->mi11<<" "<<gd->mi22<<" "<<gd->mi12<<"\n");

    // if(C < 0)
    //   TBOX_ERROR("C is smaller than 0 "<<gd->mi11<<" "<<gd->mi22<<" "<<gd->mi32<<" gd->mi12"<<" "<<gd->mi23<<"\n");

    // if(D < 0)
    //   TBOX_ERROR("D is smaller than 0 "<<gd->mi11<<" "<<gd->mi22<<" "<<gd->mi23<<"\n");
    
    
    // real_t exx = 1 / sqrt(gd->mi11);
    // real_t eyx = -gd->mi12 / sqrt(B);
    // real_t eyy = -gd->mi11 / sqrt(B);
    // real_t exz = (gd->mi13 * gd->mi22 - gd->mi12 * gd->mi23) / C / D;
    // real_t eyz = (gd->mi11 * gd->mi23 - gd->mi12 * gd->mi13) / C / D;;
    // real_t ezz = (pw2(gd->mi12) - gd->mi11 * gd->mi22) / C / D;

    real_t exx = sqrt(gd->m11);
    real_t exy = gd->m12 / sqrt(gd->m11);
    real_t eyy = sqrt(gd->m22 - exy * exy);
    real_t exz = gd->m13 / sqrt(gd->m11);
    real_t eyz = (gd->m23 - gd->m12 * gd->m13 / gd->m11) / eyy;
    real_t ezz = sqrt(gd->m33 - exz * exz - eyz * eyz);
    
    
    ParticleContainer *pc = new ParticleContainer(*ic);
    int id[PARTICLE_INT_PROPERTIES] = {0};
    double p_info_main[PARTICLE_NUMBER_OF_STATES] = {0};
    double p_info_ep1[PARTICLE_NUMBER_OF_STATES] = {0};
    double p_info_ep2[PARTICLE_NUMBER_OF_STATES] = {0};
    
    p_info_main[0] = p_info_ep1[0] = p_info_ep2[0] = ox;
    p_info_main[1] = p_info_ep1[1] = p_info_ep2[1] = oy;
    p_info_main[2] = p_info_ep1[2] = p_info_ep2[2] = oz;
    int cnt = 0;
    for(int i = 0; i < dirx.size(); i++)
    {

      real_t vhat[3] = {
        exx * dirx[i],
        exy * dirx[i] + eyy * diry[i],
        exz * dirx[i] + eyz * diry[i] + ezz * dirz[i]
      };
      real_t magvhat = vhat[0] * vhat[0] * gd->mi11 + vhat[1] * vhat[1] * gd->mi22
        + vhat[2] * vhat[2] * gd->mi33 + 2.0 * vhat[0] * vhat[1] * gd->mi12
        + 2.0 * vhat[0] * vhat[2] * gd->mi13 + 2.0 * vhat[1] * vhat[2] * gd->mi23;

      vhat[0] /= magvhat, vhat[1] /= magvhat, vhat[2] /= magvhat;
      p_info_main[3] = vhat[0], p_info_main[4] = vhat[1], p_info_main[5] = vhat[2];

      real_t s1[3] = {0, 1, 0};
      real_t s2[3] = {1, 0, 0};
      
      // special case if in x- or y-direction; or x-y plane
      if(std::abs(vhat[2]) <= 1e-10 && std::abs(vhat[1]) <= 1e-10)
      {
        s2[2] = 1.0;
        s2[0] = 0.0;
      }
      else if(std::abs(vhat[2]) <= 1e-10 && std::abs(vhat[0]) <= 1e-10)
      {
        s1[2] = 1.0;
        s1[1] = 0.0;
      }
      else
      {
        // if V is in x-y plane, pick vectors out of that plane
        if(std::abs(vhat[2]) <= 1e-10) { s1[2] = 1.0; s2[2] = 1.0; }
        // if V is in x-z plane, pick vectors out of that plane
        if(std::abs(vhat[1]) <= 1e-10) { s2[1] = 1.0; }
        // if V is in y-z plane, pick vectors out of that plane
        if(std::abs(vhat[0]) <= 1e-10) { s1[0] = 1.0; }
      }

      real_t s1hat[3] = {0}, s2hat[3] = {0};

      real_t s1dotvhat = GEODESIC_DOT_COV_VECTORS(s1, vhat);
      for(int j = 0; j < DIM; j++)
        s1hat[j] = s1[j] - s1dotvhat * vhat[j];

      real_t mags1hat = GEODESIC_DOT_COV_VECTORS(s1hat, s1hat);

      for(int j = 0; j < 3; j++)
        s1hat[j] /= mags1hat;

      real_t s2dotvhat = GEODESIC_DOT_COV_VECTORS(s2, vhat);
      real_t s2dots1hat = GEODESIC_DOT_COV_VECTORS(s2, s1hat);
      for(int j = 0; j < 3; j++)
        s2hat[j] = s2[j] - s2dotvhat * vhat[j] - s2dots1hat * s1hat[j];
      real_t mags2hat = GEODESIC_DOT_COV_VECTORS(s2hat, s2hat);
      for(int j = 0; j < 3; j++)
        s2hat[j] /= mags2hat;

      for(int j = 0; j < 3; j++)
      {
        p_info_ep1[j+3] = vhat[j] + epsilon * s1hat[j];
        p_info_ep2[j+3] = vhat[j] + epsilon * s2hat[j];
      }
      id[0] = p_id++;
      RKParticle *main_p = new RKParticle(p_info_main, id);
      id[0] = p_id++;      
      RKParticle *ep1_p = new RKParticle(p_info_ep1, id);
      id[0] = p_id++;      
      RKParticle *ep2_p = new RKParticle(p_info_ep2, id);

      //      ParticleContainer *pc = new ParticleContainer(*ic);
      cnt++;
      ParticleContainer *pc = pc_pdata->getItem(*ic);
      
      if(pc == NULL)
      {
        pc = new ParticleContainer(*ic);
        pc->addParticle(*main_p);
        pc->addParticle(*ep1_p);
        pc->addParticle(*ep2_p);
        pc_pdata->appendItem(*ic, *pc);
      }
      else
      {
        pc->addParticle(*main_p);
        pc->addParticle(*ep1_p);
        pc->addParticle(*ep2_p);
      }
      // std::cout<<"Inserting ";
      // pc->print();
      // std::cout<<std::flush<<"\n";
      
      
    }
    
  }

  
  void Geodesic::geodesic_ic_set_uniform_rays(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    BSSN *bssn)
  {
    const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
    int N_origins = cosmo_geodesic_db->getInteger("N_origins");
    std::vector<real_t> origins_x, origins_y, origins_z;
    origins_x.resize(N_origins);
    origins_y.resize(N_origins);
    origins_z.resize(N_origins);
    cosmo_geodesic_db->getDoubleArray("origins_x", origins_x.data(), N_origins);
    cosmo_geodesic_db->getDoubleArray("origins_y", origins_y.data(), N_origins);
    cosmo_geodesic_db->getDoubleArray("origins_z", origins_z.data(), N_origins);
    
    // reading direction info from file generated by healpy
    std::ifstream infile("healpix_data_NSIDE_1.dat");
    std::string line;
    std::vector<std::vector<real_t>> healpy_data;
    while(std::getline(infile, line))
    {
      std::istringstream iss(line);
      real_t d;
      std::vector<real_t> v;
      while(iss >> d)
      {
        v.push_back(d);
      }
      healpy_data.push_back(v);
    }

    if(healpy_data.size()!=3)
      TBOX_ERROR("Healpix input has wrong array size! Expecting 3 but having "<<healpy_data.size()<<"\n");
    
    std::vector<std::vector<real_t>> dirs_x, dirs_y, dirs_z;
    dirs_x.resize(N_origins);
    dirs_y.resize(N_origins);
    dirs_z.resize(N_origins);

    for(int i = 0; i < N_origins; i++)
    {
      dirs_x[i] = healpy_data[0];
      dirs_y[i] = healpy_data[1];
      dirs_z[i] = healpy_data[2];
    }
    geodesic_ic_set_bundles(
      hierarchy,
      origins_x, origins_y, origins_z, dirs_x, dirs_y, dirs_z, bssn);
    printAll(hierarchy, pc_idx);
    //    TBOX_ERROR("Here");
  }
  
  void Geodesic::geodesic_ic_set_bundles(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    std::vector<real_t> &origins_x, std::vector<real_t> &origins_y, std::vector<real_t> &origins_z,
    std::vector<std::vector<real_t>> &dirs_x,
    std::vector<std::vector<real_t>> &dirs_y,
    std::vector<std::vector<real_t>> &dirs_z,
    BSSN *bssn)
  {
    const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
    double epsilon = cosmo_geodesic_db->getDoubleWithDefault("epsilon", 0.01);
    std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
        hierarchy->getGridGeometry()));

    geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

    const double * domain_lower = &grid_geometry.getXLower()[0];
    const double * domain_upper = &grid_geometry.getXUpper()[0];

    // setting number of origins
    int n_origins = origins_x.size();
    // setting particle ids
    int p_id = 0;
    
    if(n_origins != dirs_x.size())
      TBOX_ERROR("Size of origins and dirs not match!");

    for(int o_idx = 0; o_idx < n_origins; o_idx++)
    {
      // finding highest patch that covers this origin
      int has_inserted = 0;
      for(int ln = hierarchy->getNumberOfLevels() - 1; ln>=0; ln --)
      {
        std::shared_ptr<hier::PatchLevel> level(
          hierarchy->getPatchLevel(ln));
        hier::PatchLevel::iterator ip(level->begin());

        
        // loop over patches on level
        for (;ip != level->end(); ++ip)
        {
          std::shared_ptr<hier::Patch> patch(*ip);

          const hier::Box& box = patch->getBox();
        
          const int * lower = &box.lower()[0];
          const int * upper = &box.upper()[0];

          std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
              patch->getPatchGeometry()));

          TBOX_ASSERT(patch_geometry);

          // initializing dx
          const double* dx = patch_geometry->getDx();

         
          int i0 = floor((origins_x[o_idx] - domain_lower[0] ) / dx[0] );
          int j0 = floor((origins_y[o_idx] - domain_lower[1] ) / dx[1] );
          int k0 = floor((origins_z[o_idx] - domain_lower[2] ) / dx[2] );
          std::cout<<i0<<" "<<j0<<" "<<k0<<"\n";
          // have found the patch that covers the starting point
          if( i0 >= lower[0] && i0 <= upper[0]
              && j0 >= lower[1] && j0 <= upper[1]
              && k0 >= lower[2] && k0 <= upper[2])
          {
            has_inserted = 1;

            //            for(int d_idx; d_idx < dir_x[o_idx].size(); d_idx++)
            {
              geodesic_ic_transform_inertial_vectors(
                patch, bssn, origins_x[o_idx], origins_y[o_idx], origins_z[o_idx],
                dirs_x[o_idx], dirs_y[o_idx], dirs_z[o_idx], dx, epsilon, p_id);
              p_id+=dirs_x[o_idx].size() * 3; 
            }
            break;
          }

          
        
        }
        // have to barrier here to ensure the continuiaty of p_id
        mpi.Barrier();
        if(ip != level->end())
          break;
      }
    }
    
  }
  
  // integrating a single geodesic in Shschwzchild blackhole 
  void Geodesic::geodesic_ic_Schwarzchild_test(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    std::shared_ptr<tbox::Database> cosmo_geodesic_db)
  {
    const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
    std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
        hierarchy->getGridGeometry()));

    geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

    const double * domain_lower = &grid_geometry.getXLower()[0];
    const double * domain_upper = &grid_geometry.getXUpper()[0];

    double p_info[PARTICLE_NUMBER_OF_STATES];
    
    p_info[0] = cosmo_geodesic_db->getDoubleWithDefault("x0", 0);
    p_info[1] = cosmo_geodesic_db->getDoubleWithDefault("x1", 0);
    p_info[2] = cosmo_geodesic_db->getDoubleWithDefault("x2", 0);

    p_info[3] = cosmo_geodesic_db->getDoubleWithDefault("q0", 0);
    p_info[4] = cosmo_geodesic_db->getDoubleWithDefault("q1", 0);
    p_info[5] = cosmo_geodesic_db->getDoubleWithDefault("q2", 0);

#if EVOLVE_LAMBDA
    p_info[6] = cosmo_geodesic_db->getDoubleWithDefault("lambda_0", 0);
#endif
    
    
    int has_inserted = 0;
    for(int ln = hierarchy->getNumberOfLevels() - 1; ln>=0; ln --)
    {
      std::shared_ptr<hier::PatchLevel> level(
        hierarchy->getPatchLevel(ln));
      hier::PatchLevel::iterator ip(level->begin());
      // loop over patches on level
      for (;ip != level->end(); ++ip)
      {
        std::shared_ptr<hier::Patch> patch(*ip);

        initPData(patch);
        const hier::Box& box = patch->getBox();
        
        const int * lower = &box.lower()[0];
        const int * upper = &box.upper()[0];

        
        std::shared_ptr<pdat::IndexData<ParticleContainer, pdat::CellGeometry>>
          & pc = pc_pdata;
      
        std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
          SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch->getPatchGeometry()));

        TBOX_ASSERT(patch_geometry);

        // initializing dx
        const double* dx = patch_geometry->getDx();

        int i0 = floor((p_info[0] - domain_lower[0] ) / dx[0] );
        int j0 = floor((p_info[1] - domain_lower[1] ) / dx[1] );
        int k0 = floor((p_info[2] - domain_lower[2] ) / dx[2] );

        // have found the patch that covers the starting point
        if( i0 >= lower[0] && i0 <= upper[0]
            && j0 >= lower[1] && j0 <= upper[1]
            && k0 >= lower[2] && k0 <= upper[2])
            {
              hier::Index *ic = new hier::Index(i0, j0, k0);
              ParticleContainer *pc = new ParticleContainer(*ic);
              int id[1] = {0};
              RKParticle *temp_p = new RKParticle(p_info, id);
              pc->addParticle(*temp_p);
              pc_pdata->appendItem(*ic,*pc);
              delete  pc;
              has_inserted = 1;
              break;
            }
        }
      //      if(ip != level->end()) break;
      if (mpi.getSize() > 1) {
        mpi.AllReduce(&has_inserted, 1, MPI_MAX);
      }

      if(has_inserted > 0)
        break;
      }
      

    }

// test paper 1611.09275
void Geodesic::geodesic_ic_face_null_test(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    std::shared_ptr<tbox::Database> cosmo_geodesic_db)
  {
    std::shared_ptr<geom::CartesianGridGeometry> grid_geometry_(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
        hierarchy->getGridGeometry()));

    geom::CartesianGridGeometry& grid_geometry = *grid_geometry_;

    const double * domain_lower = &grid_geometry.getXLower()[0];
    const double * domain_upper = &grid_geometry.getXUpper()[0];

    double p_info[PARTICLE_NUMBER_OF_STATES];
    RKParticle *input_p[6];

    double eps = 0.001;
    
    p_info[0] = p_info[1] = p_info[2] = 0;
    p_info[3] = 1;
    p_info[4] = p_info[5] = 0;

    int id[1] = {0};
    input_p[0] = new RKParticle(p_info, id);

    p_info[4] = eps, p_info[5] = 0, id[0] = 1;
    input_p[1] = new RKParticle(p_info, id);

    p_info[5] = eps, p_info[4] = 0, id[0] = 2;
    input_p[2] = new RKParticle(p_info, id);

    p_info[3] = 1.0/sqrt(2);
    p_info[4] = 1.0/sqrt(2);
    p_info[5] = 0;
    id[0] = 3;
    input_p[3] = new RKParticle(p_info, id);

    p_info[3] = (1.0 - eps) / sqrt(2);
    p_info[4] = (1.0 + eps) / sqrt(2);
    p_info[5] = 0;
    id[0] = 4;
    input_p[4] = new RKParticle(p_info, id);

    p_info[3] = (1.0) / sqrt(2);
    p_info[4] = (1.0) / sqrt(2);
    p_info[5] = eps;
    id[0] = 5;
    input_p[5] = new RKParticle(p_info, id);


    /********Finish creating particles********/

    
    
#if EVOLVE_LAMBDA
    p_info[6] = cosmo_geodesic_db->getDoubleWithDefault("lambda_0", 0);
#endif
    
    
    for(int i = 0; i < 6; i ++)
    {
      for(int ln = hierarchy->getNumberOfLevels() - 1; ln>=0; ln --)
      {

        std::shared_ptr<hier::PatchLevel> level(
          hierarchy->getPatchLevel(ln));
        hier::PatchLevel::iterator ip(level->begin());
        // loop over patches on level
        for (;ip != level->end(); ++ip)
        {
          std::shared_ptr<hier::Patch> patch(*ip);
          initPData(patch);
          const hier::Box& box = patch->getBox();
        
          const int * lower = &box.lower()[0];
          const int * upper = &box.upper()[0];

        
          std::shared_ptr<pdat::IndexData<ParticleContainer, pdat::CellGeometry>>
            & pc = pc_pdata;
      
          std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
              patch->getPatchGeometry()));

          TBOX_ASSERT(patch_geometry);

          // initializing dx
          const double* dx = patch_geometry->getDx();

          int i0 = floor((input_p[i]->x_a[0] - domain_lower[0] ) / dx[0] );
          int j0 = floor((input_p[i]->x_a[1] - domain_lower[1] ) / dx[1] );
          int k0 = floor((input_p[i]->x_a[2] - domain_lower[2] ) / dx[2] );

          std::cout<<i0<<" "<<j0<<" "<<k0<<"\n";
          // have found the patch that covers the starting point
          if( i0 >= lower[0] && i0 <= upper[0]
              && j0 >= lower[1] && j0 <= upper[1]
              && k0 >= lower[2] && k0 <= upper[2])
          {

            hier::Index *ic = new hier::Index(i0, j0, k0);

            ParticleContainer *pc = pc_pdata->getItem(*ic);

            if(pc == NULL)
            {
              ParticleContainer *p = new ParticleContainer(*ic);
              pc_pdata->appendItem((*ic), *p);
              pc = pc_pdata->getItem(*ic);
            }
            pc->addParticle(*input_p[i]);
            break;
          }
        }
        if(ip != level->end()) break;

      }
      

    }

  }
}

