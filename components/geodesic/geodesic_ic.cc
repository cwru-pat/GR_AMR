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
  // integrating a single geodesic in Shschwzchild blackhole 
  void Geodesic::geodesic_ic_Schwarzchild_test(
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
    
    p_info[0] = cosmo_geodesic_db->getDoubleWithDefault("x0", 0);
    p_info[1] = cosmo_geodesic_db->getDoubleWithDefault("x1", 0);
    p_info[2] = cosmo_geodesic_db->getDoubleWithDefault("x2", 0);

    p_info[3] = cosmo_geodesic_db->getDoubleWithDefault("q0", 0);
    p_info[4] = cosmo_geodesic_db->getDoubleWithDefault("q1", 0);
    p_info[5] = cosmo_geodesic_db->getDoubleWithDefault("q2", 0);

#if EVOLVE_LAMBDA
    p_info[6] = cosmo_geodesic_db->getDoubleWithDefault("lambda_0", 0);
#endif
    
    
    
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

        std::cout<<i0<<" "<<j0<<" "<<k0<<"\n";
        // have found the patch that covers the starting point
        if( i0 >= lower[0] && i0 <= upper[0]
            && j0 >= lower[1] && j0 <= upper[1]
            && k0 >= lower[2] && k0 <= upper[2])
            {
              pdat::CellIterator icend(pdat::CellGeometry::end(patch->getBox()));
              //              for (pdat::CellIterator ic(pdat::CellGeometry::begin(patch->getBox()));
              //   ic != icend; ++ic) {
              hier::Index *ic = new hier::Index(i0, j0, k0);
              ParticleContainer *pc = new ParticleContainer(*ic);
              int id[1] = {0};
              RKParticle *temp_p = new RKParticle(p_info, id);
              pc->addParticle(*temp_p);
              //     pc->addParticle(*temp_p);
              pc_pdata->appendItem(*ic,*pc);
              delete  pc;
              
              break;
            }
        }
      if(ip != level->end()) break;

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

