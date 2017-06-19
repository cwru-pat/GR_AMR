#include "static.h"
#include "../../cosmo_includes.h"
#include "../../utils/math.h"


using namespace SAMRAI;

namespace cosmo
{

Static::Static(
  const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const tbox::Dimension& dim_in,
  boost::shared_ptr<tbox::Database> database_in,
  std::ostream* l_stream_in):
  lstream(l_stream_in),
  cosmo_static_db(database_in),
  dim(dim_in)
{
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  boost::shared_ptr<hier::VariableContext> context_active(
    variable_db->getContext("ACTIVE"));

  VAR_INIT(DIFFD);

  DIFFD_a_idx =         
    variable_db->registerVariableAndContext(
      DIFFD,
      context_active, 
      hier::IntVector(dim, STENCIL_ORDER));

  
}
Static::~Static()
{
  
}

void Static::addBSSNSrc(
  BSSN *bssn,   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  for(idx_t ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    const boost::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));
    addBSSNSrc(bssn, level);
  }
}

  
void Static::addBSSNSrc(
  BSSN *bssn, const boost::shared_ptr<hier::PatchLevel> & level)
{
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const boost::shared_ptr<hier::Patch> & patch = *pit;
    addBSSNSrc(bssn, patch);
  }
}

void Static::addBSSNSrc(
  BSSN *bssn, const boost::shared_ptr<hier::Patch> & patch)
{

  bssn->initPData(patch);
  bssn->initMDA(patch);

  DIFFD_a_pdata =                               
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>( 
      patch->getPatchData(DIFFD_a_idx));
  DIFFD_a = pdat::ArrayDataAccess::access<DIM, double>(  
    DIFFD_a_pdata->getArrayData());

  
  const hier::Box& box = bssn->DIFFchi_a_pdata->getGhostBox();
  
  const int * lower = &box.lower()[0];
  const int * upper = &box.upper()[0];

  
  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(  
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch->getPatchGeometry()));

#pragma omp parallel for collapse(2)
  for(int k = lower[2]; k <= upper[2]; k++)
  {
    for(int j = lower[1]; j <= upper[1]; j++)
    {
      for(int i = lower[0]; i <= upper[0]; i++)
      {
        bssn->DIFFr_a(i, j, k) +=
          pw3(bssn->DIFFchi_a(i,j,k) + 1.0) * (DIFFD_a(i,j,k) / (1.0 + bssn->DIFFalpha_a(i,j,k)));
      }
    }
  }
  
}


} // namespace cosmo
