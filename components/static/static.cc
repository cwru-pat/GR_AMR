#include "static.h"
#include "../../cosmo_includes.h"
#include "../../utils/math.h"


using namespace SAMRAI;

namespace cosmo
{

Static::Static(
  const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
  const tbox::Dimension& dim_in,
  std::shared_ptr<tbox::Database> database_in,
  std::ostream* l_stream_in):
  lstream(l_stream_in),
  cosmo_static_db(database_in),
  dim(dim_in)
{
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  std::shared_ptr<hier::VariableContext> context_active(
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
  BSSN *bssn,   const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
  for(idx_t ln = 0; ln < hierarchy->getNumberOfLevels(); ln++)
  {
    const std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(ln));
    addBSSNSrc(bssn, level);
  }
}

  
void Static::addBSSNSrc(
  BSSN *bssn, const std::shared_ptr<hier::PatchLevel> & level)
{
  for( hier::PatchLevel::iterator pit(level->begin());
       pit != level->end(); ++pit)
  {
    const std::shared_ptr<hier::Patch> & patch = *pit;
    addBSSNSrc(bssn, patch);
  }
}

void Static::addBSSNSrc(
  BSSN *bssn, const std::shared_ptr<hier::Patch> & patch)
{

  
}


} // namespace cosmo
