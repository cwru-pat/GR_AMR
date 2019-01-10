#ifndef COSMO_STATIC
#define COSMO_STATIC

#include "../../cosmo_includes.h"
#include "../bssn/bssn.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "static_ic.h"

namespace cosmo
{

/** Static matter class **/
class Static
{
  /* Fluid field */
  // just a density variable
  VAR_CREATE(DIFFD);
public:

  RK4_IDX_CREATE(DIFFD,a);
  RK4_PDATA_CREATE(DIFFD,a);
  RK4_MDA_ACCESS_CREATE(DIFFD,a);

  std::ostream* lstream;
  std::shared_ptr<tbox::Database>& cosmo_static_db;
  const tbox::Dimension& dim;
  
  Static(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const tbox::Dimension& dim_in,
    std::shared_ptr<tbox::Database> database_in,
    std::ostream* l_stream_in);
  
  ~Static();
  
  void addBSSNSrc(
    BSSN *bssn,   const std::shared_ptr<hier::PatchHierarchy>& hierarchy);
  void addBSSNSrc(
    BSSN *bssn, const std::shared_ptr<hier::PatchLevel> & level);
  void addBSSNSrc(
    BSSN *bssn, const std::shared_ptr<hier::Patch> & patch);

  
  
};

}

#endif
