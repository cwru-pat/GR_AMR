#ifndef COSMO_BSSN_MULTIGRID_BD_FNS
#define COSMO_BSSN_MULTIGRID_BD_FNS

#include <map>

#include "../../utils/Array.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"

using namespace SAMRAI;
namespace cosmo
{

class multigridBdHandler
{
  
private:
  typedef void (multigridBdHandler::*bd_handler_func_t)(double * _array, int nx, int ny, int nz); ///< internal function pointer type

  bd_handler_func_t func;
    
  std::map<std::string, bd_handler_func_t> bd_handler_map;
  void _initBdMaps()
  {
    // Lapse functions
    bd_handler_map["periodic"] = &multigridBdHandler::periodic;
  }

  void periodic(double * _array, int nx, int ny, int nz)
  {

    for(int i = 1; i <= STENCIL_ORDER; i++) //loops for the index inside bd
        for(int j = 0 - STENCIL_ORDER; j < ny  + STENCIL_ORDER; j++)
          for(int k = 0 - STENCIL_ORDER; k < nz  + STENCIL_ORDER; k++)
          {
            _array[B_INDEX(-i, j, k, nx, ny, nz)] =
              _array[B_INDEX(-i+nx , (j+ny)%ny, (k+nz)%nz, nx, ny, nz)];

            _array[B_INDEX(nx+i - 1, j, k, nx, ny, nz)] =
              _array[B_INDEX(i -1, (j+ny)%ny, (k+nz)%nz, nx, ny, nz)];

          }
      for(int i = 0; i < nx ; i++) //loops for the index inside bd
        for(int j = 1; j <= STENCIL_ORDER; j++)
          for(int k = 0 - STENCIL_ORDER; k < nz + STENCIL_ORDER; k++)
          {
            _array[B_INDEX(i, -j, k, nx, ny, nz)] =
              _array[B_INDEX(i, ny-j, (k+nz)%nz, nx, ny, nz)];

            _array[B_INDEX(i, ny+j-1, k, nx, ny, nz)] =
              _array[B_INDEX(i, j-1, (k+nz)%nz, nx, ny, nz)];


          }
      for(int i = 0; i < nx ; i++) //loops for the index inside bd
        for(int j = 0; j < ny ; j++)
          for(int k = 1; k <= STENCIL_ORDER; k++)
          {
            _array[B_INDEX(i, j, nz + k -1, nx, ny, nz)] =
              _array[B_INDEX(i, j, k -1, nx, ny, nz)];
            
            _array[B_INDEX(i, j, -k, nx, ny, nz)] =
              _array[B_INDEX(i, j, nz-k, nx, ny, nz)];
          }
  }
public:
  multigridBdHandler(std::string bd_type)
  {
    _initBdMaps();
    if(bd_handler_map.find(bd_type)!=bd_handler_map.end())
    {
      func = bd_handler_map[bd_type];
    }
    else
    {
       TBOX_ERROR("Multigrid Bd handler cannot find boundary type:"
                 <<bd_type<<"\n");
    }
  }

  
void fillBoundary(double * _array, int nx, int ny, int nz)
{
  (*this.*func)(_array, nx, ny, nz);
}

};

}
#endif
