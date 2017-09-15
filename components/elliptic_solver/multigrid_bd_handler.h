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

  double L[3];
  double q;
  std::map<std::string, bd_handler_func_t> bd_handler_map;
  void _initBdMaps()
  {
    // Lapse functions
    bd_handler_map["periodic"] = &multigridBdHandler::periodic;
    bd_handler_map["exp_decay"] = &multigridBdHandler::exp_decay;
    bd_handler_map["continuation"] = &multigridBdHandler::continuation;
  }

  void exp_decay(double * _array, int nx, int ny, int nz)
  {
    double x, y, z, r;
    double dx = L[0] / (double) nx, dy = L[1] / (double) ny, dz = L[2] / (double) nz;
    for(int i = 1; i <= STENCIL_ORDER; i++) //loops for the index inside bd
        for(int j = 0 - STENCIL_ORDER; j < ny  + STENCIL_ORDER; j++)
          for(int k = 0 - STENCIL_ORDER; k < nz  + STENCIL_ORDER; k++)
          {
            x = (-(double)i + 0.5) * dx - L[0] / 2.0;
            y = ((double)j + 0.5) * dy - L[1] / 2.0;
            z = ((double)k + 0.5) * dz - L[2] / 2.0;
            r = sqrt(x * x + y * y + z * z);
            _array[B_INDEX(-i, j, k, nx, ny, nz)] = 1;

            x = ((double)(nx+i-1) + 0.5) * dx - L[0] / 2.0;
            y = ((double)j + 0.5) * dy - L[1] / 2.0;
            z = ((double)k + 0.5) * dz - L[2] / 2.0;
            r = sqrt(x * x + y * y + z * z);

            _array[B_INDEX(nx+i - 1, j, k, nx, ny, nz)] = 1;

          }
      for(int i = 0; i < nx ; i++) //loops for the index inside bd
        for(int j = 1; j <= STENCIL_ORDER; j++)
          for(int k = 0 - STENCIL_ORDER; k < nz + STENCIL_ORDER; k++)
          {
            x = ((double)i + 0.5) * dx - L[0] / 2.0;
            y = (-(double)j + 0.5) * dy - L[1] / 2.0;
            z = ((double)k + 0.5) * dz - L[2] / 2.0;
            r = sqrt(x * x + y * y + z * z);
            _array[B_INDEX(i, -j, k, nx, ny, nz)] = 1;

            x = ((double)i + 0.5) * dx - L[0] / 2.0;
            y = ((double)(ny + j -1) + 0.5) * dy - L[1] / 2.0;
            z = ((double)k + 0.5) * dz - L[2] / 2.0;
            r = sqrt(x * x + y * y + z * z);
            _array[B_INDEX(i, ny+j-1, k, nx, ny, nz)] = 1;

          }
      for(int i = 0; i < nx ; i++) //loops for the index inside bd
        for(int j = 0; j < ny ; j++)
          for(int k = 1; k <= STENCIL_ORDER; k++)
          {
            x = ((double)i + 0.5) * dx - L[0] / 2.0;
            y = ((double)j + 0.5) * dy - L[1] / 2.0;
            z = ((double)(nz + k -1) + 0.5) * dz - L[2] / 2.0;
            r = sqrt(x * x + y * y + z * z);
            _array[B_INDEX(i, j, nz + k -1, nx, ny, nz)] = 1;

            x = ((double)i + 0.5) * dx - L[0] / 2.0;
            y = ((double)j + 0.5) * dy - L[1] / 2.0;
            z = (-(double)k + 0.5) * dz - L[2] / 2.0;
            r = sqrt(x * x + y * y + z * z);
            _array[B_INDEX(i, j, -k, nx, ny, nz)] = 1;
          }
    
  }

  void continuation(double * _array, int nx, int ny, int nz)
  {
    double x, y, z, r;
    double dx = L[0] / (double) nx, dy = L[1] / (double) ny, dz = L[2] / (double) nz;
    for(int i = 1; i <= STENCIL_ORDER; i++) //loops for the index inside bd
        for(int j = 0 - STENCIL_ORDER; j < ny  + STENCIL_ORDER; j++)
          for(int k = 0 - STENCIL_ORDER; k < nz  + STENCIL_ORDER; k++)
          {
            _array[B_INDEX(-i, j, k, nx, ny, nz)] =
              3.0 * (_array[B_INDEX(-i+1, j, k, nx, ny, nz)]
                     - _array[B_INDEX(-i+2, j, k, nx, ny, nz)])
              + _array[B_INDEX(-i+3, j, k, nx, ny, nz)];

            _array[B_INDEX(nx+i - 1, j, k, nx, ny, nz)] =
              3.0 * (_array[B_INDEX(nx+i-2, j, k, nx, ny, nz)]
                     - _array[B_INDEX(nx+i-3, j, k, nx, ny, nz)])
              + _array[B_INDEX(nx+i-4, j, k, nx, ny, nz)];


          }
      for(int i = 0; i < nx ; i++) //loops for the index inside bd
        for(int j = 1; j <= STENCIL_ORDER; j++)
          for(int k = 0 - STENCIL_ORDER; k < nz + STENCIL_ORDER; k++)
          {
            _array[B_INDEX(i, -j, k, nx, ny, nz)] =
              3.0 * (_array[B_INDEX(i, -j+1, k, nx, ny, nz)]
                     - _array[B_INDEX(i, -j+2, k, nx, ny, nz)])
              + _array[B_INDEX(i, -j+3, k, nx, ny, nz)];

            _array[B_INDEX(i, ny+j-1, k, nx, ny, nz)] =
              3.0 * (_array[B_INDEX(i, ny+j-2, k, nx, ny, nz)]
                     - _array[B_INDEX(i, ny+j-3, k, nx, ny, nz)])
              + _array[B_INDEX(i, ny+j-4, k, nx, ny, nz)];


          }
      for(int i = 0; i < nx ; i++) //loops for the index inside bd
        for(int j = 0; j < ny ; j++)
          for(int k = 1; k <= STENCIL_ORDER; k++)
          {
            _array[B_INDEX(i, j, -k, nx, ny, nz)] =
              3.0 * (_array[B_INDEX(i, j, -k+1, nx, ny, nz)]
                     - _array[B_INDEX(i, j, -k+2, nx, ny, nz)])
              + _array[B_INDEX(i, j, -k+3, nx, ny, nz)];

            _array[B_INDEX(i, j, nz+k-1, nx, ny, nz)] =
              3.0 * (_array[B_INDEX(i, j, nz+k-2, nx, ny, nz)]
                     - _array[B_INDEX(i, j, nz+k-3, nx, ny, nz)])
              + _array[B_INDEX(i, j, nz+k-4, nx, ny, nz)];

          }
    
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
 multigridBdHandler(std::string bd_type, double L_in[3], double q_in)
  {
    _initBdMaps();
    L[0] = L_in[0], L[1] = L_in[1], L[2] = L_in[2];
    q = q_in;
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
