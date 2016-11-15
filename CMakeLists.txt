
project(COSMO)
cmake_minimum_required(VERSION 2.6.0)

# coloring for messages
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColorReset  "${Esc}[m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Cyan        "${Esc}[36m")
endif()

# optimization and debugging
include(cmake/optimization.cmake)

# check for OpenMP
include(cmake/openmp.cmake)

# check compiler support
include(cmake/mpi.cmake)

include(cmake/hdf5.cmake)

#include(cmake/hypre.cmake)

# check compiler support
include(cmake/compiler.cmake)


#include(cmake/SAMRAI.cmake)
include_directories(/home/cxt282/GR_AMR/include)





set(SAMRAI_LIB_DIR /home/cxt282/GR_AMR/lib)

message(STATUS " SAMRAI_LIBRARIES: ${SAMRAI_LIB_DIR}")
unset(COSMO_SOURCES CACHE)
file(GLOB COSMO_SOURCES *.C *.h )
add_executable(cosmo ${COSMO_SOURCES} )
target_link_libraries(cosmo  ${MPI_LIBRARIES} ${HDF5_LIBRARIES} ${SAMRAI_LIB_DIR}/libSAMRAI_hier.a ${SAMRAI_LIB_DIR}/libSAMRAI_algs.a ${SAMRAI_LIB_DIR}/libSAMRAI_appu.a ${SAMRAI_LIB_DIR}/libSAMRAI_geom.a  ${SAMRAI_LIB_DIR}/libSAMRAI_math.a ${SAMRAI_LIB_DIR}/libSAMRAI_mesh.a ${SAMRAI_LIB_DIR}/libSAMRAI_pdat.a ${SAMRAI_LIB_DIR}/libSAMRAI_solv.a ${SAMRAI_LIB_DIR}/libSAMRAI_tbox.a ${SAMRAI_LIB_DIR}/libSAMRAI_xfer.a)




