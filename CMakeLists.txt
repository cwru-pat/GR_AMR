
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

# check mpi
include(cmake/mpi.cmake)

# check hdf5
include(cmake/hdf5.cmake)

# check compiler support
include(cmake/compiler.cmake)

# check fftw
include(cmake/fftw.cmake)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(SAMRAI_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)

message(STATUS " SAMRAI_LIBRARIES: ${SAMRAI_LIB_DIR}")
unset(COSMO_SOURCES CACHE)
file(GLOB COSMO_SOURCES cosmo*.cc components/bssn/*.cc utils/*.cc components/IO/*.cc components/boundaries/*.cc sims/*.cc components/static/*.cc ICs/*.cc)

add_executable(cosmo ${COSMO_SOURCES})
target_link_libraries(cosmo ${MPI_LIBRARIES} ${HDF5_LIBRARIES} ${FFTW_LIBRARY} ${SAMRAI_LIB_DIR}/libSAMRAI_appu.a ${SAMRAI_LIB_DIR}/libSAMRAI_algs.a ${SAMRAI_LIB_DIR}/libSAMRAI_solv.a ${SAMRAI_LIB_DIR}/libSAMRAI_geom.a   ${SAMRAI_LIB_DIR}/libSAMRAI_mesh.a ${SAMRAI_LIB_DIR}/libSAMRAI_math.a  ${SAMRAI_LIB_DIR}/libSAMRAI_pdat.a ${SAMRAI_LIB_DIR}/libSAMRAI_xfer.a ${SAMRAI_LIB_DIR}/libSAMRAI_hier.a ${SAMRAI_LIB_DIR}/libSAMRAI_tbox.a)
