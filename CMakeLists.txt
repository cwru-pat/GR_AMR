
project(COSMO)
cmake_minimum_required(VERSION 2.8.3)

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

# check openmp
include(cmake/openmp.cmake)

if(EXISTS "/home/chris/.chris_pc_rec_file")

 include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

 set(SAMRAI_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)

 message("Detacting chris's pc, setting corresponding SAMRAI lib.")

elseif(EXISTS "/home/cxt282/.case_cluster_rec_file")

 include_directories(/home/cxt282/SAMRAI/include)

 set(SAMRAI_LIB_DIR /home/cxt282/SAMRAI/lib)

 message("Detacting Case cluster, setting corresponding SAMRAI lib.")

else()

 include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

 set(SAMRAI_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)

 message("NOT detacting any pre-set platform, setting SAMRAI lib to current directory!")

endif(EXISTS "/home/chris/.chris_pc_rec_file")


message(STATUS " SAMRAI_LIBRARIES: ${SAMRAI_LIB_DIR}")
unset(COSMO_SOURCES CACHE)
file(GLOB COSMO_SOURCES cosmo*.cc components/bssn/*.cc utils/*.cc components/IO/*.cc components/statistic/*.cc components/boundaries/*.cc sims/*.cc components/static/*.cc components/scalar/*.cc components/horizon/*.cc ICs/*.cc components/elliptic_solver/*.cc components/horizon/AHFD/*.cc components/horizon/AHFD/driver/BH_diagnostics.cc components/horizon/AHFD/driver/horizon_sequence.cc components/horizon/AHFD/elliptic/*.cc  components/horizon/AHFD/gr/*.cc  components/horizon/AHFD/jtutil/*.cc components/horizon/AHFD/jtutil/*.c components/horizon/AHFD/patch/*.cc components/horizon/AHFD/jtutil/interpolator/common/load.c components/horizon/AHFD/jtutil/interpolator/common/store.c components/horizon/AHFD/jtutil/interpolator/common/evaluate.c components/horizon/AHFD/jtutil/interpolator/Hermite/*.c components/horizon/AHFD/jtutil/interpolator/Hermite/*.c components/horizon/AHFD/jtutil/interpolator/Lagrange-tensor-product/*.c components/horizon/AHFD/jtutil/interpolator/Lagrange-maximum-degree/*.c components/horizon/AHFD/jtutil/interpolator/molecule_posn.c components/horizon/AHFD/jtutil/interpolator/util.c components/horizon/AHFD/jtutil/interpolator/InterpLocalUniform.c  components/horizon/AHFD/sparse-matrix/ilucg/*.f)

add_executable(cosmo ${COSMO_SOURCES})
target_link_libraries(cosmo ${MPI_LIBRARIES} ${HDF5_LIBRARIES} ${FFTW_LIBRARY} ${SAMRAI_LIB_DIR}/libSAMRAI_appu.a ${SAMRAI_LIB_DIR}/libSAMRAI_algs.a ${SAMRAI_LIB_DIR}/libSAMRAI_solv.a ${SAMRAI_LIB_DIR}/libSAMRAI_geom.a   ${SAMRAI_LIB_DIR}/libSAMRAI_mesh.a ${SAMRAI_LIB_DIR}/libSAMRAI_math.a  ${SAMRAI_LIB_DIR}/libSAMRAI_pdat.a ${SAMRAI_LIB_DIR}/libSAMRAI_xfer.a ${SAMRAI_LIB_DIR}/libSAMRAI_hier.a ${SAMRAI_LIB_DIR}/libSAMRAI_tbox.a)
