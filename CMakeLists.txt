
project(COSMO)
cmake_minimum_required(VERSION 3.10)

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
#include(cmake/openmp.cmake)

# check mpi
#include(cmake/mpi.cmake)

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

 include_directories(/scratch1/fs1/jmertens/SAMRAI/obj/include)
 #include_directories(/usr/include/RAJA)
 include_directories(/scratch1/fs1/jmertens/RAJA-v0.12.1/obj/include/RAJA/)
 include_directories(/usr/include/umpire)   

 set(SAMRAI_LIB_DIR /scratch1/fs1/jmertens/SAMRAI/obj/lib)

 message("NOT detacting any pre-set platform, setting SAMRAI lib to current directory!")

endif(EXISTS "/home/chris/.chris_pc_rec_file")


if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
 set(CMAKE_CUDA_ARCHITECTURES 70 75)
endif()
message(STATUS "CUDA architectures set to ${CMAKE_CUDA_ARCHITECTURES}")

set(BLT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/blt" CACHE PATH "")
if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/blt/SetupBLT.cmake)
         message(FATAL_ERROR "The BLT submodule is not present. \
                  If in the git repository, please run the following commands:\n \
                  git submodule init\n \
                  git submodule update")
endif()

set(ENABLE_CUDA On CACHE Bool "Enable CUDA")
set(ENABLE_MPI On CACHE Bool "Enable MPI.")
set(BLT_CXX_STD "c++11" CACHE STRING "Version of C++ standard")
set(ENABLE_FORTRAN On CACHE Bool "Enable Fortran")
set(ENABLE_OPENMP On CACHE Bool "Enable OpenMP")


set(CMAKE_CUDA_FLAGS "--expt-extended-lambda")


include(${CMAKE_CURRENT_SOURCE_DIR}/blt/SetupBLT.cmake)


message(STATUS " SAMRAI_LIBRARIES: ${SAMRAI_LIB_DIR}")
unset(COSMO_SOURCES CACHE)
file(GLOB COSMO_SOURCES cosmo*.cc components/bssn/*.cc utils/*.cc components/IO/*.cc components/statistic/*.cc components/boundaries/*.cc sims/*.cc components/static/*.cc components/scalar/*.cc components/dust_fluid/*.cc components/horizon/*.cc ICs/*.cc components/geodesic/geodesic*.cc components/elliptic_solver/*.cc components/horizon/AHFD/*.cc components/horizon/AHFD/driver/BH_diagnostics.cc components/horizon/AHFD/driver/horizon_sequence.cc components/horizon/AHFD/elliptic/*.cc  components/horizon/AHFD/gr/*.cc  components/horizon/AHFD/jtutil/*.cc components/horizon/AHFD/jtutil/*.c components/horizon/AHFD/patch/*.cc components/horizon/AHFD/jtutil/interpolator/common/load.c components/horizon/AHFD/jtutil/interpolator/common/store.c components/horizon/AHFD/jtutil/interpolator/common/evaluate.c components/horizon/AHFD/jtutil/interpolator/Hermite/*.c components/horizon/AHFD/jtutil/interpolator/Hermite/*.c components/horizon/AHFD/jtutil/interpolator/Lagrange-tensor-product/*.c components/horizon/AHFD/jtutil/interpolator/Lagrange-maximum-degree/*.c components/horizon/AHFD/jtutil/interpolator/molecule_posn.c components/horizon/AHFD/jtutil/interpolator/util.c components/horizon/AHFD/jtutil/interpolator/InterpLocalUniform.c  components/horizon/AHFD/sparse-matrix/ilucg/*.f)


#string(REGEX MATCHALL "\\..*\\.c\\/" temp ${COSMO_SOURCES} )
#string(REGEX MATCHALL "(\\..+\\.c)" temp "1d.cube.order2.c/scratch1/fs1/jmertens/GR_AMR/components/horizon/AHFD/jtutil/interpolator/Hermite/2d.cube.order3.c" )

list(APPEND BLT_C_FILE_EXTS ".cube.order2.c")
list(APPEND BLT_C_FILE_EXTS ".cube.order3.c")
list(APPEND BLT_C_FILE_EXTS ".cube.order4.c")
list(APPEND BLT_C_FILE_EXTS ".cube.order1.smooth0.c")
list(APPEND BLT_C_FILE_EXTS ".cube.order2.smooth0.c")
list(APPEND BLT_C_FILE_EXTS ".cube.order3.smooth0.c")
list(APPEND BLT_C_FILE_EXTS ".cube.order4.smooth0.c")
list(APPEND BLT_C_FILE_EXTS ".cube.order5.smooth0.c")
list(APPEND BLT_C_FILE_EXTS ".cube.order6.smooth0.c")

set (convdiff_sources ${COSMO_SOURCES})

set(SAMRAI_LIBRARIES ${MPI_LIBRARIES} ${HDF5_LIBRARIES} ${FFTW_LIBRARY} ${SAMRAI_LIB_DIR}/libSAMRAI_appu.a ${SAMRAI_LIB_DIR}/libSAMRAI_algs.a ${SAMRAI_LIB_DIR}/libSAMRAI_solv.a ${SAMRAI_LIB_DIR}/libSAMRAI_geom.a   ${SAMRAI_LIB_DIR}/libSAMRAI_mesh.a ${SAMRAI_LIB_DIR}/libSAMRAI_math.a  ${SAMRAI_LIB_DIR}/libSAMRAI_pdat.a ${SAMRAI_LIB_DIR}/libSAMRAI_xfer.a ${SAMRAI_LIB_DIR}/libSAMRAI_hier.a ${SAMRAI_LIB_DIR}/libSAMRAI_tbox.a
)

#set(RAJA_LIBRARIES /usr/lib/libRAJA.a)
set(RAJA_LIBRARIES /scratch1/fs1/jmertens/RAJA-v0.12.1/obj/lib/libRAJA.a)
set(UMPIRE_LIBRARIES /usr/lib/libumpire.a)


#set(convdiff_depends testlib)
set(convdiff_depends ${convdiff_depends} ${SAMRAI_LIBRARIES})
set(convdiff_depends ${convdiff_depends} ${RAJA_LIBRARIES})
set(convdiff_depends ${convdiff_depends} ${UMPIRE_LIBRARIES})
set(convdiff_depends ${convdiff_depends} openmp)
set(convdiff_depends ${convdiff_depends} mpi)
set(convdiff_depends ${convdiff_depends} cuda)


blt_add_executable(
  NAME convdiff
  SOURCES ${convdiff_sources}
  DEPENDS_ON ${convdiff_depends})


# add_executable(cosmo ${COSMO_SOURCES})
# target_link_libraries(cosmo ${MPI_LIBRARIES} ${HDF5_LIBRARIES} ${FFTW_LIBRARY} ${SAMRAI_LIB_DIR}/libSAMRAI_appu.a ${SAMRAI_LIB_DIR}/libSAMRAI_algs.a ${SAMRAI_LIB_DIR}/libSAMRAI_solv.a ${SAMRAI_LIB_DIR}/libSAMRAI_geom.a   ${SAMRAI_LIB_DIR}/libSAMRAI_mesh.a ${SAMRAI_LIB_DIR}/libSAMRAI_math.a  ${SAMRAI_LIB_DIR}/libSAMRAI_pdat.a ${SAMRAI_LIB_DIR}/libSAMRAI_xfer.a ${SAMRAI_LIB_DIR}/libSAMRAI_hier.a ${SAMRAI_LIB_DIR}/libSAMRAI_tbox.a)
