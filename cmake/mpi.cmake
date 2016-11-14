# Check for openmp support
find_package(MPI)
if (MPI_FOUND)
    include_directories( ${MPI_INCLUDE_PATH} )
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OMPI_CFLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OMPI_CXXFLAGS}")
else()
  message(STATUS "MPI support is not available; not compiling with MPI support.")
endif()
