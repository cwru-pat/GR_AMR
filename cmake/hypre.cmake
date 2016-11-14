find_library(HYPRE_LIBRARY
  NAMES HYPRE
  HINTS ENV LD_LIBRARY_PATH)

find_path(HYPRE_INCLUDE_DIRS
  NAMES HYPREf.h
  HINTS ENV CPLUS_INCLUDE_PATH)

if(NOT HYPRE_LIBRARY)
  message(FATAL_ERROR "${Red}The HYPRE library and include locations were not found. Please make sure you have HYPRE installed or loaded, and that the library and include directories can be found in CPLUS_INCLUDE_PATH and LD_LIBRARY_PATH environment variables.${ColorReset}")
else()
  set(HYPRE_LIBRARIES "${HYPRE_LIBRARY}")
  if(NOT HYPRE_INCLUDE_DIRS)
    message(STATUS "HYPRE partially found. An include directory was not found for HYPRE in the CPLUS_INCLUDE_PATH environment variable.")
    message(STATUS "If compilation or linking fails, you may need to add the include path to this variable.")
    message(STATUS " HYPRE_LIBRARY: ${HYPRE_LIBRARY}")
  else()
    include_directories("${HYPRE_INCLUDE_DIRS}")
    message(STATUS "Found HYPRE:")
    message(STATUS " HYPRE_LIBRARY: ${HYPRE_LIBRARY}")
    message(STATUS " HYPRE_INCLUDE_DIRS: ${HYPRE_INCLUDE_DIRS}")
  endif()
endif()
