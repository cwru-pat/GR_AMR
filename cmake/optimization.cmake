
SET(COSMO_DEBUG FALSE CACHE STRING "Turn on debug mode (turn off optimizations, turn on valgrind")
SET(COSMO_PROFILE FALSE CACHE STRING "Turn on profiling (gprof flag)")
SET(COSMO_STATIC FALSE CACHE STRING "Static build (-static flag)")
SET(COSMO_VERBOSE_COMPILER FALSE CACHE STRING "Verbose compiler (use, eg, -ftree-vectorizer-verbose=2)")


if(COSMO_DEBUG)
  # add -g for valgrind
  message(STATUS "${Cyan} Debug build enabled (-g flag).${ColorReset}")
  set(PROFILING     "-g")
  set(OPT_LEVEL     "-O1")
  set(CC_OPTS       "")
else()
  set(PROFILING     "")
  set(OPT_LEVEL     "-O3")
  set(C_OPT_LEVEL     "-O3")
  set(F_OPT_LEVEL     "-O3")  

  # try to use some GNU compiler special options
  if(CMAKE_COMPILER_IS_GNUCXX)
    set(OPT_LEVEL     "${OPT_LEVEL} -ffast-math -march=native -fext-numeric-literals -Wno-reorder")
    set(C_OPT_LEVEL     "${C_OPT_LEVEL} -ffast-math -march=native")
    set(F_OPT_LEVEL     "${F_OPT_LEVEL} -ffast-math -march=native")   
  # try to use some Intel compiler special options
  # must be using intel here, where else is can be?
  else()
    set(OPT_LEVEL     "${OPT_LEVEL} -march=native -no-prec-div -xhost -qopenmp -parallel -qopenmp-simd")
    set(F_OPT_LEVEL     "${F_OPT_LEVEL} -march=native -no-prec-div -xhost -qopenmp -parallel -qopenmp-simd")
    set(C_OPT_LEVEL     "${C_OPT_LEVEL} -march=native -no-prec-div -xhost -qopenmp -parallel -qopenmp-simd")
  endif()
endif()

# Statically linked build?
if(COSMO_STATIC)
  message(STATUS "${Cyan} Static build enabled (-static flag).${ColorReset}")
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

# add -pg for gprof
if(COSMO_PROFILE)
  message(STATUS "${Cyan} Profiling enabled (-pg flag).${ColorReset}")
  set(PROFILING "${PROFILING} -pg")
endif()

if(COSMO_VERBOSE_COMPILER)
  if(CMAKE_COMPILER_IS_GNUCXX)
    set(CC_OPTS     "${CC_OPTS} -ftree-vectorizer-verbose=2")
  endif()
endif()

set(WARNINGS          "-pedantic -Wall")
set(CMAKE_CXX_FLAGS   "${CC_OPTS} ${OPT_LEVEL} ${WARNINGS} ${PROFILING}")
set(CMAKE_C_FLAGS   "${C_OPT_LEVEL}  ${PROFILING}")
set(CMAKE_FORTRAN_FLAGS   " ${F_OPT_LEVEL} ${WARNINGS} ${PROFILING}")
set(CMAKE_EXE_LINKER_FLAGS  "${PROFILING}")

unset(COSMO_DEBUG CACHE)
unset(COSMO_STATIC CACHE)
unset(COSMO_PROFILE CACHE)
unset(COSMO_VERBOSE_COMPILER CACHE)
