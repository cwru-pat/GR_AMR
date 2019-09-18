# GR_AMR
MPI version of CosmoGRaPH Code (https://cwru-pat.github.io/cosmograph/) plus more features:

 - Fully general relativistic BSSN code with block-structured Adaptive Mesh Refinement support powered by SAMRAI (https://github.com/LLNL/SAMRAI/);
 - Elliptical solver which supports any kind of terms and any number of equations;
 - Support different matter fields, e.g., vacuum, dust fluid and scalar field;
 - AHFinderDirect (https://arxiv.org/pdf/gr-qc/0306056.pdf) is included as apparent horizon finder; 
 - Module that calculates local measurements like spin of black hole is available;
 - General relativistic ray tracing;
 - High performance???

### Software dependencies
 - git
 - hdf5
 - cmake
 - fftw3
 - a compiler with c++11 support 
 - SAMRAI with appropriate version is included as sub-module

## Setting up SAMRAI:

Clone this repository (and SAMRAI)
```
git clone --recursive
cd GR_AMR
```

Make directory to install SAMRAI
```
mkdir obj
cd obj
```

Configure SAMRAI
`sh ../SAMRAI/configure  [--with-hdf5=/parent/path/of/hdf5.h] --with-F77=gfortran`

Build and install SAMRAI
```
make library
make tools
make install
```
Optionally build SAMRAI's documentation:
`make dox`
(documentation will be in obj/docs/samrai-dox)


## Setting up GR_AMR:

Create a build directory for GR_AMR
```
cd ..
mkdir build
cd build
```

Make the program
```
cmake ..
make
```

Run it! (Example: `./cosmo ../input/static_blackhole.input`)

