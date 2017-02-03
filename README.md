# GR_AMR
BSSN code bases on SAMRAI, now doing test for solving wave equation

In a cluster environment, the following are needed:
 - module load gcc/4.9.3
 - module load hdf5/1.8.15
 - module load boost
 - module load depends
 - module load cmake

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
`sh ../SAMRAI/configure --with-boost[=/parent/path/of/boost/include] [--with-hdf5=/parent/path/of/hdf5.h] --with-F77=gfortran`

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
cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc ..
make
```

Run it! (Example: `./cosmo ../input/static_blackhole.input`)

