# GR_AMR
BSSN code bases on SAMRAI, now doing test for solving wave equation

How to make:

git clone

cd GR_AMR

mkdir obj

cd obj

sh ../SAMRAI/configure --with-boost= "...path of boost" --with-hdf5="...path of hdf5"

make library

make tools

make install

cd ..

mkdir build

cd build

cmake ..

make
