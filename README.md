# GR_AMR
BSSN code bases on SAMRAI, now doing test for solving wave equation

How to make:

git clone

cd GR_AMR

mkdir obj

cd obj

sh ../SAMRAI/configure --with-boost= "...path of boost" --with-hdf5="...path of hdf5" --with-F77=gfortran

make library

make tools

make install

cd ..

mkdir build

cd build


(for cluster only)
module load gcc/4.9.3
module load hdf5/1.8.15
module load boost
module load depends
module load cmake



cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc ..

make


