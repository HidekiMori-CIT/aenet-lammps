# aenet-lammps
# About This package

This package contains the neural network (ANN) atomic potential module implemented in LAMMPS, patch of \ae net for the LAMMPS library, and parameter file of iron.
This package is distributed under the GNU General Public License,and there is no warranty.
If you have any troubles and questions, feel free to contact package author:Hideki Mori (morih@cit.sangitan.ac.jp)

# Installation

0. This manual is only case of gnu g++, gfortran, and openmpi under linux OS.
If you use intel compiler, replace gfortran to ifort

1. Clone this package to your appropriate directory
``` 
git clone https://github.com/HidekiMori-CIT/aenet-lammps.git
```

2. Download (and copy) LAMMPS and aenet package to current directory from:  
[LAMMPS](https://lammps.sandia.gov/) :lammps-stable.tar.gz, currently 3Mar20  
[aenet](http://ann.atomistic.net/) :aenet-2.0.3.tar.bz2  

3. Unpack the lammps package and copy this package and aenet-2.0.3.tar.bz2 to src and lib in lammps(-3Mar20) directory.
```
tar -xvzf lammps-stable.tar.gz
cp -r ./aenet-lammps/USER-AENET/ ./lammps-3Mar20/src/
cp -r ./aenet-lammps/aenet/ ./lammps-3Mar20/lib/
cp aenet-2.0.3.tar.bz2 lammps-3Mar20/lib/aenet/
```

4. Patch and compile aenet library
```
cd lammps-3Mar20/lib/aenet
tar -jvxf aenet-2.0.3.tar.bz2
patch -u -p1 -d aenet-2.0.3/ < aenet_lammps.patch
cd aenet-2.0.3/src/
make -f makefiles/Makefile.gfortran_serial lib
cd ../../
```

5. Check the below two files are created:  
library of aenet: libaenet.a  
linker of aenet and lammps: Makefile.lammps  

6. Compile LAMMPS with aenet module
```
cd ../../src/
make yes-user-aenet
make mpi
```
note: other compile option of LAMMPS, please see LAMMPS manual.

# How to run LAMMPS with ANN potential

In ANN/Fe/ directory parameter file and example input file for LAMMPS
Parameter file of ANN potential for BCC iron : Fe.10tw-10tw.ann
Example input file for LAMMPS (1): in.aenet_mm, lattice optimization of BCC iron
Example input file for LAMMPS (2): in.aenet_nve, run nve MD of BCC iron

To check the LAMMPS work properly, for example, in ANN/Fe/ directory
```
../../../lammps-3Mar20/src/lmp_mpi -i in.aenet_mm
```
After job done, please check lattice constant (lx) in log.lammps and log.lammps.g++_mm

In in.aenet_mm(nve) the following two line activate the ANN potential:
pair_style      aenet
pair_coeff      * * v01 Fe 10tw-10tw.ann Fe

The pair_style line might be always same.
In pair_coeff line, v01 mean use our original version. 
If you chose v00, you can use parameter file from original \ae net package.
The Fe betwenn v01 and 10tw-10tw.ann set element(s).
The term:10tw-10tw.ann set name of parameter file, in this case, file name is set as Fe.10tw-10tw.ann.
The final Fe assign elements to atom type.
If you have multi atom type such as atom type 1:Fe(free), 2:Fe(fix), set pair_style line as follow:  
```
pair_coeff      * * v01 Fe 10tw-10tw.ann Fe Fe 
```
