# About this package

This package provides the interface module between aenet [1] and LAMMPS [2], patch of aenet for the LAMMPS library, and Artificial Neural network (ANN) potential parameter file of BCC iron.
This package is distributed under the GNU General Public License, and there is no warranty.
If you have any troubles and questions, feel free to contact package author.

# Installation

0. This manual explains only clean install case using gnu g++, gfortran, and openmpi under linux OS. 
  
1. Clone this package to your appropriate directory.
``` 
git clone https://github.com/HidekiMori-CIT/aenet-lammps.git
```

2. Download (and copy) LAMMPS and aenet package to same directory of aenet-lammps from:  
[LAMMPS](https://lammps.sandia.gov/) :lammps-stable.tar.gz, currently 3Mar20  
[aenet](http://ann.atomistic.net/) :aenet-2.0.3.tar.bz2  

3. Unpack the lammps package and copy USER-AENET/ and aenet/ in aenet-lammps/ to src/ and lib/ in lammps(-3Mar20), respectively, and copy aenet-2.0.3.tar.bz2 to lib/aenet/ in lammps(-3Mar20) directory.
```
tar -xvzf lammps-stable.tar.gz
cp -r ./aenet-lammps/USER-AENET/ ./lammps-3Mar20/src/
cp -r ./aenet-lammps/aenet/ ./lammps-3Mar20/lib/
cp aenet-2.0.3.tar.bz2 lammps-3Mar20/lib/aenet/
```

4. Patch and compile aenet library.
```
cd lammps-3Mar20/lib/aenet/
tar -jvxf aenet-2.0.3.tar.bz2
patch -u -p1 -d aenet-2.0.3/ < aenet_lammps.patch
cd aenet-2.0.3/src/
make -f makefiles/Makefile.gfortran_serial lib
cd ../../
```
note: If you use intel compiler, replace gfortran to ifort.  

5. Check the below two files are created:  
(a) libaenet.a: library of aenet  
(b) Makefile.lammps: linker of aenet and lammps  

6. Compile LAMMPS with aenet module.
```
cd ../../src/
make yes-user-aenet
make mpi
```
note: Other compile option of LAMMPS, please see LAMMPS manual.

# How to run LAMMPS with ANN potential

In ANN/Fe/ directory parameter file and example input file for LAMMPS  
(a) Fe.10tw-10tw.ann: parameter file of ANN potential for BCC iron   
(b) in.aenet_mm: example input file for LAMMPS, lattice optimization of BCC iron  
(c) in.aenet_nve: example input file for LAMMPS, run nve MD of BCC iron  

To check the LAMMPS work properly, for example, in ANN/Fe/ directory
```
../../../lammps-3Mar20/src/lmp_mpi -i in.aenet_mm
```
After job done, please check lattice constant (lx) in log.lammps and log.lammps.g++_mm

In in.aenet_mm(nve) the following two lines activate the ANN potential:
```
pair_style      aenet
pair_coeff      * * v01 Fe 10tw-10tw.ann Fe
```

The pair_style line might be always same.
In pair_coeff line, v01 means to use our original version. 
**_If you chose v00, you can use parameter file from original aenet package._**
The Fe between v01 and 10tw-10tw.ann set element(s).
The term:10tw-10tw.ann is **_file mask_** for name of parameter file, in this case, file name is set as Fe.10tw-10tw.ann.
The final Fe assign elements to atom type.
If you have multi atom type such as atom type 1:Fe(free), 2:Fe(fix), set pair_style line as follow:  
```
pair_coeff      * * v01 Fe 10tw-10tw.ann Fe Fe 
```

# Citing of this package and ANN potential
Please use this bibtex,  
@article{mori2020neural,  
  title={Neural network atomic potential to investigate the dislocation dynamics in bcc iron},  
  author={Mori, Hideki and Ozaki, Taisuke},  
  journal={Physical Review Materials},  
  volume={4},  
  number={4},  
  pages={040601},  
  year={2020},  
  publisher={APS}  
}

# Reference
[1] N. Artrith and A. Urban, Comput. Mater. Sci. 114, 135 (2016).  
[2] S. Plimpton, J. Comput. Phys. 117, 1 (1995).  
 

# Author & contact information
Author: Hideki Mori, College of Industrial Technology, Japan  
E-mail: morih@cit.sangitan.ac.jp

