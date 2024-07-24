# n2p2-4G-D3-Rep-External

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This repository provides a modified version of [n2p2-4G](https://github.com/CompPhysVienna/n2p2/tree/4G-HDNNP-training?tab=readme-ov-file), a software for high-dimensional neural network potentials in computational physics and chemistry. 

Added to the original verion :
- dftd3 grimme correction to the energie
- repulsion potential correction to the energie
- external (any) potential correction to the energie

## Installation

Using git,  Type : 
```console
git clone https://github.com/allouchear/n2p2-4G-D3-Rep-External

```
You can also download the .zip file of n2p2-4G-D3-Rep-External : Click on Code and Download ZIP

The code is interfaced with [LAMMPS](https://www.lammps.org/#gsc.tab=0). To compile LAMMPS for n2p2-4G-D3-Rep-External, see README file in n2p2-4G-D3-Rep-External/src/interface

## How to use it 

See examples/\* folders. 

## Contributors
This software was originally created by Andreas Singraber during his PhD programme at the University of Vienna. Code development started in 2015 and the first release to the public was in 2018.

Contributions are much appreciated and will be recorded here:

 - CabanaMD interface (libnnpif/CabanaMD): Saaketh Desai and Sam Reeve
 - Polynomial symmetry functions (libnnp/Sym(Fnc/Grp)(Base)Comp...): Martin P. Bircher and Andreas Singraber
 - [DFTD3](https://github.com/dftbplus/dftd3-lib) library writen by S. Grimme and coworkers
 - Abdul-Rahman Allouche add dftd3 grimme correction, repulsion potential and  external (any) potential corrections
