# Use n2p2-4G-D3-Rep-External via LAMMPS

## Set environement variables
Edit env.sh to set paths to your configuration (openmi,...)

## Compile LAMMPS with nnp
```console
cp env.sh <path-to-LAMMPS>/
cp xcompWithNNP <path-to-LAMMPS>/
cd <path-to-LAMMPS>/src
./xcompWithNNP
```
