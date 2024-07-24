#!/bin/bash

source /home/theochem/allouche/Softwares/dftbplus-23.1.x86_64-linux/env.sh

fileIn=$1
fileOut=$2

fham=${fileIn%.inp}_ham.txt

function createHamiltonian()
{
	printf "%s\n" "Hamiltonian = xTB {"
	printf "%s\n" "  Method=\"GFN1-xTB\""
	printf "%s\n" "  kPointsAndWeights = SuperCellFolding {"
	printf "%s\n" "       2   0   0"
	printf "%s\n" "       0   2   0"
	printf "%s\n" "       0   0   2"
	printf "%s\n" "       0.5 0.5 0.5"
	printf "%s\n" "  }"
	printf "%s\n" "}"
}
createHamiltonian > $fham

source /home/theochem/allouche/shell/tensorFlowEnv
python3 $N2P2DIR/bin/runDFTBP.py $fileIn $fileOut $fham

