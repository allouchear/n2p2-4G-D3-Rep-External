// n2p2 - A neural network potential package
// Copyright (C) 2018 Andreas Singraber (University of Vienna)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// Created by A.R. Allouche
// Refs :
// Eckstein, W., Hackel, S., Heinemann, D. et al. Influence of the interaction potential on simulated sputtering and reflection data. Z Phys D - Atoms, Molecules and Clusters 24, 171–176 (1992). https://doi.org/10.1007/BF01426703
// https://docs.lammps.org/pair_gromacs.html for S(r) function

#include "Repulsion.h"
#include "utility.h"
#include <cstdlib>   // atof, atoi
#include <iostream> 
#include <vector> 


using namespace std;
using namespace nnp;

 void Repulsion::setnormalization(double mEnergy, double cEnergy, double cLength, double cCharge)
{
    	normalize = true;
	meanEnergy = mEnergy;
	convEnergy = cEnergy;
	convLength = cLength;
	convCharge = cCharge;
}
void Repulsion::setcontoatomicunit(double cEnergy, double cLength)
{
    	// Convert From Phys to Atomic unit
	convToHartree = cEnergy;
	convToBohr = cLength;
}

void Repulsion::calculateEnergy(Structure& structure) const // structure in normalized unit
{
	if(type==RT_UNKNOWN)
	{
		structure.energyRep = 0.0;
		return;
	}
	//cout<<"structure.hasRepE="<<structure.hasRepE<<endl;
	if(structure.hasRepE) return;

	vector<double> eners(structure.atoms.size());
	if(normalize) structure.toPhysicalUnits(meanEnergy, convEnergy, convLength, convCharge);
	// Loop over all atoms, center atom i (ai).
#ifdef _OPENMP
	#pragma omp parallel
	{
	#pragma omp for
#endif
	for (size_t i = 0; i < structure.atoms.size(); ++i)
	{
		// Set pointer to atom.
		Atom &ai = structure.atoms.at(i);
		eners[i] = energy(structure, ai);
	}
#ifdef _OPENMP
	}
#endif
	structure.energyRep = 0.0;
	for (size_t i = 0; i < structure.atoms.size(); ++i)
		structure.energyRep += eners[i];
	structure.energyRep /= 2;
	structure.hasRepE=true;
	if(normalize) structure.toNormalizedUnits(meanEnergy,  convEnergy,  convLength,  convCharge);
}
void Repulsion::calculateForces(Structure& structure) const // structure in normalized  unit
{
	if(type==RT_UNKNOWN)
	{
		structure.energyRep = 0.0;
		int nAtoms = structure.numAtoms;
		for(int i=0;i<3;i++)
		for(int ia=0;ia<nAtoms;ia++)
			structure.atoms[ia].fRep[i] = 0.0;
		return;
	}
	if(structure.hasRepF) return;
	if(normalize) structure.toPhysicalUnits(meanEnergy, convEnergy, convLength, convCharge);
	// Loop over all atoms, center atom i (ai).
#ifdef _OPENMP
	#pragma omp parallel
	{
	#pragma omp for
#endif
	for (size_t i = 0; i < structure.atoms.size(); ++i)
	{
		// Set pointer to atom.
		Atom &ai = structure.atoms.at(i);
		computeForce(structure, ai);
	}
#ifdef _OPENMP
	}
#endif
	structure.hasRepF=true;
	if(normalize) structure.toNormalizedUnits(meanEnergy,  convEnergy,  convLength,  convCharge);
	calculateEnergy(structure);
}
void Repulsion::printForces(Structure& structure, Log log) const
{
    	// Print forces for control
	log << "\n";
	log << "Repulsive potential forces:\n";
	for (size_t i = 0; i < structure.atoms.size(); ++i)
	{
		Atom &ai = structure.atoms.at(i);
		std::string s = structure.elementMap.symbol(ai.element);
        	log << strpr("%10d %2s %16.8E %16.8E %16.8E\n",
                                i+1,
                                s.c_str(),
                                ai.fRep[0],
                                ai.fRep[1],
                                ai.fRep[2]
                                );
	}
	log << "-----------------------------------------"
               "--------------------------------------\n";
	log << "Reppulsion potential energy : ";
	log << strpr("%16.8E\n", structure.energyRep);
	log << "========================================="
               "======================================\n";
}
