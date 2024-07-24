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
#include "External.h"
#include "utility.h"
#include <cstdlib>   // atof, atoi
#include <iostream> 
#include <vector> 
#include <string> 
#include <iomanip> 
#include <sstream>


using namespace std;
using namespace nnp;

std::string addInt(std::string s, int v)
{
	std::stringstream stream;
	stream<<v;
	std::string sv;
	stream>>sv;
	return s+sv;
}
 void External::setnormalization(double mEnergy, double cEnergy, double cLength, double cCharge)
{
    	normalize = true;
	meanEnergy = mEnergy;
	convEnergy = cEnergy;
	convLength = cLength;
	convCharge = cCharge;
}
void External::setcontoatomicunit(double cEnergy, double cLength)
{
    	// Convert From Phys to Atomic unit
	convToHartree = cEnergy;
	convToBohr = cLength;
}

std::string External::createInputFile(Structure& structure, int type) const// type = 0 => energy only, type = 1 => energy + forces
{
	std::string fileName = addInt(prefixFileName,  structure.index);
	fileName += ".inp";
	ofstream f(fileName);
	int mulplicity=1; // To change when n2p2 will support multiplicity
	f<<type<<" "<<int(structure.charge)<<" "<<mulplicity<<endl;
	f.setf(ios::scientific);
	for (size_t i = 0; i < structure.atoms.size(); ++i)
	{
		// Set pointer to atom.
		Atom &ai = structure.atoms.at(i);
		std::string s = structure.elementMap.symbol(ai.element);
		f<<s;
		for(int j=0;j<3;j++)
			f<<" "<<setprecision(20)<<ai.r[j]*convToBohr;
		f<<endl;
	}
	if(structure.isPeriodic)
	{
		for(int i=0;i<3;i++)
		{
			f<<"Tv";
			for(int j=0;j<3;j++)
				f<<" "<<setprecision(20)<<structure.box[i][j]*convToBohr;
			f<<endl;
		}
	}
	f.close();
	//std::string  c = std::string("cat ")+fileName; system(c.c_str());
	return fileName;
}
void External::getResults(Structure& structure, int type) const
{
	std::string fileName = addInt(prefixFileName,  structure.index);
	fileName += ".out";
	ifstream f(fileName);
	if(f.fail())
	{
		cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		cerr<<"I cannot open "<<fileName<<endl;
		cerr<<"Please check your script "<<script<<endl;
		cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	}
	f>>structure.energyExternal; // Energy in Hartree
	structure.energyExternal /= convToHartree; // convert in internal unit ( Phys or normalized)
	structure.hasExternalE=true;
	if(type>0) 
	{
		double  convForce=convToBohr/convToHartree;
		for (size_t i = 0; i < structure.atoms.size(); ++i)
		{	
			Atom &ai = structure.atoms.at(i);
			for(int j=0;j<3;j++)
				f>>ai.fExternal[j];
			for(int j=0;j<3;j++)
				ai.fExternal[j] *= convForce; // convert in internal unit ( Phys or normalized)
		}
		structure.hasExternalF=true;
	}
	f.close();
}
void External::compute(Structure& structure, int type) const
{
	if (normalize) structure.toPhysicalUnits(meanEnergy, convEnergy, convLength, convCharge);
	std::string fileNameIn  = createInputFile(structure, type);
	std::string fileNameOut = addInt(prefixFileName,  structure.index);
	fileNameOut += ".out";
	std::string command = script + std::string(" ") + fileNameIn + std::string(" ") + fileNameOut;
	system(command.c_str());
	getResults(structure, type);
	if(normalize) structure.toNormalizedUnits(meanEnergy,  convEnergy,  convLength,  convCharge);
}

void External::calculateEnergy(Structure& structure) const
{
	if(script=="unknown")
	{
		structure.energyExternal = 0.0;
		return;
	}
	//cout<<"structure.hasExternalE="<<structure.hasExternalE<<endl;
	if(structure.hasExternalE) return;
	compute(structure,0);
}
void External::calculateForces(Structure& structure) const
{
	if(script=="unknown")
	{
		structure.energyExternal = 0.0;
		int nAtoms = structure.numAtoms;
		for(int i=0;i<3;i++)
		for(int ia=0;ia<nAtoms;ia++)
			structure.atoms[ia].fExternal[i] = 0.0;
		return;
	}
	if(structure.hasExternalF) return;
	compute(structure,1);
}
void External::printForces(Structure& structure, Log log) const
{
    	// Print forces for control
	log << "\n";
	log << "External forces:\n";
	for (size_t i = 0; i < structure.atoms.size(); ++i)
	{
		Atom &ai = structure.atoms.at(i);
		std::string s = structure.elementMap.symbol(ai.element);
        	log << strpr("%10d %2s %16.8E %16.8E %16.8E\n",
                                i+1,
                                s.c_str(),
                                ai.fExternal[0],
                                ai.fExternal[1],
                                ai.fExternal[2]
                                );
	}
	log << "-----------------------------------------"
               "--------------------------------------\n";
	log << "External energy : ";
	log << strpr("%16.8E\n", structure.energyExternal);
	log << "========================================="
               "======================================\n";
}
