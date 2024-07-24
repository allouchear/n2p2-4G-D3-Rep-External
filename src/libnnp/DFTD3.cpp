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
#include "DFTD3.h"
#include "utility.h"
#include <cstdlib>   // atof, atoi
#include <iostream> 

using namespace std;
using namespace nnp;

 void DFTD3::setnormalization(double mEnergy, double cEnergy, double cLength, double cCharge)
{
    	normalize = true;
	meanEnergy = mEnergy;
	convEnergy = cEnergy;
	convLength = cLength;
	convCharge = cCharge;
}
void DFTD3::setcontoatomicunit(double cEnergy, double cLength)
{
    	// Convert From Phys to Atomic unit
	convToHartree = cEnergy;
	convToBohr = cLength;
}
void DFTD3::calculateEnergy(Structure& structure) const
{
	if(name=="no")
	{
		structure.energyDFTD3 = 0.0;
		return;
	}
	//cout<<"structure.hasDFTD3E="<<structure.hasDFTD3E<<endl;
	if(structure.hasDFTD3E) return;
	int v = version;
	char* nam = strdup(name.c_str());
	double edisp = 0.0;
	int nAtoms;
	int* atnum;
	double latVecs[9];
	double* coords;
	double P[5];

	if(normalize) structure.toPhysicalUnits(meanEnergy, convEnergy, convLength, convCharge);
	nAtoms = structure.numAtoms;
	atnum = new int [nAtoms];
        for (int ia = 0; ia < nAtoms; ia++)
	{
		atnum[ia] = int(structure.elementMap.atomicNumber(structure.atoms[ia].element));
	}
//////////////////////////
	/*
        for (int ia = 0; ia < structure.numAtoms; ia++)
        {
		cout<<atnum[ia]<<" ";
		for(int i=0;i<3;i++)
			cout<<structure.atoms[ia].r[i]<<" ";
		cout<<endl;
	}
	*/
//////////////////////////

	coords = new double [3*nAtoms];

	int k=0;
	for(int i=0;i<3;i++)
        for (int ia = 0; ia < nAtoms; ia++)
        {
      		coords [k] = structure.atoms[ia].r[i]*convToBohr;
		k++;
	}
    	if (!structure.isPeriodic)
	{
		if(!strcmp(nam,"custom"))
		{
			for(int i=0;i<5;i++) P[i] = pars[i];
			if(verbose>0) cout<<"call dftd3cepars_ "<<endl;
			dftd3cepars_(P,&v, &nAtoms, atnum, coords, &edisp);
		}
		else
		{
			if(verbose>0) cout<<"call dftd3ce_ "<<endl;
			dftd3ce_(nam,&v, &nAtoms, atnum, coords, &edisp);
		}
	}
	else
	{
		int k = 0;
		for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
		{
			latVecs [k] = structure.box[i][j]*convToBohr;
			k++;
		}
		if(!strcmp(nam,"custom"))
		{
			for(int i=0;i<5;i++) P[i] = pars[i];
			if(verbose>0) cout<<"call dftd3cpbcepars_ "<<endl;
			dftd3cpbcepars_(P, &v, &nAtoms, atnum, coords, latVecs, &edisp);
		}
		else
		{
			if(verbose>0) cout<<"call dftd3cpbce_ "<<endl;
			dftd3cpbce_(nam,&v, &nAtoms, atnum, coords, latVecs, &edisp);
		}
	}
	if(verbose>0) cout<<"edisp = "<<edisp<<endl;
	edisp /=  convToHartree; // dfdt3 energy are given in hartree
	structure.energyDFTD3 = edisp;
	delete [] coords;
	delete [] atnum;
	free(nam);
	structure.hasDFTD3E=true;
	if(normalize) structure.toNormalizedUnits(meanEnergy,  convEnergy,  convLength,  convCharge);
}
void DFTD3::calculateForces(Structure& structure) const
{
	if(name=="no")
	{
		structure.energyDFTD3 = 0.0;
		int nAtoms = structure.numAtoms;
		for(int i=0;i<3;i++)
		for(int ia=0;ia<nAtoms;ia++)
			structure.atoms[ia].fDFTD3[i] = 0.0;
		return;
	}
	if(structure.hasDFTD3F) return;
	int v = version;
	char* nam = strdup(name.c_str());
	double edisp = 0.0;
	int nAtoms;
	int* atnum;
 	double stress[9];
	double latVecs[9];
	double* coords;
	double* grads;
	double P[5];
	if(normalize) structure.toPhysicalUnits(meanEnergy, convEnergy, convLength, convCharge);

	nAtoms = structure.numAtoms;
	atnum = new int [nAtoms];
        for (int ia = 0; ia < nAtoms; ia++)
	{
		atnum[ia] = int(structure.elementMap.atomicNumber(structure.atoms[ia].element));
	}

	coords = new double [3*nAtoms];
	grads = new double [3*nAtoms];

	int k=0;
	for(int i=0;i<3;i++)
        for (int ia = 0; ia < nAtoms; ia++)
        {
      		coords [k] = structure.atoms[ia].r[i]*convToBohr;
		k++;
	}
    	if (!structure.isPeriodic)
	{
		if(!strcmp(nam,"custom"))
		{
			for(int i=0;i<5;i++) P[i] = pars[i];
			if(verbose>0) cout<<"call dftd3cpars_ "<<endl;
      			dftd3cpars_(P,&v, &nAtoms, atnum, coords, &edisp, grads);
		}
		else
		{
			if(verbose>0) cout<<"call dftd3c_ "<<endl;
      			dftd3c_(nam,&v, &nAtoms, atnum, coords, &edisp, grads);
		}
	}
	else
	{
		int k = 0;
		for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
		{
			latVecs [k] = structure.box[i][j]*convToBohr;
			k++;
		}
		if(!strcmp(nam,"custom"))
		{
			for(int i=0;i<5;i++) P[i] = pars[i];
			if(verbose>0) cout<<"call dftd3cpbcpars_ "<<endl;
      			dftd3cpbcpars_(P,&v, &nAtoms, atnum, coords, latVecs, &edisp, grads, stress);
		}
		else
		{
			if(verbose>0) cout<<"call dftd3cpbc_ "<<endl;
      			dftd3cpbc_(nam,&v, &nAtoms, atnum, coords, latVecs, &edisp, grads, stress);
		}
	}
	if(verbose>0) cout<<"edisp = "<<edisp<<endl;
	edisp /= convToHartree;
	structure.energyDFTD3 = edisp;
	structure.hasDFTD3E=true;// if calculateEnergy called , dftd3 not run because we have already calculate energy
	k = 0;
	for(int i=0;i<3;i++)
	for(int ia=0;ia<nAtoms;ia++)
	{
		structure.atoms[ia].fDFTD3[i] = -grads[k]/convToHartree*convToBohr; 
		k++;
	}
	delete [] coords;
	delete [] grads;
	delete [] atnum;
	free(nam);
	structure.hasDFTD3F=true;
	if(normalize) structure.toNormalizedUnits(meanEnergy,  convEnergy,  convLength,  convCharge);
}
void DFTD3::printForces(Structure& structure, Log log) const
{
    	// Print forces for control
	log << "\n";
	log << "DFTD3 forces:\n";
	for (size_t i = 0; i < structure.atoms.size(); ++i)
	{
		Atom &ai = structure.atoms.at(i);
		std::string s = structure.elementMap.symbol(ai.element);
        	log << strpr("%10d %2s %16.8E %16.8E %16.8E\n",
                                i+1,
                                s.c_str(),
                                ai.fDFTD3[0],
                                ai.fDFTD3[1],
                                ai.fDFTD3[2]
                                );
	}
	log << "-----------------------------------------"
               "--------------------------------------\n";
	log << "DFTD3 energy : ";
	log << strpr("%16.8E\n", structure.energyDFTD3);
	log << "========================================="
               "======================================\n";
}
