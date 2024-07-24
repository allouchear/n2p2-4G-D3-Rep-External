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
#ifndef DFTD3_H
#define DFTD3_H

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "Structure.h"
#include "Log.h"
#include "utility.h"

extern"C" {
void  dftd3c_(char* func, int* version, int* nAtoms, int* atnum, double* coords, double* edisp, double* grads);
void  dftd3cpbc_(char* func, int* version, int* nAtoms, int* atnum, double* coords, double* latVecs, double* edisp, double* grads, double* stress);
void  dftd3ce_(char* func, int* version, int* nAtoms, int* atnum, double* coords, double* edisp);
void  dftd3cpbce_(char* func, int* version, int* nAtoms, int* atnum, double* coords, double* latVecs, double* edisp);

// pars : # Parameters are as follows: (s6, rs6, s18, rs18, alp6)
//        # otherwise known as         (s6, a1,  s8,  a2,   alp6)
void  dftd3cpars_(double* pars, int* version, int* nAtoms, int* atnum, double* coords, double* edisp, double* grads);
void  dftd3cpbcpars_(double* pars, int* version, int* nAtoms, int* atnum, double* coords, double* latVecs, double* edisp, double* grads, double* stress);
void  dftd3cepars_( double* pars, int* version, int* nAtoms, int* atnum, double* coords, double* edisp);
void  dftd3cpbcepars_(double* pars, int* version, int* nAtoms, int* atnum, double* coords, double* latVecs, double* edisp);

}

namespace nnp {

/// Setup data for DFTD3
class DFTD3 {
private:
    int verbose;
    std::string name;
    int version;// 2=> DFT-D2 , 3=> FT-D3 zero-damp, 4=>  DFT-D3(BJ), 5=> DFT-D3 M(zero), 6=>  DFT-D3M(BJ)
    double pars[5]; //double s6,rs6,s18,rs18,alp; // used if name=custom : s6,rs6,s18,rs18,alp
    // Values for version = 2 : s6=1.1 alp =20.0 rs6=0.75      s18=0.0       rs18=0.0  
    // Values for version = 3 : s6=1.0 alp =14.0 rs6=1.217     s18=0.722     rs18=1.0 
    // Values for version = 4 : s6=1.0 alp =14.0 rs6 =0.4289   s18 =0.7875   rs18=4.4407
    // Values for version = 5 : s6=1.0 alp =14.0 rs6 =2.340218 s18 =0.000000 rs18=0.129434
    // Values for version = 6 : s6=1.0 alp =14.0 rs6 =0.012092 s18 =0.358940 rs18=5.938951 
    // Convert from Phys to Normalized unit
    bool normalize;
    double meanEnergy;
    double convEnergy;
    double convLength;
    double convCharge;
    // Convert From Phys to Atomic unit
    double convToHartree;
    double convToBohr;

    void setDefaultConv()
    {
	normalize                 = false;
	meanEnergy                = 0.0;
	convEnergy                = 1.0;
	convLength                = 1.0;
	convCharge                = 1.0;
	convToHartree             = 1.0;
	convToBohr                = 1.0;
    }



    void setdefpars() { 
	if(version==2)
	{
		sets6(1.1);
		setrs6(0.75);
		sets18(0.0);
		setrs18(0.0);
		setalp6(20);
	}
	else if(version==3)
	{
		sets6(1.0);
		setrs6(1.217);
		sets18(0.722);
		setrs18(1.0);
		setalp6(14.0);
	}
	else if(version==4)
	{
		sets6(1.0);
		setrs6(0.4289);
		sets18(0.7875);
		setrs18(4.4407);
		setalp6(14.0);
	}
	else if(version==5)
	{
		sets6(1.0);
		setrs6(2.340218);
		sets18(0.0);
		setrs18(0.12943);
		setalp6(14.0);
	}
	else if(version==6)
	{
		sets6(1.0);
		setrs6(0.012092);
		sets18(0.358940);
		setrs18(5.938951);
		setalp6(14.0);
	}
    };


public:
    /// Default constructor.
    DFTD3(): verbose(0), name("no"), version(4) { setdefpars(); setDefaultConv();}
    DFTD3(std::string nm,int v,int vrb):verbose(vrb), name(nm), version(v){
	transform(name.begin(), name.end(), name.begin(), ::tolower);
	setdefpars();
	setDefaultConv();
	};

    int getVerbose() const { return verbose; };

    //bool isNo() const { return name=="no"; };
    void sets6(double s6) { pars[0] = s6; };
    void setrs6(double rs6) { pars[1] = rs6; };
    void sets18(double s18) { pars[2] = s18; };
    void setrs18(double rs18) { pars[3] = rs18; };
    void setalp6(double alp6) { pars[4] = alp6; };


    void setnormalization(double mEnergy, double cEnergy, double cLength, double cCharge);
    void setcontoatomicunit(double cEnergy, double cLength);
    void calculateEnergy(Structure& structure) const;
    void calculateForces(Structure& structure) const;
    void printAvailableVersions(Log log) const {
	log<<" Available versions:\n 2=> DFT-D2 , 3=> FT-D3 zero-damp, 4=>  DFT-D3(BJ), 5=> DFT-D3 M(zero), 6=>  DFT-D3M(BJ)\n";
	};
    void checkVersion(Log log) const {
	if (version<2 || version>6)
        {
		log<<" Unkown DFTD3 version\n";
    		printAvailableVersions(log);
		log<<" Calaculation stopped\n";
		exit(EXIT_FAILURE);
        }
    }
    void printParameterNames(Log log){
	log << "Custom parameters : \n";
	log << strpr("\tdftd3_s6   = %f\n", pars[0]);
	log << strpr("\tdftd3_rs6  = %f\n", pars[1]);
	log << strpr("\tdftd3_s18  = %f\n", pars[2]);
	log << strpr("\tdftd3_rs18 = %f\n", pars[3]);
	log << strpr("\tdftd3_alp6 = %f\n", pars[4]);
	}
    void printForces(Structure& structure, Log log) const;

};

}
#endif //DFTD3_H
