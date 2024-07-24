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

#ifndef REPULSION
#define REPULSION

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "Structure.h"
#include "utility.h"
#include "Log.h"

namespace nnp {

/// Setup data for Repulsive potential
// V = E(r) + S(r) = Z1 Z2 /r phi(r/a) + S(r)
// phi(r) = sum c_i exp(-d_i r/a) ; i=1..4
// a  = .8853 (Z1**alpha + Z2**alpha)**beta
// S(r) = 
//      C if r<rIn
//      A/3(r-rIn)**3 + B/4(r-rIn)**4 + C rIn<r<rOut
//      0 if r>rOut
//     A = ( -3 E'(rOut) +(rOut-rIn)*E"(rOut) ) / (rOut-rIn)**2 
//     B = (  2 E'(rOut) -(rOut-rIn)*E"(rOut) ) / (rOut-rIn)**3 
//     C =  -E(rOut) + 1/2.0 (rOut-rIn)* E'(rOut) -1/12.0 (rOut-rIn)**2*E"(rOut) 

class Repulsion {
    /** Enumerates different potential types 
     */
    enum RepType
    {
        /** Rep type not assigned yet.
         */
        RT_UNKNOWN,
        /**ZBL Potential
         */
        RT_ZBL,
        /**Moliere  Potential
         */
        RT_MOLIERE,
        /**KrC  Potential
         */
        RT_KRC,
    };

private:
    double rIn, rOut;
    int nc = 4;
    double c[4];
    double d[4];
    double alpha, beta;// a  = .8853 (Z1**alpha + Z2**alpha)**beta
    RepType type;
    std::string name;
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
	normalize = false;
	meanEnergy                =0.0;
	convEnergy                =1.0;
	convLength                =1.0;
	convCharge                =1.0;
	convToHartree             =1.0;
	convToBohr                =1.0;
    }
    void convrinrout()
    {
	rIn *= convToBohr;
	rOut *= convToBohr;
    }

    void setZBL()
    {
	name = "ZBL";
	type = RT_ZBL;
	nc = 4;
    	alpha = 0.23;
	beta  = -1.0;
	c[0]  = 0.02817;
	c[1]  = 0.28022;
	c[2]  = 0.50986;
	c[3] = 1.0;
	for(int i=0;i<3;i++) c[3] -= c[i];
	d[0] = 0.20162;
	d[1] = 0.40290;
	d[2] = 0.94229;
	d[3] = 3.1998;
    }
    void setKrC()
    {
	name = "KrC";
	type = RT_KRC;
	nc = 3;
    	alpha = 0.5;
	beta  = -2.0/3.0;
	c[0]  = 0.190945;
	c[1]  = 0.473674;
	c[2]  = 0.335381;
	c[3] = 0.0;
	d[0] = 0.278544;
	d[1] = 0.637174;
	d[2] = 1.919249;
	d[3] = 0.0;
    }
    void setMoliere()
    {
	type = RT_MOLIERE;
	name = "Moliere";
	nc = 3;
    	alpha = 0.5;
	beta  = -2.0/3.0;
	c[0]  = 0.35;
	c[1]  = 0.55;
	c[2]  = 0.10;
	c[3] = 0.0;

	d[0] = 0.3;
	d[1] = 1.2;
	d[2] = 6.0;
	d[3] = 0.0;
    }
    void setUnknown()
    {
	type = RT_UNKNOWN;
	nc=0;
	for(int i=0;i<4;i++) c[0] = 0;
	for(int i=0;i<4;i++) d[0] = 0;
    }
    double a(double Zi, double Zj) const
    {
	return 0.8853*pow(pow(Zi,alpha) + pow(Zj,alpha), beta);
    }
    double phi(double x) const
    {
	double p = 0;
	for (int i=0;i<nc;i++)
		p += c[i]*exp(-d[i]*x);
	return p;
    }
    double dphi(double x) const
    {
	double dp = 0;
	for (int i=0;i<nc;i++)
		dp += -d[i]*c[i]*exp(-d[i]*x);
	return dp;
    }
    double d2phi(double x) const
    {
	double d2p = 0;
	for (int i=0;i<nc;i++)
		d2p += d[i]*d[i]*c[i]*exp(-d[i]*x);
	return d2p;
    }
    double E(double Zi, double Zj, double rij) const
    {
	double aij = a(Zi,Zj);
	double x = rij/aij;
	double pij = Zi*Zj/aij; // E = pij*phi(x) / x
	return pij*phi(x)/x;
    }
    double dE(double Zi, double Zj, double rij) const // dE/dr
    {
	double aij = a(Zi,Zj);
	double x = rij/aij;
	double pij = Zi*Zj/aij; // E = pij * phi(x) / x
	
	return pij/aij*(dphi(x)-phi(x)/x)/x;
    }
    double d2E(double Zi, double Zj, double rij) const // d2E/dr^2
    {
	double aij = a(Zi,Zj);
	double x = rij/aij;
	double pij = Zi*Zj/aij; // E = pij *1/x * phi(x)
	
	return pij/aij/aij* (d2phi(x)-2*dphi(x)/x + 2*phi(x)/x/x)/x;
    }
    double S(double Zi, double Zj, double rij) const
    {
	if( rij>=rOut) return 0.0;
	double E0 = E(Zi,Zj,rOut);
	double Ep = dE(Zi,Zj,rOut);
	double Es = d2E(Zi,Zj,rOut);
	double r12 = (rOut-rIn);

	double C = -E0 + 0.5*r12*Ep - 1/12.0*r12*r12*Es;

	if( rij<rIn) return C;

	double rr = 1.0/r12/r12;
	double A = (-3*Ep+r12*Es)*rr;
	double B = (2*Ep-r12*Es)*rr/r12;
	double r = (rij-rIn);
	double r3 = r*r*r;
	double r4 = r3*r;
	double s = A/3*r3+B/4*r4+C;

	return s;
    }
    double dS(double Zi, double Zj, double rij) const // dS/dr
    {


	if( rij<rIn) return 0.0;
	if( rij>=rOut) return 0.0;

	double Ep = dE(Zi,Zj,rOut);
	double Es = d2E(Zi,Zj,rOut);
	double r12 = (rOut-rIn);

	double rr = 1.0/r12/r12;
	double A = (-3*Ep+r12*Es)*rr;
	double B = (2*Ep-r12*Es)*rr/r12;
	double r = (rij-rIn);
	double r2 = r*r;
	double r3 = r2*r;
	double ds = A*r2+B*r3;

	return ds;
    }
    double energy(double Zi, double Zj, double rij) const
    {
	return E(Zi,Zj,rij)+S(Zi,Zj,rij);
    }
    double dEnergy(double Zi, double Zj, double rij)  const // d/dr
    {
	return dE(Zi,Zj,rij)+dS(Zi,Zj,rij);
    }
/*
    double energy(Structure& structure, Atom& ai, Atom& aj) const
    {
	double Zi = structure.elementMap.atomicNumber(ai.element);
	double Zj = structure.elementMap.atomicNumber(ai.element);
        double const rij = (ai.r - aj.r).norm();
	return energy(Zi,Zj,rij);
    }
    Vec3D getForce(Structure& structure, const Atom& ai, const Atom& aj) const
    {
	double Zi = structure.elementMap.atomicNumber(ai.element);
	double Zj = structure.elementMap.atomicNumber(ai.element);
	Vec3D v = ai.r - aj.r;
        double const rij = v.norm();
	v *= -dEnergy(Zi,Zj,rij)/rij;
	return v;
    }
*/
    double energy(Structure& structure, Atom& ai) const // structure in phys unit
    {
	double Zi = structure.elementMap.atomicNumber(ai.element);
	size_t numNeighbors = ai.numNeighbors;

	double e = 0;
	// convert distance in Bohr
	for (size_t j = 0; j < numNeighbors; j++)
	{
        	Atom::Neighbor& nj = ai.neighbors[j];
		if( ai.index==nj.index) continue;
        	double rij = nj.d;
		rij *= convToBohr;
		
		if (rij < rOut) // rOut is already in Bohr
		{
			double Zj = structure.elementMap.atomicNumber(nj.element);
			e += energy(Zi,Zj,rij);
		}
	}
	// convert energy from Hartree to phys unit
	e /= convToHartree;
	return e;
    }

    void computeForce(Structure& structure, Atom& ai) const  // structure in phys unit
    {
	double convForce = 1.0/convToHartree*convToBohr; // convert forces to phys unit
	double Zi = structure.elementMap.atomicNumber(ai.element);
	size_t numNeighbors = ai.numNeighbors;
        // Reset forces.
        ai.fRep = Vec3D{};
	for (size_t j = 0; j < numNeighbors; j++)
	{
        	Atom::Neighbor& nj = ai.neighbors[j];
		if( ai.index==nj.index) continue;
        	double rij = nj.d;
		rij *= convToBohr;
		if (rij < rOut) // rOut is already in Bohr
		{
			double Zj = structure.elementMap.atomicNumber(nj.element);
			Vec3D v = nj.dr; // ai.r - aj.r;
			v *= convToBohr;
			v *= -dEnergy(Zi,Zj,rij)/rij;
			ai.fRep += v;
		}
	}
	ai.fRep *= convForce; // convert forces to phys unit 
    }

public:
    /// Default constructor.
    Repulsion(): type(RT_UNKNOWN){ 
    		setDefaultConv();
		setUnknown(); 
		rIn=0; 
		rOut=10.0; 
		name="Unknown";
		};
    Repulsion(std::string rt ) {
    	setDefaultConv();
	transform(rt.begin(), rt.end(), rt.begin(), ::tolower);
	if(rt=="zbl") setZBL();
	else if(rt=="moliere") setMoliere();
	else if(rt=="krc") setKrC();
	else if(rt=="no")
	{
			//std::cout<<" No repulsive potential"<<std::endl;
	}
	else
	{
		std::cout<<" Unknown potential"<<std::endl;
		std::cout<<" The implemented potential are : ZBL, Moliere and KrC"<<std::endl;
		exit(EXIT_FAILURE);
	}
	rIn=0;
	rOut=10.0;
	convrinrout();
	}
    void setrIn(double r){ rIn=r; rIn  *= convToBohr; };// We assume that r is in Physical unit, not yet normalized
    void setrOut(double r){ rOut=r; rOut *= convToHartree; };

    void setnormalization(double mEnergy, double cEnergy, double cLength, double cCharge);
    void setcontoatomicunit(double cEnergy, double cLength);
    void calculateEnergy(Structure& structure) const;
    void calculateForces(Structure& structure) const;
    void printAvailables(Log log) const {
		log<<" Available potential: ZBL, Moliere, KrC\n";
	};
    void check(Log log) const {
		if(rIn>rOut)
		{	
			log<<" rIn must be > rOut for repulsive potential\n";
			exit(EXIT_FAILURE);
		}
	};

    void printParameters(Log log) const{
	if(type == RT_UNKNOWN) return;
	log << strpr("\t type  = %s\n",name);
	log << strpr("\t alpha = %lf\n",alpha);
	log << strpr("\t beta  = %lf\n",beta);
	log << strpr("\t c     ="); for(int i=0;i<nc;i++) log << strpr(" %lf",c[i]); log<<"\n";
	log << strpr("\t d     ="); for(int i=0;i<nc;i++) log << strpr(" %lf",d[i]); log<<"\n";
	log << "\tV = E(r) + S(r) = Z1 Z2 /r phi(r/a) + S(r)\n";
	log << "\tphi(r) = sum c_i exp(-d_i r/a) ; i=1..4\n";
	log << "\t\ta  = .8853 (Z1**alpha + Z2**alpha)**beta\n";
	log<<"\tS(r) = \n";
	log<<"\t\tC if r<rIn\n";
	log<< "\t\tA/3(r-rIn)**3 + B/4(r-rIn)**4 + C rIn<r<rOut\n";
	log << "\t\t0 if r>rOut\n";
	log<< "\t\tA = ( -3 E'(rOut) +(rOut-rIn)*E\"(rOut) ) / (rOut-rIn)**2\n";
	log<< "\t\tB = (  2 E'(rOut) -(rOut-rIn)*E\"(rOut) ) / (rOut-rIn)**3\n";
	log<<"\t\t C =  -E(rOut) + 1/2.0 (rOut-rIn)* E'(rOut) -1/12.0 (rOut-rIn)**2*E\"(rOut)\n";
    }
    void printForces(Structure& structure, Log log) const;

};

}
#endif // REPULSION
