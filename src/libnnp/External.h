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

#ifndef EXTERNAL
#define EXTERNAL

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "Structure.h"
#include "Log.h"

namespace nnp {

class External {

private:
    std::string prefixFileName;
    std::string script;
    // Convert from Phys to Normalized unit
    bool normalize;
    double meanEnergy;
    double convEnergy;
    double convLength;
    double convCharge;
    // Convert From Phys to Atomic unit
    double convToHartree;
    double convToBohr;

    std::string createInputFile(Structure& structure, int type) const;
    void getResults(Structure& structure, int type) const;
    void compute(Structure& structure, int type) const;

public:
    /// Default constructor.
    External() : prefixFileName           ("external"),
		script                    ("unknown" ),
		meanEnergy                (0.0       ),
		convEnergy                (1.0       ),
		convLength                (1.0       ),
		convCharge                (1.0       ),
		convToHartree             (1.0       ),
		convToBohr                (1.0       )
    { 
    };
    void setnormalization(double mEnergy, double cEnergy, double cLength, double cCharge);
    void setcontoatomicunit(double cEnergy, double cLength);
    void setscript(std::string name){ script = name; };
    void setprefixfilename(std::string name){ prefixFileName = name; };
    void calculateEnergy(Structure& structure) const;
    void calculateForces(Structure& structure) const;
    void printForces(Structure& structure, Log log) const;
};

}
#endif // EXTERNAL
