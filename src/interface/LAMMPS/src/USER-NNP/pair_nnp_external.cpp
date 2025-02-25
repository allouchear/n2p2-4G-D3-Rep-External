// Copyright 2018 Andreas Singraber (University of Vienna)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <mpi.h>
#include <cstdlib>
#include <fstream>
#include <string.h>
#include <vector>
#include "pair_nnp_external.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "utils.h"
#include "Atom.h" // nnp::Atom
#include "utility.h" // nnp::

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */

PairNNPExternal::PairNNPExternal(LAMMPS *lmp) : Pair(lmp)
{
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairNNPExternal::~PairNNPExternal()
{
}

/* ---------------------------------------------------------------------- */

void PairNNPExternal::compute(int eflag, int vflag)
{
  if(eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  // Clear structure completely.
  structure.reset();

  // Set simulation box.
  if (domain->nonperiodic < 2) structure.isPeriodic = true;
  else structure.isPeriodic = false;
  if (structure.isPeriodic) structure.isTriclinic = (bool)domain->triclinic;
  structure.box[0][0] = domain->h[0] * cflength;
  structure.box[0][1] = 0.0;
  structure.box[0][2] = 0.0;
  structure.box[1][0] = domain->h[5] * cflength;
  structure.box[1][1] = domain->h[1] * cflength;
  structure.box[1][2] = 0.0;
  structure.box[2][0] = domain->h[4] * cflength;
  structure.box[2][1] = domain->h[3] * cflength;
  structure.box[2][2] = domain->h[2] * cflength;

  // Fill structure with atoms.
  for (int i = 0; i < atom->nlocal; ++i)
  {
    structure.atoms.push_back(nnp::Atom());
    nnp::Atom& a = structure.atoms.back();
    a.r[0] = atom->x[i][0] * cflength;
    a.r[1] = atom->x[i][1] * cflength;
    a.r[2] = atom->x[i][2] * cflength;
    a.element = atom->type[i] - 1;
    structure.numAtoms++;
  }

  // Write "input.data" file to disk.
  structure.writeToFile(string(directory) + "input.data");

  // Run external command and throw away stdout output.
  string com = "cd " + string(directory) + "; " + string(command) + " > external.out";
  system(com.c_str());

  // Read back in total potential energy.
  ifstream f;
  f.open(string(directory) + "energy.out");
  if (!f.is_open()) error->all(FLERR,"Could not open energy output file");
  string line;
  double energy = 0.0;
  getline(f, line); // Ignore first line, RuNNer and n2p2 have header here.
  while (getline(f, line))
  {
    if ((line.size() > 0) && (line.at(0) != '#')) // Ignore n2p2 header.
    {
      // 4th columns contains NNP energy.
      sscanf(line.c_str(), "%*s %*s %*s %lf", &energy);
    }
  }
  f.close();
  energy /= cfenergy;

  // Add energy contribution to total energy.
  if (eflag_global)
     ev_tally(0,0,atom->nlocal,1,energy,0.0,0.0,0.0,0.0,0.0);

  // Read back in forces.
  f.open(string(directory) + "nnforces.out");
  if (!f.is_open()) error->all(FLERR,"Could not open force output file");
  int c = 0;
  double const cfforce = cfenergy / cflength;
  getline(f, line); // Ignore first line, RuNNer and n2p2 have header here.
  while (getline(f, line))
  {
    if ((line.size() > 0) && (line.at(0) != '#')) // Ignore n2p2 header.
    {
      if (c > atom->nlocal - 1) error->all(FLERR,"Too many atoms in force file.");
      double fx;
      double fy;
      double fz;
      sscanf(line.c_str(), "%*s %*s %*s %*s %*s %lf %lf %lf", &fx, &fy, &fz);
      atom->f[c][0] = fx / cfforce;
      atom->f[c][1] = fy / cfforce;
      atom->f[c][2] = fz / cfforce;
      c++;
    }
  }

  f.close();

  // If virial needed calculate via F dot r.
  // TODO: Pressure calculation is probably wrong anyway, tell user only to use
  // in NVE, NVT ensemble.
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairNNPExternal::settings(int narg, char **arg)
{
  int iarg = 0;

  // elements list is mandatory
  if (narg < 1) error->all(FLERR,"Illegal pair_style command : narg<1");

  // element list, mandatory
  int len = strlen(arg[iarg]) + 1;
  elements = new char[len];
  sprintf(elements, "%s", arg[iarg]);
  iarg++;
  em.registerElements(elements);
  vector<FILE*> fp;
  if (screen) fp.push_back(screen);
  if (logfile) fp.push_back(logfile);
  structure.setElementMap(em);

  // default settings
  len = strlen("nnp/") + 1;
  directory = new char[len];
  strcpy(directory,"nnp/");
  len = strlen("nnp-predict 0") + 1;
  command = new char[len];
  strcpy(command,"nnp-predict 0");
  cflength = 1.0;
  cfenergy = 1.0;

  while(iarg < narg) {
    // set NNP directory
    if (strcmp(arg[iarg],"dir") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair_style command : dir");
      delete[] directory;
      len = strlen(arg[iarg+1]) + 2;
      directory = new char[len];
      sprintf(directory, "%s/", arg[iarg+1]);
      iarg += 2;
    // set external prediction command
    } else if (strcmp(arg[iarg],"command") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair_style command: command");
      delete[] command;
      len = strlen(arg[iarg+1]) + 1;
      command = new char[len];
      sprintf(command, "%s", arg[iarg+1]);
      iarg += 2;
    // length unit conversion factor
    } else if (strcmp(arg[iarg],"cflength") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair_style command : cflength");
      cflength = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    // energy unit conversion factor
    } else if (strcmp(arg[iarg],"cfenergy") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair_style command : cfenergy");
      cfenergy = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else 
	{
	string s = "Illegal pair_style command : other";
	s += " arg =";
	s +=string(arg[iarg]);
	//error->all(FLERR,"Illegal pair_style command : other");
	error->all(FLERR,s) ;
	}
  }

  for (auto f : fp)
  {
      fprintf(f, "*****************************************"
                 "**************************************\n");
      fprintf(f, "pair_style nnp/external settings:\n");
      fprintf(f, "---------------------------------\n");
      fprintf(f, "elements = %s\n", elements);
      fprintf(f, "dir      = %s\n", directory);
      fprintf(f, "command  = %s\n", command);
      fprintf(f, "cflength = %16.8E\n", cflength);
      fprintf(f, "cfenergy = %16.8E\n", cfenergy);
      fprintf(f, "*****************************************"
                 "**************************************\n");
      fprintf(f, "CAUTION: Explicit element mapping is not available for nnp/external,\n");
      fprintf(f, "         please carefully check whether this map between LAMMPS\n");
      fprintf(f, "         atom types and element strings is correct:\n");
      fprintf(f, "---------------------------\n");
      fprintf(f, "LAMMPS type  |  NNP element\n");
      fprintf(f, "---------------------------\n");
      int lammpsNtypes = em.size();
      for (int i = 0; i < lammpsNtypes; ++i)
      {
          fprintf(f, "%11d <-> %2s (%3zu)\n",
                  i, em[i].c_str(), em.atomicNumber(i));
      }
      fprintf(f, "*****************************************"
                 "**************************************\n");
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairNNPExternal::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3) error->all(FLERR,"Incorrect args for pair coefficients");

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  maxCutoffRadius = utils::numeric(FLERR,arg[2],false,lmp);

  // TODO: Check how this flag is set.
  int count = 0;
  for(int i=ilo; i<=ihi; i++) {
    for(int j=MAX(jlo,i); j<=jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairNNPExternal::init_style()
{
  if (comm->nprocs > 1)
  {
    error->all(FLERR,"MPI is not supported for this pair_style");
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairNNPExternal::init_one(int i, int j)
{
  // TODO: Check how this actually works for different cutoffs.
  return maxCutoffRadius;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairNNPExternal::write_restart(FILE *fp)
{
    return;
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairNNPExternal::read_restart(FILE *fp)
{
    return;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairNNPExternal::write_restart_settings(FILE *fp)
{
    return;
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairNNPExternal::read_restart_settings(FILE *fp)
{
    return;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairNNPExternal::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}
