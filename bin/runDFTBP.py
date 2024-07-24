import sys

import numpy as np
import os

class Geometry():
	def __str__(self):
		return "Geometry"

	def __init__(self,
		):
		self._ctype=0
		self._charge=0
		self._multiplicity=1
		self._symbols = []
		self._positions = []
		self._lattice = []
		self._forces = []
		self._typeatoms = [] # num of atom type : integer : 1, 2 ...
		self._types = [] # list of atom type (string)
		self._energy=None

	def readDetailed(self,fname):
		check_file = os.path.isfile(fname)
		if not check_file:
			return False, "I cannot open "+fname+" file"
		f = open(fname, "r")
		lines = f.readlines()
		rforces=False
		self._forces = []
		for line in lines:
			# Total energy:                      -92.5951510969 H        -2519.6423 eV
			if line.find("Total energy:") != -1: 
				self._energy = np.double(line.split()[2])
			if rforces: 
				if len(line.split())==0: 
					rforces = False 
				elif len(line.split())!=4: 
					print("Error in reading of forces in ",fname," file", file=sys.stderr)
					sys.exit(1)
				else: 
					x = np.double(line.split()[1])
					y = np.double(line.split()[2])
					z = np.double(line.split()[3])
					self._forces += [[x,y,z]]
			if line.find("Total Forces") != -1 and len(line.split())==2: 
				rforces=True
		

	def write(self,f):
		spr =""
		for   s in  self._types:
			spr += " " + s
		print("Geometry = GenFormat {",file=f)
		if len(self._lattice)>0:
			print(len( self._symbols)," S",file=f)
		else:
			print(len( self._symbols)," C",file=f)
		print(spr,file=f)
		print("#  Index Type  Coordinates ",file=f)
		for k, (s, t, pos) in enumerate(zip(self._symbols, self._typeatoms, self._positions)):
			print('{:4d} {:4d} {:20.10e} {:20.10e} {:20.10e}'.format(k+1, t, pos[0], pos[1], pos[2]),file=f)
		if len(self._lattice)>0:
			print('{:20.10e} {:20.10e} {:20.10e}'.format(0.0,0.0,0.0),file=f)
			for v in self._lattice:
				print('{:20.10e} {:20.10e} {:20.10e}'.format(v[0], v[1], v[2]),file=f)
			
		print("}",file=f)
		if self.ctype == 1:
			print("Analysis {",file=f)
			print(" CalculateForces = Yes",file=f)
			print("}",file=f)
	def writeEnergyAndForces(self,fname):
		f = open(fname, "w")
		print('{:20.10e}'.format(self._energy),file=f)
		for k, fo in enumerate(self._forces):
			print('{:20.10e} {:20.10e} {:20.10e}'.format(fo[0], fo[1], fo[2]),file=f)
		f.close()

	def settypes(self):
		self._types =  set(self._symbols)
		self._typeatoms = [1]*len(self._symbols)
		i=0
		for t in self._types:
			i += 1
			for k, s in enumerate(self._symbols):
				if t==s:
					self._typeatoms[k] = i
					
	def convertInAng(self):
		ca = 0.529177
		self._lattice =  np.array(self._lattice)*ca
		self._positions =  np.array(self._positions)*ca

	def read(self, fname):
		check_file = os.path.isfile(fname)
		if not check_file:
			return False, "I cannot open "+fname+" file"
		f = open(fname, "r")
		lines = f.readlines()
		if len(lines[0].split()) != 3:
			return False, "First line in "+fname+" must contain 3 values : type, charge , multiplicity"
		self._ctype = int(lines[0].split()[0])
		self._charge = np.double(lines[0].split()[1])
		self._multiplicity = int(lines[0].split()[2])
		for i in range(1,len(lines)):
			if len(lines[i].split()) != 4:
				return False, "Each line in "+fname+" must contain symbol, x, y, z in atomic unit"
			s = lines[i].split()[0]
			x = np.double(lines[i].split()[1])
			y = np.double(lines[i].split()[2])
			z = np.double(lines[i].split()[3])
			if s.lower() =='tv':
				self._lattice += [[x,y,z]]
			else:
				self._symbols += [s]
				self._positions += [[x,y,z]]
				self._forces += [[0.0,0.0,0.0]]
				
		self.settypes()
		self.convertInAng()
				
	@property
	def ctype(self):
		return self._ctype

	@property
	def charge(self):
		return self._charge

	@property
	def multiplicity(self):
		return self._multiplicity

	@property
	def symbol(self):
		return self._symbols

	@property
	def positions(self):
		return self._positions

	@property
	def lattice(self):
		return self._lattice

	@property
	def forces(self):
		return self._forces

	@property
	def typeatoms(self):
		return self._typeatoms

	@property
	def types(self):
		return self._types

def setFiles(fname):
	base=os.path.basename(fname)
	base=os.path.splitext(base)[0]
	directory=os.path.dirname(fname)+"dftbp_"+base
	#print("creating directory...")
	if not os.path.exists(directory):
		os.makedirs(directory)

	hsd = os.path.join(directory, 'dftb_in.hsd')
	#output = os.path.join(directory, 'output.txt')
	output = os.path.join('output.txt')
	detailed = os.path.join(directory, 'detailed.out')

	return directory, hsd, output, detailed

'''
def addHamiltonian(f, method):
	print("Hamiltonian = xTB {", file=f)
	print("  Method=\""+method+"\"",file=f)
	print("  kPointsAndWeights = SuperCellFolding {",file=f)
	print("       2   0   0",file=f)
	print("       0   2   0",file=f)
	print("       0   0   2",file=f)
	print("       0.5 0.5 0.5",file=f)
	print("  }", file=f)
	print("}", file=f)
'''

def addHamiltonian(f, hamiltonianfname):
	check_file = os.path.isfile(hamiltonianfname)
	if not check_file:
		print("Error : I cannnot read Hamiltonian from ",fname," file", file=sys.stderr)
		sys.exit(1)
	fh= open( hamiltonianfname, "r")
	lines = fh.readlines()
	for line in lines:
		f.write(line)
	fh.close()

if len(sys.argv) <4:
	print("ERROR : runDFTBP You mut give the name of inputfil, name of outputfile and name of file containing the Hamiltonian ",file=sys.stderr)
	sys.exit(1)

inputfname=sys.argv[1]
outputfname=sys.argv[2]
hamiltonianfname=sys.argv[3]

directory, hsd, output , detailed = setFiles(inputfname)
#print("directory = ", directory)
#print("hsd = ", hsd)
#print("output = ", output)
#print("detailed = ", detailed)

geom = Geometry()
geom.read(inputfname)

fhsd= open(hsd,'w')
geom.write(fhsd)
addHamiltonian(fhsd, hamiltonianfname)
fhsd.close()

command ="cd " + directory
command +=" ; dftb+ >" + output
#print("command=", command)
os.system(command)

geom.readDetailed(detailed)
geom.writeEnergyAndForces(outputfname)
#print("see file = ", outputfname)

# remover directory
#import shutil
#shutil.rmtree(directory)


