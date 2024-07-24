#!/usr/bin/env python

# Written by ALLOUCHE, modified version of generate-wcsf-angular*.py
import numpy as np
import sys

# Settings
elements = ["H","O","Al","Si"]
mode     = "center"
#c=0.529177
c=1.0
r_0      = 1.0/c
#r_c      =  15.1178150222 
r_c      =  10.0
r_N      = r_c - 0.5/c

N        = 3
zetas    = [1.0]

grid = np.linspace(r_0, r_N, N)
dr = (r_N - r_0) / (N - 1)

sys.stdout.write("# Generating narrow angular symmetry function set:\n")
sys.stdout.write("# mode  = {0:9s}\n".format(mode))
sys.stdout.write("# r_0   = {0:9.3E}\n".format(r_0))
sys.stdout.write("# r_c   = {0:9.3E}\n".format(r_c))
sys.stdout.write("# r_N   = {0:9.3E}\n".format(r_N))
sys.stdout.write("# N     = {0:9d}\n".format(N))
sys.stdout.write("# grid  = " + " ".join(str(r) for r in grid) + "\n")
sys.stdout.write("# zetas = " + " ".join(str(z) for z in zetas) + "\n")

if mode == "center":
    eta_grid = 1.0 / (2.0 * grid**2)
    rs_grid = [0.0] * N
elif mode == "shift":
    eta_grid = [1.0 / (2.0 * dr * dr)] * N
    rs_grid = grid

for e in elements:
    sys.stdout.write("# Weighted Narrow angular symmetry functions for element {0:2s}\n".format(e))
    for e1 in elements:
        elements_reduced = elements[elements.index(e1):]
        for (eta, rs) in zip(eta_grid, rs_grid):
            for zeta in zetas:
                for lambd in [-1.0, 1.0]:
		     #symfunction_short <element-central> 13 <eta> <rshift> <lambda> <zeta> <rcutoff>
                     sys.stdout.write("symfunction_short {:2s} 13 {:14.8f} {:3.1f} {:>4.1f} {:>4.1f} {:14.8f}\n".format(e, eta, 0.0, lambd, zeta, r_c))
        sys.stdout.write("\n")
