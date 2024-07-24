#!/usr/bin/env python

# Written by ALLOUCHE, modified version of generate-wcsf-radial.py

import numpy as np
import sys

# Settings
elements = ["H","O","Al", "Si"]
#mode     = "center"
mode     = "shift"
#c=0.529177
c=1.0
r_0      = 1.0/c
#r_c      =  15.1178150222 
r_c      =  10.0
r_N      = r_c - 0.5/c
N        = 8

grid = np.linspace(r_0, r_N, N)
dr = (r_N - r_0) / (N - 1)

sys.stdout.write("# Generating radial symmetry function set:\n")
sys.stdout.write("# mode  = {0:9s}\n".format(mode))
sys.stdout.write("# r_0   = {0:9.3E}\n".format(r_0))
sys.stdout.write("# r_c   = {0:9.3E}\n".format(r_c))
sys.stdout.write("# r_N   = {0:9.3E}\n".format(r_N))
sys.stdout.write("# N     = {0:9d}\n".format(N))
sys.stdout.write("# grid  = " + " ".join(str(r) for r in grid) + "\n")

if mode == "center":
    eta_grid = 1.0 / (2.0 * grid**2)
    rs_grid = [0.0] * N
elif mode == "shift":
    eta_grid = [1.0 / (2.0 * dr * dr)] * N
    rs_grid = grid

for e in elements:
    sys.stdout.write("# Radial symmetry functions for element {0:2s}\n".format(e))
    for (eta, rs) in zip(eta_grid, rs_grid):
        sys.stdout.write("symfunction_short {:2s} 12 {:9.3E} {:9.3E} {:9.3E}\n".format(e, eta, rs, r_c))
    sys.stdout.write("\n")
