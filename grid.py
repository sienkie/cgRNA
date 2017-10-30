#!/usr/bin/env python

import sys
from CABS.cabsDock.cabs import CabsLattice

cl = CabsLattice(grid_spacing=float(sys.argv[1]), r12=(3.37, 4.34), r13=(4.97, 7.29))
#for v in cl.vectors:
#	print v
print(len(cl.vectors))
