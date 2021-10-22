#!/usr/bin/env python
"""
models.py

NIRS optimization models
@author Domenico Salvagnin dominiqs@gmail.com
"""

from __future__ import print_function

import sys
sys.path.append('../python')

from batchUtils import locate
import cplex
import time
from instances import *
from nirs import *
from models import *
from solve import *

instances = [f for f in locate("ROI*.txt", "./data/down2/")]
nSources = [1, 2, 4, 8, 16]
nDetectors = [1, 2, 4, 8, 16]
covWeight = [0.0, 1.0, 2.0, 3.0, 5.0, 10.0]
#minRhoOpt = [10, 15, 20, 25, 30]
algos = ['oldMachado', 'newMachado', 'bigModel']

algoDict = {
	'oldMachado': solveOldMachado,
	'newMachado': solveNewMachado,
	'bigModel': solveBigModel,
	'bigIter': solveBigIterative
}


if __name__ == '__main__':
	outfile = open("results.csv", "w", buffering=1)
	print("instance;nSources;nDetectors;covWeight;method;sources;detectors;time;coverage;signal;primal;dual", file=outfile)
	data = NIRSData()
	data.read(scalpFile, voxelFile, allPMDFfile, wFile)
	for f in instances:
		for nS in nSources:
			for nD in nDetectors:
				if (nS < nD):
					continue
				for cw in covWeight:
					instance = readFromFile(data, nS, nD, f)
					instance.covWeight = cw
					for a in algos:
						func = algoDict[a]
						try:
							start = time.time()
							solution = func(data, instance)
							end = time.time()
							line = "%s;%i;%i;%g;%s;%s;%s;%g;%g;%g;%g;%s" % (f, nS, nD, cw, a,
								str(solution.sources), str(solution.detectors), end-start,
								solution.coverage, solution.signal, solution.objval, str(solution.dualBound))
						except:
							line = "%s;%i;%i;%g;%s;ERR;ERR;ERR;ERR;ERR;ERR;ERR" % (f, nS, nD, cw, a)
						print("#####", line)
						print(line, file=outfile)
