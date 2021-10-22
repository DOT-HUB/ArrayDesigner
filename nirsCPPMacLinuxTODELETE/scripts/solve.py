#!/usr/bin/env python
"""
models.py

NIRS optimization models
@author Domenico Salvagnin dominiqs@gmail.com
"""

from __future__ import print_function

import sys
sys.path.append('../python')

import cplex
import math
import h5py
import numpy as np
import random
from nirs import *
from models import *


def solve(lp, timelimit, polishafter):
	lp.parameters.timelimit.set(timelimit)
	lp.parameters.mip.tolerances.mipgap.set(0.01)
	lp.parameters.mip.polishafter.time.set(polishafter)
	lp.solve()


def getSolution(cmodel, lp):
	mipmodel = cmodel.mipmodel
	sources = []
	detectors = []
	try:
		values = lp.solution.get_values()
		objval = lp.solution.get_objective_value()
	except:
		return NIRSSolution()
	x = mipmodel.getVars('x')
	y = mipmodel.getVars('y')
	for p in x:
		if values[x[p].modelindex] > 0.5:
			sources.append(p)
	for p in y:
		if values[y[p].modelindex] > 0.5:
			detectors.append(p)
	return NIRSSolution(sources, detectors, objval)


def solveModel(cmodel, lp, timelimit, polishafter):
	solve(lp, timelimit, polishafter)
	solution = getSolution(cmodel, lp)
	cmodel.instance.evalCoverageAndSignal(cmodel.data, solution)
	solution.dualBound = lp.solution.MIP.get_best_objective()
	return solution


def solveOldMachado(data, instance):
	cmodel = NIRSModel(data, instance)
	mipmodel = cmodel.generateOldCompactModel()
	cmodel.mipmodel = mipmodel
	lp = mipmodel.exportCPLEX()
	# lp.write("oldmachado.lp")
	return solveModel(cmodel, lp, 1800, 1500)


def solveNewMachado(data, instance):
	cmodel = NIRSModel(data, instance)
	mipmodel = cmodel.generateCompactModel()
	cmodel.mipmodel = mipmodel
	lp = mipmodel.exportCPLEX()
	# lp.write("newmachado.lp")
	return solveModel(cmodel, lp, 1800, 1500)


def solveBigModel(data, instance):
	cmodel = NIRSModel(data, instance)
	# first solve Machado model in order to be able to normalize signal
	print("First solve Machado model:")
	machado = cmodel.generateCompactModel()
	cmodel.mipmodel = machado
	lp = machado.exportCPLEX()
	# lp.write("big.lp")
	sol = solveModel(cmodel, lp, 300, 250)
	instance.sigNorm = sol.signal
	print("Signal norm set to:", instance.sigNorm)
	print("Now solve bigModel:")
	# then generate bigModel
	mipmodel = cmodel.generateBigModel()
	cmodel.mipmodel = mipmodel
	lp = mipmodel.exportCPLEX()
	# lp.write("big.lp")
	return solveModel(cmodel, lp, 1500, 1200)


def getSolutionEx(cmodel, lp):
	mipmodel = cmodel.mipmodel
	sources = []
	detectors = []
	try:
		values = lp.solution.get_values()
		objval = lp.solution.get_objective_value()
	except:
		return NIRSSolution()
	x = mipmodel.getVars('x')
	y = mipmodel.getVars('y')
	for p in x:
		if values[x[p].modelindex] > 0.5:
			sources.append(p)
	for p in y:
		if values[y[p].modelindex] > 0.5:
			detectors.append(p)
	sol = NIRSSolution(sources, detectors, objval)
	# get coverage value
	w = mipmodel.getVars('w')
	cov = 0
	for v in w:
		if values[w[v].modelindex] > 0.5:
			cov += 1
	sol.coverage = cov
	# get signal value
	s = mipmodel.getVars('s')
	tot = 0.0
	for v in s:
		tot += values[s[v].modelindex]
	sol.signal = tot
	return sol


def filterTooClose(sources, elem, data, minDist):
	""" remove from list 'sources' all elements closer to 'elem' than 'minDist' """
	return [s for s in sources if data.getPosDistance(s, elem) > minDist]


def installSolution(cmodel, lp, solution):
	print("Installing solution:")
	solution.printQuality()
	mipmodel = cmodel.mipmodel
	# fix x values
	x = mipmodel.getVars('x')
	for p in solution.sources:
		lp.variables.set_lower_bounds(x[p].modelindex, 1.0)
	# fix y values
	y = mipmodel.getVars('y')
	for p in solution.detectors:
		lp.variables.set_lower_bounds(y[p].modelindex, 1.0)
	# quick solve
	lp.solve()
	# unfix x
	for p in solution.sources:
		lp.variables.set_lower_bounds(x[p].modelindex, x[p].lb)
	# unfix y
	for p in solution.detectors:
		lp.variables.set_lower_bounds(y[p].modelindex, y[p].lb)
	# quick solve
	lp.parameters.dettimelimit.set(1000)
	lp.solve()


def randomConstruction(cmodel, lp):
	data = cmodel.data
	nSources = cmodel.instance.nSources
	candidates = cmodel.instance.positions
	mipmodel = cmodel.mipmodel
	sources = []
	for i in range(nSources):
		# break on empty list
		if not candidates:
			break
		# pick element at random
		elem = random.choice(candidates)
		sources.append(elem)
		# shrink list
		candidates = filterTooClose(candidates, elem, data, cmodel.instance.minDist)
	print("Initial sources:", sources)
	# fix x values and optimize
	x = mipmodel.getVars('x')
	for p in sources:
		lp.variables.set_lower_bounds(x[p].modelindex, 1.0)
	lp.solve()
	# unfix x
	for p in x:
		lp.variables.set_lower_bounds(x[p].modelindex, x[p].lb)
	# quick solve
	lp.parameters.dettimelimit.set(1000)
	lp.solve()
	getSolutionEx(cmodel, lp).printQuality()


def findInitialSolution(cmodel, lp):
	lp.parameters.dettimelimit.set(10000)
	#lp.parameters.mip.strategy.search.set(1)
	lp.parameters.mip.tolerances.mipgap.set(0.01)
	lp.parameters.mip.limits.eachcutlimit.set(0)
	lp.parameters.mip.limits.cutpasses.set(-1)
	randomConstruction(cmodel, lp)


def alternateHeuristic(cmodel, lp):
	"""
	Fixes sources/detectors positions in turns.
	Leads to massive simplification of the model
	(all z_pq variables are aggregated out).
	However, if all sources (resp. detectors) are in the 'wrong' place,
	then there is no way for this heuristic to fix it!
	"""
	mipmodel = cmodel.mipmodel
	solution = getSolutionEx(cmodel, lp)
	lp.parameters.dettimelimit.set(10000)
	while True:
		lastobj = solution.objval
		# fix x values and reoptimize
		#print("####### Fix X")
		x = mipmodel.getVars('x')
		for p in solution.sources:
			lp.variables.set_lower_bounds(x[p].modelindex, 1.0)
		lp.solve()
		#lp.write("A.lp")
		solution = getSolutionEx(cmodel, lp)
		# unfix x
		for p in x:
			lp.variables.set_lower_bounds(x[p].modelindex, x[p].lb)
		# fix y values and reoptimize
		#print("####### Fix Y")
		y = mipmodel.getVars('y')
		for p in solution.detectors:
			lp.variables.set_lower_bounds(y[p].modelindex, 1.0)
		lp.solve()
		#lp.write("B.lp")
		solution = getSolutionEx(cmodel, lp)
		# unfix y
		for p in y:
			lp.variables.set_lower_bounds(y[p].modelindex, y[p].lb)
		solution.printQuality()
		if solution.objval <= lastobj + 1e-2:
			break
	lp.parameters.dettimelimit.set(1000)
	lp.solve()
	getSolutionEx(cmodel, lp).printQuality()


def localSearch(cmodel, lp, maxTries=5):
	"""
	Unfortunately, this is going to be rather expensive,
	as it does not trigger the massive model simplifications
	of the alternate heuristic.
	"""
	mipmodel = cmodel.mipmodel
	solution = getSolutionEx(cmodel, lp)
	lp.parameters.dettimelimit.set(20000)
	for k in range(maxTries):
		print("### Local Search Iteration")
		lastobj = solution.objval
		sources_to_fix = random.sample(solution.sources, len(solution.sources)-1)
		detectors_to_fix = random.sample(solution.detectors, len(solution.detectors)-1)
		# fix
		x = mipmodel.getVars('x')
		for p in sources_to_fix:
			lp.variables.set_lower_bounds(x[p].modelindex, 1.0)
		y = mipmodel.getVars('y')
		for p in detectors_to_fix:
			lp.variables.set_lower_bounds(y[p].modelindex, 1.0)
		# solve
		#lp.write("LS_%i.lp" % k)
		lp.solve()
		# get new solution
		solution = getSolutionEx(cmodel, lp)
		solution.printQuality()
		# unfix
		for p in sources_to_fix:
			lp.variables.set_lower_bounds(x[p].modelindex, x[p].lb)
		for p in detectors_to_fix:
			lp.variables.set_lower_bounds(y[p].modelindex, y[p].lb)
	lp.parameters.dettimelimit.set(1000)
	lp.solve()
	getSolutionEx(cmodel, lp).printQuality()


def solveMIP(cmodel, lp):
	lp.solve()
	status = lp.solution.get_status()
	return status == lp.solution.status.MIP_optimal or status == lp.solution.status.optimal_tolerance


def heuristicSearch(cmodel, lp, initialConstr=True, ls=False, maxDives=5):
	bestsol = getSolutionEx(cmodel, lp)
	print("### Initial solution:")
	bestsol.printQuality()
	for k in range(maxDives):
		print("################ Heuristic run:", k)
		if initialConstr:
			findInitialSolution(cmodel, lp)
		alternateHeuristic(cmodel, lp)
		if ls:
			localSearch(cmodel, lp)
		lp.parameters.dettimelimit.set(20000)
		opt = solveMIP(cmodel, lp)
		solution = getSolutionEx(cmodel, lp)
		if solution.objval > bestsol.objval + 1e-2:
			bestsol = solution
		if opt:
			break
	installSolution(cmodel, lp, bestsol)
	getSolutionEx(cmodel, lp).printQuality()


def solveBigIterative(data, instance):
	cmodel = NIRSModel(data, instance)
	mipmodel = cmodel.generateBigModel()
	cmodel.mipmodel = mipmodel
	print("#vars =", mipmodel.numVars)
	print("#constrs =", mipmodel.numConstraints)
	lp = mipmodel.exportCPLEX()
	# lp.write("combined.lp")
	# first quick solve
	# what about solving the Machado model to find a first solution?
	lp.parameters.dettimelimit.set(10000)
	opt = solveMIP(cmodel, lp)
	if opt is False:
		lp.set_results_stream(None)
		# heuristic solve first
		heuristicSearch(cmodel, lp)
		# truly final solve
		lp.set_results_stream(sys.stdout)
		lp.parameters.dettimelimit.set(10000000)
		lp.parameters.timelimit.set(1800)
		lp.parameters.mip.limits.eachcutlimit.set(1000000)
		lp.parameters.mip.limits.cutpasses.set(0)
		lp.solve()
	solution = getSolution(cmodel, lp)
	cmodel.instance.evalCoverageAndSignal(cmodel.data, solution)
	solution.dualBound = lp.solution.MIP.get_best_objective()
	return solution



if __name__ == '__main__':
	import sys
	from instances import *
	data = NIRSData()
	data.read(scalpFile, voxelFile, allPMDFfile, wFile)
	try:
		n = int(sys.argv[2])
	except:
		n = 8
	try:
		which = int(sys.argv[1])
		mixture = mixtures[which]
		instance = generateFromGaussianMixture(data, n, n, mixture)
	except:
		instance = readFromFile(data, n, n, sys.argv[1])
		try:
			instance.covWeight = float(sys.argv[3])
		except:
			instance.covWeight = 10.0
	random.seed(20070512)
	# solve
	solution = solveBigModel(data, instance)
	# print solution
	solution.printQuality()
	instance.plot(data, solution)
	instance.plotSolValues(data, solution)
