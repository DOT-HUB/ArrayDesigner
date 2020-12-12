#!/usr/bin/env python
"""
nirs.py

NIRS optimization data
@author Domenico Salvagnin dominiqs@gmail.com
"""

from __future__ import print_function
import math
import h5py
from itertools import product, combinations
import numpy as np
from scipy.linalg import norm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def drawCube(ax, center, size, color):
	r = [-1, 1]
	for s, e in combinations(np.array(list(product(r,r,r))), 2):
		if np.sum(np.abs(s-e)) == r[1]-r[0]:
			ax.plot3D(*zip(center + 0.5 * size * s, center + 0.5 * size * e), color=color)


def getDistance(pt1, pt2):
	return np.linalg.norm(pt1-pt2)


class NIRSData(object):
	"""Global data for NIRS"""
	def __init__(self):
		super(NIRSData, self).__init__()
		self.positions = None
		self.voxels = None
		self.nPositions = 0
		self.nVoxels = 0
		self.pmdfsource = None
		self.pmdfpairs = {}
		self.weights = None

	def read(self, posfile, voxfile, pmdffile, weights=None):
		f = h5py.File(posfile, "r")
		positions = f['posScalp']
		self.positions = np.array(positions)
		print("positions:", self.positions.shape)
		f = h5py.File(voxfile, "r")
		voxels = f['posGM']
		self.voxels = np.array(voxels)
		print("voxels:", self.voxels.shape)
		self.pmdfsource = h5py.File(pmdffile, "r")
		self.nPositions = len(self.positions[0])
		self.nVoxels = len(self.voxels[0])
		self.pmdfpairs = {}
		for k in self.pmdfsource.keys():
			tokens = k.split('_')
			row = int(tokens[2])
			col = int(tokens[4])
			self.pmdfpairs[(row-1,col-1)] = k
		if weights:
			self.weights = np.loadtxt(weights)
			# make it symmetric (in the file only the upper triangular part is nonzero)
			#upper = np.triu(self.weights)
			#self.weights += upper.transpose()
			print("weights:", self.weights.shape)

	def getPosDistance(self, p, q):
		pos1 = self.positions[:,p]
		pos2 = self.positions[:,q]
		return getDistance(pos1, pos2)

	def plotSets(self):
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		ax.scatter(self.positions[0], self.positions[1], self.positions[2], c='r')
		ax.scatter(self.voxels[0], self.voxels[1], self.voxels[2], c='b')
		plt.show()


class NIRSSolution(object):
	"""docstring for Instance"""
	def __init__(self, sources=[], detectors=[], objval=-1e20):
		super(NIRSSolution, self).__init__()
		self.sources = sources
		self.detectors = detectors
		self.objval = objval
		self.coverage = 0.0
		self.signal = 0.0
		self.dualBound = None
	def printQuality(self):
		print(">>> Solution: sources=%s detectors=%s coverage=%i signal=%g primal=%g dual=%s"
			% (self.sources, self.detectors, self.coverage, self.signal, self.objval, str(self.dualBound)))


class GaussianMixture(object):
	"""Data needed to generate a NIRS instance in a gaussian mixture model"""
	def __init__(self, roiCenters, gamma, threshold):
		super(GaussianMixture, self).__init__()
		self.roiCenters = roiCenters
		self.gamma = gamma
		self.threshold = threshold
		self.maxRho = 45 #mm
		self.gmDepth = 15 #mm


class Instance(object):
	"""docstring for Instance"""
	def __init__(self, nSources, nDetectors):
		super(Instance, self).__init__()
		self.nSources = nSources
		self.nDetectors = nDetectors
		self.positions = None
		self.voxels = None
		self.minRho = 15 # minimum distance between a source-detector pair (in mm)
		self.maxRho = 1000 # maximum distance between a source-detector pair (in mm)
		self.minDist = 10 # minimum distance between any two optodes (in mm)
		self.covNorm = 1.0 # denominator to normalize coverage
		self.sigNorm = 1.0 # denominator to normalize signal
		self.covWeight = 10.0 # weight of coverage in the objective
		self.covThr = 0.1528 # hard coded threshold for the current mesh

	def computeAggregate(self, data):
		aggr = {}
		for p in self.positions:
			for q in self.positions:
				if (p >= q):
					continue
				if (p,q) not in data.pmdfpairs:
					aggr[(p,q)] = 0.0
					aggr[(q,p)] = 0.0
					continue
				weight = 1.0
				if data.weights is not None:
					weight = data.weights[p,q]
				pmdf = data.pmdfsource[data.pmdfpairs[(p,q)]]
				arr = np.array(pmdf)
				val = weight * np.sum(arr[:,self.voxels])
				aggr[(p,q)] = val
				aggr[(q,p)] = val
			aggr[(p,p)] = 0
		return aggr

	def evalCoverageAndSignal(self, data, solution):
		sources = solution.sources
		detectors = solution.detectors
		if not sources or not detectors:
			solution.signal = 0.0
			solution.coverage = 0.0
			return
		# compute signal in each voxel of the ROI
		cum = np.zeros((1,data.nVoxels))
		for s in sources:
			for d in detectors:
				pmdf = None
				weight = 1.0
				if (s,d) in data.pmdfpairs:
					pmdf = data.pmdfsource[data.pmdfpairs[(s,d)]]
					if data.weights is not None:
						weight = data.weights[s,d]
				elif (d,s) in data.pmdfpairs:
					pmdf = data.pmdfsource[data.pmdfpairs[(d,s)]]
					if data.weights is not None:
						weight = data.weights[s,d]
				if pmdf is not None:
					cum += weight * np.array(pmdf)
		projected = cum[:,self.voxels]
		# total signal
		solution.signal = np.sum(projected)
		# coverage
		cov = 0
		for v in self.voxels:
			sigv = cum[0,v]
			if sigv >= self.covThr -1e-6:
				cov += 1
		solution.coverage = cov


	def plot(self, data, solution=None):
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		if self.voxels:
			voxels = data.voxels[:,self.voxels]
			ax.scatter(voxels[0], voxels[1], voxels[2], c='r', alpha=0.6, linewidth=0)
		if self.positions:
			positions = data.positions[:,self.positions]
			ax.scatter(positions[0], positions[1], positions[2], c='b', alpha=0.6, linewidth=0)
		if solution:
			sources = solution.sources
			detectors = solution.detectors
			for s in sources:
				drawCube(ax, data.positions[:,s], 5, 'k')
			for d in detectors:
				drawCube(ax, data.positions[:,d], 5, 'g')
		plt.show()

	def plotSolValues(self, data, solution):
		sources = solution.sources
		detectors = solution.detectors
		if not sources or not detectors:
			return
		# compute signal in each voxel of the ROI
		cum = np.zeros((1,data.nVoxels))
		for s in sources:
			for d in detectors:
				pmdf = None
				weight = 1.0
				if (s,d) in data.pmdfpairs:
					pmdf = data.pmdfsource[data.pmdfpairs[(s,d)]]
					if data.weights is not None:
						weight = data.weights[s,d]
				elif (d,s) in data.pmdfpairs:
					pmdf = data.pmdfsource[data.pmdfpairs[(d,s)]]
					if data.weights is not None:
						weight = data.weights[s,d]
				if pmdf is not None:
					cum += weight * np.array(pmdf)
		cum = cum[:,self.voxels]
		print("Objective value: coverage=%f signal=%f" % (solution.coverage,solution.signal))
		ax1 = plt.subplot2grid((2,2), (0,0))
		ax2 = plt.subplot2grid((2,2), (0,1))
		ax3 = plt.subplot2grid((2,2), (1,0), colspan=2)
		ax1.boxplot(cum)
		ax2.bar(range(len(self.voxels)), np.sort(sum(cum)), log=True)
		ax3.hist(cum.transpose(), 20)
		plt.show()


def generateFromGaussianMixture(data, nSources, nDetectors, mixture):
	instance = Instance(nSources, nDetectors)
	# select voxels
	nVoxels = data.nVoxels
	nCenters = len(mixture.roiCenters)
	distances = np.zeros((nVoxels,nCenters))
	instance.voxels = []
	for v in range(nVoxels):
		voxel = data.voxels[:,v]
		roi = 0.0
		for c in range(nCenters):
			center = mixture.roiCenters[c]
			dist = math.pow(getDistance(voxel, center), 2)
			roi += math.exp(-mixture.gamma * dist)
		if roi > mixture.threshold:
			instance.voxels.append(v)
	print("Selected %i voxels in ROI" % len(instance.voxels))
	# select positions
	hyp = math.sqrt(mixture.maxRho*mixture.maxRho + mixture.gmDepth*mixture.gmDepth)
	voxels = data.voxels[:,instance.voxels]
	nPositions = data.nPositions
	instance.positions = []
	for p in range(nPositions):
		position = data.positions[:,p]
		mindist = 1e20
		for v in instance.voxels:
			voxel = data.voxels[:,v]
			dist = getDistance(voxel, position)
			mindist = min(mindist, dist)
		if mindist < hyp:
			instance.positions.append(p)
	print("Selected %i possible positions" % len(instance.positions))
	return instance


def readFromFile(data, nSources, nDetectors, filename):
	instance = Instance(nSources, nDetectors)
	fp = open(filename, 'r')
	# first line contains the voxels in the ROI
	instance.voxels = [int(x)-1 for x in fp.readline().split()]
	print("Selected %i voxels in ROI" % len(instance.voxels))
	# second line contains the possible positions
	instance.positions = [int(x)-1 for x in fp.readline().split()]
	print("Selected %i possible positions" % len(instance.positions))
	# compute normalizations
	instance.covNorm = float(len(instance.voxels))
	instance.sigNorm = 1.0
	return instance



if __name__ == '__main__':
	import sys
	from instances import *
	data = NIRSData()
	data.read(scalpFile, voxelFile, allPMDFfile, wFile)
	instanceFile = sys.argv[1]
	n = int(sys.argv[2])
	instance = readFromFile(data, n, n, sys.argv[1])
	try:
		solution = NIRSSolution(eval(sys.argv[3]), eval(sys.argv[4]))
		instance.evalCoverageAndSignal(data, solution)
		solution.printQuality()
	except:
		solution = None
	# plot
	instance.plot(data, solution)
	instance.plotSolValues(data, solution)
