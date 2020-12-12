#!/usr/bin/env python
"""
models.py

NIRS optimization models
@author Domenico Salvagnin dominiqs@gmail.com
"""

import sys
sys.path.append('../python')

import math
import random
import h5py
import numpy as np
from nirs import *
import model

eps = 1e-9
bigM = 1000.0

class NIRSModel(object):
	"""docstring for NIRSModel"""
	def __init__(self, data, instance):
		super(NIRSModel, self).__init__()
		self.data = data
		self.instance = instance
		self.mipmodel = None

	def createXY(self, mipmodel):
		positions = self.instance.positions
		x = {}
		y = {}
		for p in positions:
			x[p] = model.Variable("x_%i" % p, 'B', obj=0.0, lb=0.0, ub=1.0, index=p)
			y[p] = model.Variable("y_%i" % p, 'B', obj=0.0, lb=0.0, ub=1.0, index=p)
		mipmodel.addVars(x, 'x')
		mipmodel.addVars(y, 'y')
		return (x,y)

	def basicConstraints(self, x, y, mipmodel):
		positions = self.instance.positions
		sumsources = model.Constraint('sumsources', 'L', self.instance.nSources)
		sumdetectors = model.Constraint('sumdetectors', 'L', self.instance.nDetectors)
		for p in positions:
			sumsources.row.append((x[p], 1))
			sumdetectors.row.append((y[p], 1))
			onlyone = model.Constraint('onlyone_%i' % p, 'L', 1, index=p)
			onlyone.row.append((x[p], 1))
			onlyone.row.append((y[p], 1))
			mipmodel.addConstraint(onlyone)
		mipmodel.addConstraint(sumsources)
		mipmodel.addConstraint(sumdetectors)

	def distanceConstraints(self, x, y, z, mipmodel):
		minRho = self.instance.minRho
		maxRho = self.instance.maxRho
		minDist = self.instance.minDist
		for p in self.instance.positions:
			for q in self.instance.positions:
				if p >= q:
					continue
				dist = self.data.getPosDistance(p,q)
				# source-detectors distance constraints
				if dist < max(minRho,minDist) or dist > maxRho:
					incxy = model.Constraint('incxy_%i_%i' % (p,q), 'L', 1, index=(p,q))
					incxy.row.append((x[p], 1))
					incxy.row.append((y[q], 1))
					mipmodel.addConstraint(incxy)
					incyx = model.Constraint('incyx_%i_%i' % (p,q), 'L', 1, index=(p,q))
					incyx.row.append((y[p], 1))
					incyx.row.append((x[q], 1))
					mipmodel.addConstraint(incyx)
					# we can set the corresponding z[(p,q)] to 0
					# This is implied by x_p + y_q <= 1, but it strenghten the LP relaxation
					if z is not None:
						z[(p,q)].ub = 0.0
						z[(q,p)].ub = 0.0
				# any two optodes
				if dist < minDist:
					incxx = model.Constraint('incxx_%i_%i' % (p,q), 'L', 1, index=(p,q))
					incxx.row.append((x[p], 1))
					incxx.row.append((x[q], 1))
					mipmodel.addConstraint(incxx)
					incyy = model.Constraint('incyy_%i_%i' % (p,q), 'L', 1, index=(p,q))
					incyy.row.append((y[p], 1))
					incyy.row.append((y[q], 1))
					mipmodel.addConstraint(incyy)

	def generateOldCompactModel(self):
		# original version of the Machado model
		aggr = self.instance.computeAggregate(self.data)
		mipmodel = model.Model()
		mipmodel.minimize = False
		positions = self.instance.positions
		# declare variables
		x,y = self.createXY(mipmodel)
		w = {}
		for p in positions:
			w[p] = model.Variable("w_%i" % p, 'C', obj=1.0, index=p)
		mipmodel.addVars(w, 'w')
		# add basic constraints
		self.basicConstraints(x, y, mipmodel)
		# bigM constraints
		for p in positions:
			vub = model.Constraint('vub_%i' % p, 'L', 0, index=p)
			vub.row.append((w[p], 1))
			vub.row.append((y[p], -bigM))
			mipmodel.addConstraint(vub)
		# constraints that define w variables
		for q in positions:
			wdef = model.Constraint('wdef_%i' % q, 'L', 0, index=q)
			wdef.row.append((w[q], 1))
			for p in positions:
				coef = aggr[(p,q)]
				if coef > eps:
					wdef.row.append((x[p], -coef))
			mipmodel.addConstraint(wdef)
		# distance constraints
		self.distanceConstraints(x, y, None, mipmodel)
		return mipmodel

	def generateCompactModel(self):
		# improved version of the Machado model
		aggr = self.instance.computeAggregate(self.data)
		mipmodel = model.Model()
		mipmodel.minimize = False
		positions = self.instance.positions
		# declare variables
		x,y = self.createXY(mipmodel)
		z = self.createZ(mipmodel)
		# add basic constraints
		self.basicConstraints(x, y, mipmodel)
		# define z
		self.zDef(x, y, z, mipmodel, usercuts=True)
		# distance constraints
		self.distanceConstraints(x, y, z, mipmodel)
		# objective
		for p in positions:
			for q in positions:
				if p == q:
					continue
				objcoef = aggr[(p,q)]
				if objcoef > eps:
					z[(p,q)].obj = objcoef
		return mipmodel

	def createZ(self, mipmodel):
		positions = self.instance.positions
		z = {}
		for p in positions:
			for q in positions:
				if p == q:
					continue
				z[(p,q)] = model.Variable("z_%i_%i" % (p,q), 'B', obj=0.0, lb=0.0, ub=1.0, index=(p,q))
		mipmodel.addVars(z, 'z')
		return z

	def createS(self, mipmodel):
		voxels = self.instance.voxels
		s = {}
		objval = 1.0 / self.instance.sigNorm
		for v in voxels:
			s[v] = model.Variable("s_%i" % v, 'C', obj=objval, lb=0.0, ub=None, index=v)
		mipmodel.addVars(s, 's')
		return s

	def setObjCoef(self, vars, value):
		for v in vars:
			vars[v].obj = value

	def zDef(self, x, y, z, mipmodel, usercuts=False):
		positions = self.instance.positions
		for p in positions:
			xp = x[p]
			c = model.Constraint('zsumq_%i' % p, 'L', 0, index=p)
			for q in positions:
				if p == q:
					continue
				zpq = z[(p,q)]
				c.row.append((zpq,1))
			c.row.append((xp,-self.instance.nDetectors))
			mipmodel.addConstraint(c)
		for q in positions:
			yq = y[q]
			c = model.Constraint('zsump_%i' % q, 'L', 0, index=q)
			for p in positions:
				if p == q:
					continue
				zpq = z[(p,q)]
				c.row.append((zpq,1))
			c.row.append((yq,-self.instance.nSources))
			mipmodel.addConstraint(c)
		if usercuts:
			for p in positions:
				for q in positions:
					if p == q:
						continue
					zpq = z[(p,q)]
					xp = x[p]
					yq = y[q]
					c1 = model.Constraint('zx_%i_%i' % (p,q), 'L', 0, index=(p,q))
					c1.row.append((zpq,1))
					c1.row.append((xp,-1))
					mipmodel.addConstraint(c1, cut=True)
					c2 = model.Constraint('zy_%i_%i' % (p,q), 'L', 0, index=(p,q))
					c2.row.append((zpq,1))
					c2.row.append((yq,-1))
					mipmodel.addConstraint(c2, cut=True)

	def signalDef(self, s, z, mipmodel, sense='G'):
		voxels = [s[v].index for v in s]
		positions = self.instance.positions
		# create empty constraints
		sdef = {}
		for v in voxels:
			c = model.Constraint('sdef_%i' % v, sense, 0, index=v)
			c.row.append((s[v], -1))
			sdef[v] = c
		# fill up coefficients
		for p in positions:
			for q in positions:
				if p >= q:
					continue
				# if the pair is not in the PMDF file, then it means that it is all zeros
				if (p,q) not in self.data.pmdfpairs:
					continue
				weight = 1.0
				if self.data.weights is not None:
					weight = self.data.weights[p,q]
				# same if the weight is zero
				if weight <= eps:
					continue
				pmdf = self.data.pmdfsource[self.data.pmdfpairs[(p,q)]]
				arr = np.array(pmdf)
				for v in voxels:
					val = weight * arr[0,v]
					if abs(val) > eps:
						sdef[v].row.append((z[(p,q)], val))
						sdef[v].row.append((z[(q,p)], val))
		# add constraints to model
		for v in voxels:
			mipmodel.addConstraint(sdef[v])

	def createW(self, mipmodel):
		voxels = self.instance.voxels
		w = {}
		objval = self.instance.covWeight / self.instance.covNorm
		for v in voxels:
			w[v] = model.Variable("w_%i" % v, 'B', obj=objval, lb=0.0, ub=1.0, index=v)
		mipmodel.addVars(w, 'w')
		return w

	def linkSW(self, s, w, mipmodel):
		voxels = self.instance.voxels
		for v in voxels:
			c = model.Constraint('link_sw_%i' % v, 'G', 0, index=v)
			c.row.append((s[v], 1))
			c.row.append((w[v], -self.instance.covThr))
			mipmodel.addConstraint(c)

	def generateBigModel(self):
		mipmodel = model.Model()
		mipmodel.minimize = False
		# declare variables
		x,y = self.createXY(mipmodel)
		z = self.createZ(mipmodel)
		s = self.createS(mipmodel)
		w = self.createW(mipmodel)
		# add basic constraints
		self.basicConstraints(x, y, mipmodel)
		# constraints linking z and x,y
		self.zDef(x, y, z, mipmodel, True)
		# constraints that define s[v]
		self.signalDef(s, z, mipmodel, 'E')
		self.linkSW(s, w, mipmodel)
		# distance constraints
		self.distanceConstraints(x, y, z, mipmodel)
		return mipmodel



if __name__ == '__main__':
	import sys
	from instances import *
	data = NIRSData()
	data.read(scalpFile, voxelFile, allPMDFfile, wFile)
	try:
		which = int(sys.argv[1])
		mixture = mixtures[which]
		instance = generateFromGaussianMixture(data, 5, 5, mixture)
	except:
		instance = readFromFile(data, 5, 5, sys.argv[1])
	cmodel = NIRSModel(data, instance)
	cmodel.generateBigModel().exportCPLEX().write(sys.argv[2])
