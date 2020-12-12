#!/usr/bin/env python
# encoding: utf-8

import sys
sys.path.append('../python')
import math
import numpy as np
import random
import itertools
import model


def powerset(n):
	s = list(range(n))
	return list(itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1)))


def pointsPorta(n,m):
	left = powerset(n)
	right = powerset(m)
	print "DIM =", n + m + n*m
	print "CONV_SECTION"
	for x in left:
		xind = np.zeros(n, dtype=np.int8)
		for elem in x:
			xind[elem] = 1
		for y in right:
			yind = np.zeros(m, dtype=np.int8)
			for elem in y:
				yind[elem] = 1
			zind = np.zeros((n,m), dtype=np.int8)
			for i in x:
				for j in y:
					zind[i,j] = 1
			zind = zind.flatten()
			point = np.concatenate((xind,yind,zind))
			print point
	print "END"


def basicPortaFormulation(n,m):
	dim = n + m + n*m
	print "DIM =", dim
	print "VALID"
	print ("0 " * dim).strip()
	print "LOWER_BOUNDS"
	print ("0 " * dim).strip()
	print "UPPER_BOUNDS"
	print ("1 " * dim).strip()
	print "INEQUALITIES_SECTION"
	# bounds
	for i in range(1,dim+1):
		print "x%i >= 0" % i
		print "x%i <= 1" % i
	# McCormick
	for i in range(n):
		for j in range(m):
			xind = i+1
			yind = n+j+1
			zind = n+m+i*n+j+1
			print "x%i - x%i >= 0" % (xind,zind)
			print "x%i - x%i >= 0" % (yind,zind)
			print "x%i + x%i - x%i <= 1" % (xind,yind,zind)
	print "END"


def makeRandomBBQP(n, m):
	mip = model.Model()
	mip.minimize = True
	xpos = range(1,n+1)
	ypos = range(1,m+1)
	# x,y vars
	x = {}
	for p in xpos:
		x[p] = model.Variable("x_%i" % p, 'B', obj=0.0, lb=0.0, ub=1.0, index=p)
	y = {}
	for q in ypos:
		y[q] = model.Variable("y_%i" % q, 'B', obj=0.0, lb=0.0, ub=1.0, index=q)
	mip.addVars(x, 'x')
	mip.addVars(y, 'y')
	# z vars
	z = {}
	for p in xpos:
		for q in ypos:
			#if p == q:
			#	continue
			z[(p,q)] = model.Variable("z_%i_%i" % (p,q), 'B', obj=random.uniform(-10,10), lb=0.0, ub=1.0, index=(p,q))
	mip.addVars(z, 'z')
	# linking z to x,y
	for p in xpos:
		for q in ypos:
			#if p == q:
			#	continue
			zpq = z[(p,q)]
			xp = x[p]
			yq = y[q]
			c1 = model.Constraint('zx_%i_%i' % (p,q), 'L', 0, index=(p,q))
			c1.row.append((zpq,1))
			c1.row.append((xp,-1))
			mip.addConstraint(c1)
			c2 = model.Constraint('zy_%i_%i' % (p,q), 'L', 0, index=(p,q))
			c2.row.append((zpq,1))
			c2.row.append((yq,-1))
			mip.addConstraint(c2)
			c3 = model.Constraint('zxy_%i_%i' % (p,q), 'L', 1, index=(p,q))
			c3.row.append((xp,1))
			c3.row.append((yq,1))
			c3.row.append((zpq,-1))
			mip.addConstraint(c3)
	return mip

if __name__ == '__main__':
	n = int(sys.argv[1])
	m = int(sys.argv[2])
	#pointsPorta(n,m)
	#basicPortaFormulation(n,m)
	m = makeRandomBBQP(n,m)
	m.output()
