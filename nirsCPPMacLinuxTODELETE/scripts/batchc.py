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
import subprocess
import os.path

instances = [f for f in locate("ROI*.txt", "./data/down/")]
nSources = [2, 4, 8, 16]
minRhoOpt = [10, 15, 20, 25, 30]


if __name__ == '__main__':
	for f in instances:
		for n in nSources:
			for minDist in minRhoOpt:
				outputdir = "%s_%i_%i" % (os.path.basename(f)[:-4], n, minDist)
				subprocess.call("./nirsmain.opt -C Globals.outputDir=%s %s %i %i" % (outputdir, f, n, minDist), shell=True)
