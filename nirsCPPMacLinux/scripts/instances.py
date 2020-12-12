#!/usr/bin/env python
"""
instances.py

NIRS toy instances
@author Domenico Salvagnin dominiqs@gmail.com
"""

import nirs

scalpFile   = "data/down2/scalp.h5"
voxelFile   = "data/down2/voxels.h5"
allPMDFfile = "data/down2/allpmdfs.h5"
wFile = "data/PMDF_MASK_30-60mmWeight_ToastCorrect.txt"


mixtures = [nirs.GaussianMixture(([-38.0, -33.0, 71.0], ), 0.01, 0.02),
	nirs.GaussianMixture(([-38.0, -33.0, 71.0], ), 0.01, 0.02),
	nirs.GaussianMixture(([-38.0, -33.0, 71.0], [-18.0, -48.0, 78.0], [-53.0, -23.0, 59.0], [-65.0, -12.0, 37.0]), 0.01, 0.08),
	nirs.GaussianMixture(([-38.0, -33.0, 71.0], [38.0, -33.0, 71.0]), 0.025, 0.1)]
