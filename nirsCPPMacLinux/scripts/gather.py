#!/usr/bin/env python

# gatherUtils.py
#
# @author Domenico Salvagnin <salvagni at math dot unipd dot it>
# 2017

from __future__ import print_function

import sys
sys.path.append('../python')
import xml.etree.cElementTree as et
import os.path
#import re

#from functools import partial
from batchUtils import locate
from table import Table


def getValue(doc, address, defValue):
	""" get the value of an XML node given its path (or defValue in case of error) """
	try:
		return doc.find(address).text
	except:
		return defValue


def parseDoc(filename, rows):
	""" parse a single log file, adding rows to rows """
	try:
		doc = et.parse(filename).getroot()
		common = {}
		common["problem"] = getValue(doc, './/instance', "unknown")
		common["nSources"] = int(getValue(doc, './/nSources', 0))
		common["nDetectors"] = int(getValue(doc, './/nDetectors', 0))
		common["minDist"] = float(getValue(doc, './/minDist', 0.0))
		resNode = doc.find('.//results')
		for algo in resNode.findall('algorithm'):
			local = {}
			local["method"] = getValue(algo, 'name', "unknown")
			local["time"] = float(getValue(algo, 'time', 0.0))
			local["coverage"] = float(getValue(algo, 'coverage', 0.0))
			local["signal"] = float(getValue(algo, 'signal', 0.0))
			local["primal"] = float(getValue(algo, 'primal', 0.0))
			local["sources"] = getValue(algo, 'sources', "[]")
			local["detectors"] = getValue(algo, 'detectors', "[]")
			row = common.copy()
			row.update(local)
			rows.append(row)
	except:
		print("Error parsing: %s" % filename)


def parseFolder(folder):
	""" parse all log files in a folder, returns table with results """
	t = Table()
	# add columns
	t.addColumn("problem")
	t.addColumn("nSources")
	t.addColumn("nDetectors")
	t.addColumn("minDist")
	t.addColumn("method")
	t.addColumn("time")
	t.addColumn("coverage")
	t.addColumn("signal")
	t.addColumn("primal")
	t.addColumn("sources")
	t.addColumn("detectors")
	# parse all matching logs in folder
	rows = []
	for f in locate("run.xml", folder):
		parseDoc(f, rows)
	t.setRows(rows)
	return t



if __name__ == '__main__':
	folder = sys.argv[1]
	csvfile = os.path.basename(folder) + ".csv"

	t = parseFolder(folder)
	t.sort(['problem', 'method', 'nSources', 'minDist'])
	t.write(csvfile)
	t.report()
