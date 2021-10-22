/**
 * @file nirs.cpp
 * @brief NIRS Data Structures
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2019 Domenico Salvagnin
 */

#include "nirs/nirs.h"
#include "nirs/asserter.h"
#include "nirs/str_utils.h"
#include <iostream>
#include <fstream>
#include <numeric>

using namespace dominiqs;


static void readPoint3dVector(const std::string& file, std::vector<Point3D>& v)
{
	std::ifstream in(file);
	double x,y,z;
	while(in >> x >> y >> z)
	{
		v.emplace_back(x, y, z);
	}
}


/* Some comments on the pmdf file format.
 *
 * For each pair of optodes positions (p,q), we need the signal val at each voxel v.
 * This is very sparse, and we exploit this to optimize storage and read times.
 * The format is split into two files: an index file and a data file.
 *
 * The index file (a text file) gives the number of voxels with a nonzero signal
 * for each pair (p,q) with p < q. In other words, it is a list of [n(n-1)/2 - n]
 * non-negative integers, where n is the number of possible positions for an optode,
 * sorted by p and then by q (basically an upper triangular matrix in row-major).
 *
 * The data file (a binary file) is basically a long stream of pairs (v,val),
 * where v is a signed 32-bit integer, and val is a double (64-bit).
 *
 * Example for n=3 and 4 values in total
 *
 * Index file:
 * 2 1
 * 1
 *
 * Data file (binary file of size 48 bytes):
 * |v1|val1|v2|val2|v3|val3|v4|val4|
 *
 * Where v1...val2 gives the values for the pair (0,1), v3..val3 for the pair (0,2)
 * and v4..val4 for the pair (1,2).
 */

void NIRSData::read(const std::string& posfile,
                    const std::string& voxfile,
                    const std::string& pmdffileidx,
                    const std::string& pmdffiledata,
                    const std::string& wfile)
{
	std::cout << "NIRSDATA:" << std::endl;
	// read positions from HDF5 file
	DOMINIQS_ASSERT(posfile.size());
	readPoint3dVector(posfile, positions);
	std::cout << "# Scalp positions = " << positions.size() << std::endl;

	// read voxels from HDF5 file
	DOMINIQS_ASSERT(voxfile.size());
	readPoint3dVector(voxfile, voxels);
	std::cout << "# GM nodes = " << voxels.size() << std::endl;

	// read all PMDFs
	size_t numpositions = positions.size();
	allpmdfs.clear();
	allpmdfs.resize(numpositions * numpositions);
	std::ifstream inidx(pmdffileidx);
	std::vector<int> counts(numpositions * numpositions, 0);
	size_t totcnt = 0;
	for (size_t p = 0; p < numpositions; p++)
	{
		for (size_t q = p+1; q < numpositions; q++)
		{
			inidx >> counts[p*numpositions+q];
			totcnt += counts[p*numpositions+q];
		}
	}
	std::ifstream indata(pmdffiledata, std::ios::binary);
	size_t numvalues = 0;
	for (size_t p = 0; p < numpositions; p++)
	{
		for (size_t q = p+1; q < numpositions; q++)
		{
			int cnt = counts[p*numpositions+q];
			if (!cnt)  continue;
			PMDFPtr pmdf = std::make_shared<PMDF>(p,q);
			for (int k = 0; k < cnt; k++)
			{
				int32_t vox;
				double val;
				indata.read((char*)&vox, sizeof(vox));
				indata.read((char*)&val, sizeof(val));
				pmdf->data.emplace_back(vox, val);
				numvalues++;
			}
			allpmdfs[p*numpositions+q] = pmdf;
		}
	}
	std::cout << "# AllPMDF values = " << numvalues << std::endl;
	DOMINIQS_ASSERT( numvalues == totcnt );
	// read weights
	weights.clear();
	weights.reserve(numpositions * numpositions);
	if (wfile.size())
	{
		std::ifstream inw(wfile);
		while (inw)
		{
			double w;
			inw >> w;
			weights.push_back(w);
		}
		DOMINIQS_ASSERT(weights.size() >= (numpositions * numpositions));
		// on input the weight matrix might be upper triangular: we have to make symmetric now
		for (size_t i = 0; i < numpositions; i++)
		{
			for (size_t j = 0; j < i; j++)
			{
				weights[i*numpositions+j] = weights[j*numpositions+i];
			}
		}
	}
	else
	{
		weights.resize(numpositions * numpositions, 1.0);
	}
}


Instance::~Instance()
{
	pmdfpairs.clear();
}


void Instance::read(size_t ns, size_t nd, const std::string& filename)
{
	std::cout << "INSTANCE:" << std::endl;
	nSources = ns;
	nDetectors = nd;
	std::cout << "# sources = " << nSources << std::endl;
	std::cout << "# detectors = " << nDetectors << std::endl;
	std::ifstream in(filename);

	// first line contains the voxels in the ROI
	std::string line;
	std::getline(in, line);
	voxels = split<int>(line, " ");
	// subtract -1 from all values (again, Matlab indexing)
	std::transform(voxels.begin(), voxels.end(), voxels.begin(), [](int x) { return x-1; });
	std::cout << "# voxels = " << voxels.size() << std::endl;
	// second line contains the possible positions
	std::getline(in, line);
	positions = split<int>(line, " ");
	// subtract -1 from all values (again, Matlab indexing)
	std::transform(positions.begin(), positions.end(), positions.begin(), [](int x) { return x-1; });
	std::cout << "# positions = " << positions.size() << std::endl;

	// compute matrix of pair-wise positions
	size_t numpositions = gdata.positions.size();
	distances.resize(numpositions * numpositions);
	std::fill(distances.begin(), distances.end(), 0.0);
	for (int p: positions)
	{
		for (int q: positions)
		{
			distances[p*numpositions+q] = getDistance(gdata.positions[p], gdata.positions[q]);
		}
	}

	// get only PMDFs of interest, scaling signal by weight
	copyPMDFs();

	// fill sparse representation (basically a bucket sort)
	std::vector<size_t> counts(numpositions, 0);
	start.resize(numpositions+1);
	std::fill(start.begin(), start.end(), 0);
	for (PMDFPtr pmdf: pmdfpairs)
	{
		if (!pmdf) continue;
		counts[pmdf->p]++;
	}
	size_t tot = 0;
	for (int p = 0; p < numpositions; p++)
	{
		start[p] = tot;
		tot += counts[p];
	}
	start[numpositions] = tot;
	pmdfs.resize(tot);
	for (PMDFPtr pmdf: pmdfpairs)
	{
		if (!pmdf) continue;
		int p = pmdf->p;
		pmdfs[start[p] + (--counts[p])] = pmdf;
	}

	// compute normalizations
	covNorm = double(voxels.size());
	sigNorm = 1.0;
}



PMDFPtr Instance::getPMDF(int p, int q) const
{
	return pmdfpairs[p*offset+q];
}


double Instance::getDist(int p, int q) const
{
	return distances[p*offset+q];
}


void Instance::copyPMDFs()
{
	size_t numpositions = gdata.positions.size();
	pmdfpairs.clear();
	pmdfpairs.resize(offset * offset);
	std::vector<double> weighted(gdata.voxels.size(), 0.0);

	int numpairs = 0;
	for (int p: positions)
	{
		for (int q: positions)
		{
			if (p >= q) continue;
			if (!gdata.allpmdfs[p*offset+q])  continue;

			PMDFPtr origpmdf = gdata.allpmdfs[p*offset+q];
			DOMINIQS_ASSERT( origpmdf->p == p );
			DOMINIQS_ASSERT( origpmdf->q == q );

			// scatter pmdf->data into weighted
			for (const auto& item: origpmdf->data)  weighted[item.voxel] = item.value;

			// create copy of PMDF, with only the voxels in the ROI of the current instance,
			// and values multiplied by the current weight
			PMDFPtr pmdf = std::make_shared<PMDF>(p, q);
			numpairs++;
			// note: we only look at the voxels in the ROI defined by the current instance
			pmdf->sum = 0.0;
			for (int v: voxels)
			{
				// the value stored is already multiplied by the corresponding weight
				// in addition, small values are already discarded
				double coeff = weighted[v];
				coeff *= gdata.weights[p*numpositions+q];
				if (coeff >= 1e-9)
				{
					pmdf->data.emplace_back(v, coeff);
					pmdf->sum += coeff;
				}
			}
			// Set a PMDF for both (p,q) and (q,p)
			pmdfpairs[p*offset+q] = pmdf;
			pmdfpairs[q*offset+p] = pmdf;

			// reset weighted
			for (const auto& item: origpmdf->data)  weighted[item.voxel] = 0.0;
		}
	}
	std::cout << "# PMDFs = " << numpairs << std::endl;
}


void Instance::evalObjectives(const std::vector<double>& cum, double& cov, double& sig, double& obj) const
{
	// compute total signal and coverage
	double _sig = 0.0;
	size_t _cov = 0;
	for (auto v: voxels)
	{
		double value = cum[v];
		double thr = coverageThr;
		_sig += value;
		if (value >= thr - 1e-6) _cov++;
	}
	sig = _sig;
	cov = double(_cov);
	// compute weighted normalized objective
	obj = (sig / sigNorm) + coverageWeight * (cov / covNorm);
}


std::ostream& operator<<(std::ostream& os, const Point3D& p)
{
	os << "(" << p.x() << "," << p.y() << "," << p.z() << ")";
	return os;
}


BSolution& BSolution::operator=(const BSolution& other)
{
	if (this != &other)
	{
		DOMINIQS_ASSERT(&instance == &(other.instance));
		sources = other.sources;
		detectors = other.detectors;
		coverage = other.coverage;
		signal = other.signal;
		objective = other.objective;
	}
	return *this;
}


void BSolution::reset()
{
	sources.clear();
	detectors.clear();
	coverage = 0.0;
	signal = 0.0;
	objective = 0.0;
}


void BSolution::computeSignal(std::vector<double>& cum) const
{
	std::fill(cum.begin(), cum.end(), 0.0);
	for (auto s: sources)
	{
		if (s == -1) continue;
		for (auto d: detectors)
		{
			if (d == -1) continue;
			PMDFPtr pmdf = instance.getPMDF(s, d);
			if (!pmdf) continue;
			for (const PMDF::Item& item: pmdf->data) cum[item.voxel] += item.value;
		}
	}
}


void BSolution::eval()
{
	coverage = 0.0;
	signal = 0.0;
	objective = 0.0;
	if (sources.empty() || detectors.empty()) return;

	// compute total signal in each voxel
	std::vector<double> cum(instance.gdata.voxels.size(), 0.0);
	computeSignal(cum);

	// evaluate objective (that depends on the ROI)
	instance.evalObjectives(cum, coverage, signal, objective);
}


bool BSolution::checkSignal(const std::vector<double>& cum) const
{
	std::vector<double> cumok(instance.gdata.voxels.size(), 0.0);
	computeSignal(cumok);
	for (auto v: instance.voxels)
	{
		if (fabs(cum[v]-cumok[v]) > 1e-9)
		{
			std::cout << "Check failed for voxel " << v << ": " << cum[v] << " != " << cumok[v] << std::endl;
			return false;
		}
	}
	return true;
}


bool BSolution::checkCompatibilities() const
{
	double thr = instance.minDist;
	for (auto p1: sources)
	{
		for (auto p2: sources)
		{
			if (p1 == p2) continue;
			if (getDistance(instance.gdata.positions[p1], instance.gdata.positions[p2]) < thr) return false;
		}
	}
	for (auto q1: detectors)
	{
		for (auto q2: detectors)
		{
			if (q1 == q2) continue;
			if (getDistance(instance.gdata.positions[q1], instance.gdata.positions[q2]) < thr) return false;
		}
	}
	thr = std::max(instance.minRho, instance.minDist);
	for (auto p: sources)
	{
		for (auto q: detectors)
		{
			if (p == q) return false;
			if (getDistance(instance.gdata.positions[p], instance.gdata.positions[q]) < thr) return false;
		}
	}
	return true;
}


void BSolution::print(bool optodes) const
{
	std::cout << "Solution: signal=" << signal << " [" << (signal / instance.sigNorm) << "] coverage="
	          << coverage << " [" << (coverage / instance.covNorm) << "] obj=" << objective << std::endl;
	if (optodes)
	{
		std::cout << "Sources:";
		printList(std::cout, sources.begin(), sources.end());
		std::cout << std::endl;
		std::cout << "Detectors:";
		printList(std::cout, detectors.begin(), detectors.end());
		std::cout << std::endl;
	}
}

/* Output format:
 * list of sources positions (Python format, 0-based index)
 * list of detectors positions (Python format, 0-based index)
 * signal
 * signal percentage
 * coverage
 * coverage percentage
 */
void BSolution::save(std::ostream& out) const
{
	printList(out, sources.begin(), sources.end(), 10, true);
	printList(out, detectors.begin(), detectors.end(), 10, true);
	out << signal << std::endl;
	out << (signal / instance.sigNorm) << std::endl;
	out << coverage << std::endl;
	out << (coverage / instance.covNorm) << std::endl;
}


Solution& Solution::operator=(const Solution& other)
{
	if (this != &other)
	{
		DOMINIQS_ASSERT(&instance == &(other.instance));
		sensors = other.sensors;
		coverage = other.coverage;
		signal = other.signal;
		objective = other.objective;
	}
	return *this;
}


void Solution::reset()
{
	sensors.clear();
	coverage = 0.0;
	signal = 0.0;
	objective = 0.0;
}


void Solution::computeSignal(std::vector<double>& cum) const
{
	std::fill(cum.begin(), cum.end(), 0.0);
	for (int s: sensors)
	{
		if (s == -1) continue;
		for (auto d: sensors)
		{
			if (d == -1) continue;
			if (d >= s) continue; //< do not count the same pair twice
			PMDFPtr pmdf = instance.getPMDF(s, d);
			if (!pmdf) continue;
			for (const PMDF::Item& item: pmdf->data) cum[item.voxel] += item.value;
		}
	}
}


void Solution::eval()
{
	coverage = 0.0;
	signal = 0.0;
	objective = 0.0;
	if (sensors.empty()) return;

	// compute total signal in each voxel
	std::vector<double> cum(instance.gdata.voxels.size(), 0.0);
	computeSignal(cum);

	// evaluate objective (that depends on the ROI)
	instance.evalObjectives(cum, coverage, signal, objective);
}


bool Solution::checkSignal(const std::vector<double>& cum) const
{
	std::vector<double> cumok(instance.gdata.voxels.size(), 0.0);
	computeSignal(cumok);
	for (auto v: instance.voxels)
	{
		if (fabs(cum[v]-cumok[v]) > 1e-9)
		{
			std::cout << "Check failed for voxel " << v << ": " << cum[v] << " != " << cumok[v] << std::endl;
			return false;
		}
	}
	return true;
}


bool Solution::checkCompatibilities() const
{
	double thr = instance.minDist;
	//< or should it be this one??? thr = std::max(instance.minRho, instance.minDist);
	for (auto p: sensors)
	{
		for (auto q: sensors)
		{
			if (p == q) continue;
			if (getDistance(instance.gdata.positions[p], instance.gdata.positions[q]) < thr) return false;
		}
	}
	return true;
}


void Solution::print(bool optodes) const
{
	std::cout << "Solution: signal=" << signal << " [" << (signal / instance.sigNorm) << "] coverage="
	          << coverage << " [" << (coverage / instance.covNorm) << "] obj=" << objective << std::endl;
	if (optodes)
	{
		std::cout << "Sensors:";
		printList(std::cout, sensors.begin(), sensors.end());
		std::cout << std::endl;
	}
}


/* Output format:
 * list of sensor positions (Python format, 0-based index)
 * signal
 * signal percentage
 * coverage
 * coverage percentage
 */
void Solution::save(std::ostream& out) const
{
	printList(out, sensors.begin(), sensors.end(), 10, true);
	out << signal << std::endl;
	out << (signal / instance.sigNorm) << std::endl;
	out << coverage << std::endl;
	out << (coverage / instance.covNorm) << std::endl;
}

