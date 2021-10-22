/**
 * @file nirs.h
 * @brief NIRS Data Structures
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2017 Domenico Salvagnin
 */

#ifndef NIRS_H
#define NIRS_H

#include <vector>
#include <memory>
#include <list>
#include <tuple>
#include <map>
#include <string>
#include <cmath>

class Instance;
class Solution;

const double COVERAGE_WEIGHT = 10.0;
const double COVERAGE_THR    =  0.1;


/* Point in 3D space (coordinates are in mm) */
class Point3D
{
public:
	Point3D(double _x, double _y, double _z) : data(_x, _y, _z) {}
	inline double x() const { return std::get<0>(data); }
	inline double y() const { return std::get<1>(data); }
	inline double z() const { return std::get<2>(data); }
	inline void set_x(double x) { std::get<0>(data) = x; }
	inline void set_y(double y) { std::get<1>(data) = y; }
	inline void set_z(double z) { std::get<2>(data) = z; }
protected:
	std::tuple<double, double, double> data;
};

std::ostream& operator<<(std::ostream& os, const Point3D& p);

/* Euclidian distance between two 3D points (in mm) */
inline double getDistance(const Point3D& p1, const Point3D& p2)
{
	double ret = 0.0;
	ret += (p1.x() - p2.x())*(p1.x() - p2.x());
	ret += (p1.y() - p2.y())*(p1.y() - p2.y());
	ret += (p1.z() - p2.z())*(p1.z() - p2.z());
	return std::sqrt(ret);
}


/* @brief A PMDF gives the signal measured in the voxels for a given
 * source/detector pair in positions (p,q).
 */
class PMDF
{
public:
	PMDF(int _p, int _q) : p(_p), q(_q) {}
	// data
	int p; //< position of source (index into NIRSData::positions)
	int q; //< position of detector (index into NIRSData::positions)
	struct Item
	{
		Item(int v = 0, double val = 0.0) : voxel(v), value(val) {}
		int voxel;
		double value; //< signal value at a given voxel
	};
	std::vector<Item> data;
	// work data
	double sum = 0.0; //< sum of signal values for the ROI Instance::voxels
};

typedef std::shared_ptr<PMDF> PMDFPtr;


/* @brief Global data about NIRS models.
 *
 * It defines the global set of positions/voxels given the current mesh.
 */
class NIRSData
{
public:
	// data
	std::vector<Point3D> positions; //< size P
	std::vector<Point3D> voxels; //< size V
	std::vector<PMDFPtr> allpmdfs; //<< size NIRSData::positions x NIRSData::positions
	std::vector<double> weights; //< size PxP
	// methods
	void read(const std::string& posfile,
				const std::string& voxfile,
				const std::string& pmdffileidx,
				const std::string& pmdffiledata,
				const std::string& wfile);
};


/** @brief: An instance of the NIRS problem
 *
 * A instance defines a subset of the voxels (a.k.a. as ROI) we are interested in.
 * Accordingly it also defines a subset of positions we can place optodes in.
 * It also defines the number of sources/detectors we want to use.
 */
class Instance
{
public:
	Instance(const NIRSData& data) : gdata(data), offset(gdata.positions.size()) {}
	~Instance();
	// data
	const NIRSData& gdata;
	size_t nSources = 0;
	size_t nDetectors = 0;
	std::vector<int> positions; //< subset of NIRSData::positions which are candidate for this instance
	std::vector<int> voxels; //< subset of NIRSData::voxels that define the ROI for this instance
	double minRho = 15; //< minimum distance between a source-detector pair (in mm)
	// double maxRho = 1000; //< maximum distance between a source-detector pair (in mm)
	double minDist = 10; //< minimum distance between any two optodes (in mm)
	double covNorm = 1.0; //< denominator to normalize coverage value
	double sigNorm = 1.0; //< denominator to normalize signal value
	double coverageThr = COVERAGE_THR;
	double coverageWeight = COVERAGE_WEIGHT;
	// Read from file
	void read(size_t ns, size_t nd, const std::string& filename);
	// Getters
	PMDFPtr getPMDF(int p, int q) const;
	double getDist(int p, int q) const;
	// Sparse representation of PMDF pairs
	std::vector<PMDFPtr> pmdfs;
	std::vector<size_t> start;
	// Evaluate objectives value given the signal on all voxels
	void evalObjectives(const std::vector<double>& cum, double& cov, double& sig, double& obj) const;
protected:
	void copyPMDFs(); //< helper
	size_t offset; //< cache offset to speed-up getPMDF/getDistance
	std::vector<PMDFPtr> pmdfpairs; //<< size NIRSData::positions x NIRSData::positions
	std::vector<double> distances; //<< size NIRSData::positions x NIRSData::positions
};


/** @brief: A bipartite solution to an instance of NIRS
 *
 * A solution is defined as the positions of the optodes (sources/detectors).
 * Both need to be subsets of Instance::positions, of cardinality Instance::nSources
 * and Instance::nDetectors, respectively.
 */
class BSolution
{
public:
	BSolution(const Instance& inst) : instance(inst) {}
	BSolution& operator=(const BSolution& other);
	const Instance& instance;
	std::vector<int> sources; //< subset of Instance::positions, of cardinality Instance::nSources
	std::vector<int> detectors; //< subset of Instance::positions, of cardinality Instance::nDetectors
	double coverage = 0.0; //< how many voxels of Instance::voxels are covered
	double signal = 0.0; //< cumulative signal in Instance::voxels
	double objective = 0.0; //< weighted objective value
	void reset();
	void computeSignal(std::vector<double>& cum) const;
	void eval();
	bool checkSignal(const std::vector<double>& cum) const;
	bool checkCompatibilities() const;
	void print(bool optodes=true) const;
	void save(std::ostream& out) const;
};


/** @brief: A non-bipartite solution to an instance of NIRS
 *
 * A solution is defined as the positions of the optodes.
 * Optodes act as both sources and detectors: their position need to be a subset of Instance::positions
 * of cardinality Instance::nSources = Instance::nDetectors.
 */
class Solution
{
public:
	Solution(const Instance& inst) : instance(inst) {}
	Solution& operator=(const Solution& other);
	const Instance& instance;
	std::vector<int> sensors; //< subset of Instance::positions, of cardinality Instance::nSources = Instance::nDetectors
	double coverage = 0.0; //< how many voxels of Instance::voxels are covered
	double signal = 0.0; //< cumulative signal in Instance::voxels
	double objective = 0.0; //< weighted objective value
	void reset();
	void computeSignal(std::vector<double>& cum) const;
	void eval();
	bool checkSignal(const std::vector<double>& cum) const;
	bool checkCompatibilities() const;
	void print(bool optodes=true) const;
	void save(std::ostream& out) const;
};

#endif /* end of include guard: NIRS_H */
