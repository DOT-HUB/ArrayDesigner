/**
 * @file nirsheur-bipartite.cpp
 * @brief Heuristics for bipartite NIRS
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2017 Domenico Salvagnin
 */

#include "nirs/nirsheur.h"
#include "nirs/trail.h"
#include "nirs/asserter.h"
#include "nirs/local_search.h"
#include <algorithm>
#include <iostream>
#include <random>

using namespace dominiqs;

typedef std::tuple<int,int,double> OptoPairObj;

static const double OBJEPS = 1e-6;

class UserData
{
public:
	int maxSignal = false;
};


static void filterPositions(const Instance& inst, int point, double thr, std::vector<int>& candidates)
{
	auto func = [&] (int p)
	{
		return (inst.getDist(p, point) < thr);
	};
	auto newend = std::remove_if(candidates.begin(), candidates.end(), func);
	candidates.erase(newend, candidates.end());
}


struct LocalDataWeightedObj
{
public:
	LocalDataWeightedObj(const Instance& inst)
		: instance(inst), cum(inst.gdata.voxels.size(), 0.0), tr(cum),
		sourcePositions(inst.positions), detectorPositions(inst.positions) {}
	// data
	const Instance& instance;
	std::vector<double> cum;
	trail<double> tr;
	std::vector<int> sourcePositions;
	std::vector<int> detectorPositions;
	// API
	// eval* methods use the trail
	// append/set/remove method do not
	void initWithBSolution(const BSolution& sol)
	{
		sol.computeSignal(cum);
	}
	double getCurrentObjective() const
	{
		double cov = 0.0;
		double sig = 0.0;
		double obj = 0.0;
		instance.evalObjectives(cum, cov, sig, obj);
		return obj;
	}
	// NOTE: does not use trail!
	void appendSource(BSolution& sol, int p)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		sol.sources.push_back(p);
		addSourceContribution(sol, p, false);
		filterPositions(instance, p, instance.minDist, sourcePositions);
		filterPositions(instance, p, std::max(instance.minDist, instance.minRho), detectorPositions);
	}
	// NOTE: does not use trail!
	void appendDetector(BSolution& sol, int q)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		sol.detectors.push_back(q);
		addDetectorContribution(sol, q, false);
		filterPositions(instance, q, std::max(instance.minDist, instance.minRho), sourcePositions);
		filterPositions(instance, q, instance.minDist, detectorPositions);
	}
	// evaluate objective if we added a given PMDF
	double evalPMDF(PMDFPtr pmdf)
	{
		double cov = 0.0;
		double sig = 0.0;
		double obj = 0.0;
		// add current pmdf to cum (using trail!)
		for (const PMDF::Item& item: pmdf->data) tr.set(item.voxel, cum[item.voxel] + item.value);
		// evaluate objective
		instance.evalObjectives(cum, cov, sig, obj);
		// undo
		tr.restore();
		return obj;
	}
	// evaluate objective if we added a source in position p to sol (but without actually adding it)
	double evalAddSource(BSolution& sol, int p)
	{
		double cov = 0.0;
		double sig = 0.0;
		double obj = 0.0;
		addSourceContribution(sol, p, true);
		instance.evalObjectives(cum, cov, sig, obj);
		tr.restore();
		return obj;
	}
	// evaluate objective if we added a detector in position q to sol (but without actually adding it)
	double evalAddDetector(BSolution& sol, int q)
	{
		double cov = 0.0;
		double sig = 0.0;
		double obj = 0.0;
		addDetectorContribution(sol, q, true);
		instance.evalObjectives(cum, cov, sig, obj);
		tr.restore();
		return obj;
	}
	// Remove the source in position of index i form sol
	void removeSource(BSolution& sol, int i)
	{
		int p = sol.sources[i];
		// temporarily mark p as invalid
		sol.sources[i] = -1;
		// remove contribution of source p (without trail!)
		removeSourceContribution(sol, p, false);
	}
	// Set source in position of index i to p
	void setSource(BSolution& sol, int i, int p)
	{
		DOMINIQS_ASSERT( sol.sources[i] == -1 );
		sol.sources[i] = p;
		addSourceContribution(sol, p, false);
	}
	// Remove detector in position of index j
	void removeDetector(BSolution& sol, int j)
	{
		int q = sol.detectors[j];
		// temporarily mark q as invalid
		sol.detectors[j] = -1;
		// remove contribution of detector q (without trail!)
		removeDetectorContribution(sol, q, false);
	}
	// Set detector in position of index j to q
	void setDetector(BSolution& sol, int j, int q)
	{
		DOMINIQS_ASSERT( sol.detectors[j] == -1 );
		sol.detectors[j] = q;
		addDetectorContribution(sol, q, false);
	}
	// Eval objective if we removed the pair of optodes in positions of indices (i,j)
	double evalRemovePair(BSolution& sol, int i, int j)
	{
		int p = sol.sources[i];
		int q = sol.detectors[j];
		sol.sources[i] = -1;
		removeSourceContribution(sol, p, true);
		sol.detectors[j] = -1;
		removeDetectorContribution(sol, q, true);
		double obj = getCurrentObjective();
		tr.restore();
		sol.sources[i] = p;
		sol.detectors[j] = q;
		return obj;
	}
	// Remove (source,detector) in positions of index (i,j)
	void removePair(BSolution& sol, int i, int j)
	{
		int p = sol.sources[i];
		int q = sol.detectors[j];
		sol.sources[i] = -1;
		removeSourceContribution(sol, p, false);
		sol.detectors[j] = -1;
		removeDetectorContribution(sol, q, false);
	}
	// Eval objective if we added a pair of optodes in positions (p,q) in indices (i,j)
	double evalSetPair(BSolution& sol, int i, int j, int p, int q)
	{
		DOMINIQS_ASSERT( sol.sources[i] == -1 );
		DOMINIQS_ASSERT( sol.detectors[j] == -1 );
		sol.sources[i] = p;
		addSourceContribution(sol, p, true);
		sol.detectors[j] = q;
		addDetectorContribution(sol, q, true);
		DOMINIQS_ASSERT(sol.checkSignal(cum));
		double obj = getCurrentObjective();
		sol.sources[i] = -1;
		sol.detectors[j] = -1;
		tr.restore();
		return obj;
	}
	// Set pair of optodes in positions (p,q) in indices (i,j)
	void setPair(BSolution& sol, int i, int j, int p, int q)
	{
		DOMINIQS_ASSERT( sol.sources[i] == -1 );
		DOMINIQS_ASSERT( sol.detectors[j] == -1 );
		sol.sources[i] = p;
		addSourceContribution(sol, p, false);
		sol.detectors[j] = q;
		addDetectorContribution(sol, q, false);
	}
protected:
	// helpers
	void addSourceContribution(const BSolution& sol, int p, bool useTrail)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		for (int q: sol.detectors)
		{
			if (q == -1) continue; // skip invalid detectors
			PMDFPtr pmdf = instance.getPMDF(p, q);
			if (!pmdf) continue;
			if (useTrail) for (const PMDF::Item& item: pmdf->data) tr.set(item.voxel, cum[item.voxel] + item.value);
			else for (const PMDF::Item& item: pmdf->data) cum[item.voxel] += item.value;
		}
	}
	void addDetectorContribution(const BSolution& sol, int q, bool useTrail)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		for (int p: sol.sources)
		{
			if (p == -1) continue; // skip invalid sources
			PMDFPtr pmdf = instance.getPMDF(p, q);
			if (!pmdf) continue;
			if (useTrail) for (const PMDF::Item& item: pmdf->data) tr.set(item.voxel, cum[item.voxel] + item.value);
			else for (const PMDF::Item& item: pmdf->data) cum[item.voxel] += item.value;
		}
	}
	void removeSourceContribution(const BSolution& sol, int p, bool useTrail)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		for (int q: sol.detectors)
		{
			if (q == -1) continue; // skip invalid detectors
			PMDFPtr pmdf = instance.getPMDF(p, q);
			if (!pmdf) continue;
			if (useTrail) for (const PMDF::Item& item: pmdf->data) tr.set(item.voxel, cum[item.voxel] - item.value);
			else for (const PMDF::Item& item: pmdf->data) cum[item.voxel] -= item.value;
		}
	}
	void removeDetectorContribution(const BSolution& sol, int q, bool useTrail)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		for (int p: sol.sources)
		{
			if (p == -1) continue; // skip invalid sources
			PMDFPtr pmdf = instance.getPMDF(p, q);
			if (!pmdf) continue;
			if (useTrail) for (const PMDF::Item& item: pmdf->data) tr.set(item.voxel, cum[item.voxel] - item.value);
			else for (const PMDF::Item& item: pmdf->data) cum[item.voxel] -= item.value;
		}
	}
};


struct LocalDataMaxSignal
{
public:
	LocalDataMaxSignal(const Instance& inst)
		: instance(inst), sourceDelta(inst.gdata.positions.size(), 0), detectorDelta(inst.gdata.positions.size(), 0),
		sourcePositions(inst.positions), detectorPositions(inst.positions) {}
	// data
	const Instance& instance;
	std::vector<double> sourceDelta;
	std::vector<double> detectorDelta;
	std::vector<int> sourcePositions;
	std::vector<int> detectorPositions;
	double currentObjective = 0.0;
	// API
	void initWithBSolution(const BSolution& sol)
	{
		// set up deltas
		std::fill(sourceDelta.begin(), sourceDelta.end(), 0);
		std::fill(detectorDelta.begin(), detectorDelta.end(), 0);
		// first delta for sources
		for (int p: instance.positions)
		{
			double term = 0.0;
			for (int q: sol.detectors)
			{
				if (q == -1) continue; // skip invalid detectors
				PMDFPtr pmdf = instance.getPMDF(p, q);
				if (!pmdf) continue;
				term += pmdf->sum;
			}
			sourceDelta[p] = term;
		}
		// then delta for detectors
		for (int q: instance.positions)
		{
			double term = 0.0;
			for (int p: sol.sources)
			{
				if (p == -1) continue; // skip invalid sources
				PMDFPtr pmdf = instance.getPMDF(p, q);
				if (!pmdf) continue;
				term += pmdf->sum;
			}
			detectorDelta[q] = term;
		}
		// compute current objective
		currentObjective = 0.0;
		for (int p: sol.sources)
		{
			if (p == -1) continue; // skip invalid sources
			for (int q: sol.detectors)
			{
				if (q == -1) continue; // skip invalid detectors
				PMDFPtr pmdf = instance.getPMDF(p, q);
				if (!pmdf) continue;
				currentObjective += pmdf->sum;
			}
		}
	}
	double getCurrentObjective() const
	{
		return currentObjective;
	}
	void appendSource(BSolution& sol, int p)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		sol.sources.push_back(p);
		addSourceContribution(sol, p);
		filterPositions(instance, p, instance.minDist, sourcePositions);
		filterPositions(instance, p, std::max(instance.minDist, instance.minRho), detectorPositions);
	}
	void appendDetector(BSolution& sol, int q)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		sol.detectors.push_back(q);
		addDetectorContribution(sol, q);
		filterPositions(instance, q, std::max(instance.minDist, instance.minRho), sourcePositions);
		filterPositions(instance, q, instance.minDist, detectorPositions);
	}
	void clearSources(BSolution& sol)
	{
		sol.sources.clear();
		currentObjective = 0.0;
		std::fill(detectorDelta.begin(), detectorDelta.end(), 0.0);
	}
	void clearDetectors(BSolution& sol)
	{
		sol.detectors.clear();
		currentObjective = 0.0;
		std::fill(sourceDelta.begin(), sourceDelta.end(), 0.0);
	}
	// evaluate objective if we added a given PMDF
	double evalPMDF(PMDFPtr pmdf)
	{
		return pmdf->sum;
	}
	// evaluate objective if we added a source in position p to sol (but without actually adding it)
	double evalAddSource(BSolution& sol, int p)
	{
		addSourceContribution(sol, p);
		double obj = currentObjective;
		removeSourceContribution(sol, p);
		return obj;
	}
	// evaluate objective if we added a detector in position q to sol (but without actually adding it)
	double evalAddDetector(BSolution& sol, int q)
	{
		addDetectorContribution(sol, q);
		double obj = currentObjective;
		removeDetectorContribution(sol, q);
		return obj;
	}
	// Remove the source in position of index i from sol
	void removeSource(BSolution& sol, int i)
	{
		int p = sol.sources[i];
		// temporarily mark p as invalid
		sol.sources[i] = -1;
		// remove contribution of source p
		removeSourceContribution(sol, p);
	}
	// Set source in position of index i to p
	void setSource(BSolution& sol, int i, int p)
	{
		DOMINIQS_ASSERT( sol.sources[i] == -1 );
		sol.sources[i] = p;
		addSourceContribution(sol, p);
	}
	// Remove detector in position of index j
	void removeDetector(BSolution& sol, int j)
	{
		int q = sol.detectors[j];
		// temporarily mark q as invalid
		sol.detectors[j] = -1;
		// remove contribution of detector q
		removeDetectorContribution(sol, q);
	}
	// Set detector in position of index j to q
	void setDetector(BSolution& sol, int j, int q)
	{
		DOMINIQS_ASSERT( sol.detectors[j] == -1 );
		sol.detectors[j] = q;
		addDetectorContribution(sol, q);
	}
	// Eval objective if we removed the pair of optodes in positions of indices (i,j)
	double evalRemovePair(BSolution& sol, int i, int j)
	{
		int p = sol.sources[i];
		int q = sol.detectors[j];
		removePair(sol, i, j);
		double obj = currentObjective;
		setPair(sol, i, j, p, q);
		return obj;
	}
	// Remove (source,detector) in positions of index (i,j)
	void removePair(BSolution& sol, int i, int j)
	{
		int p = sol.sources[i];
		int q = sol.detectors[j];
		sol.sources[i] = -1;
		removeSourceContribution(sol, p);
		sol.detectors[j] = -1;
		removeDetectorContribution(sol, q);
	}
	// Eval objective if we added a pair of optodes in positions (p,q) in indices (i,j)
	double evalSetPair(BSolution& sol, int i, int j, int p, int q)
	{
		DOMINIQS_ASSERT( sol.sources[i] == -1 );
		DOMINIQS_ASSERT( sol.detectors[j] == -1 );
		setPair(sol, i, j, p, q);
		double obj = currentObjective;
		removePair(sol, i, j);
		return obj;
	}
	// Set pair of optodes in positions (p,q) in indices (i,j)
	void setPair(BSolution& sol, int i, int j, int p, int q)
	{
		DOMINIQS_ASSERT( sol.sources[i] == -1 );
		DOMINIQS_ASSERT( sol.detectors[j] == -1 );
		sol.sources[i] = p;
		addSourceContribution(sol, p);
		sol.detectors[j] = q;
		addDetectorContribution(sol, q);
	}
protected:
	// helpers
	void addSourceContribution(const BSolution& sol, int p)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		// update objective
		for (int q: sol.detectors)
		{
			if (q == -1) continue; // skip invalid detectors
			PMDFPtr pmdf = instance.getPMDF(p, q);
			if (!pmdf) continue;
			currentObjective += pmdf->sum;
		}
		// update deltas for detectors
		for (int k = instance.start[p]; k < instance.start[p+1]; k++)
		{
			PMDFPtr pmdf = instance.pmdfs[k];
			DOMINIQS_ASSERT( pmdf );
			detectorDelta[pmdf->q] += pmdf->sum;
		}
	}
	void addDetectorContribution(const BSolution& sol, int q)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		// update objective
		for (int p: sol.sources)
		{
			if (p == -1) continue; // skip invalid sources
			PMDFPtr pmdf = instance.getPMDF(p, q);
			if (!pmdf) continue;
			currentObjective += pmdf->sum;
		}
		// update deltas for sources
		for (int k = instance.start[q]; k < instance.start[q+1]; k++)
		{
			PMDFPtr pmdf = instance.pmdfs[k];
			DOMINIQS_ASSERT( pmdf );
			sourceDelta[pmdf->q] += pmdf->sum; //< yes: pmdf->q and not pmdf->p!
		}
	}
	void removeSourceContribution(const BSolution& sol, int p)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		// update objective
		for (int q: sol.detectors)
		{
			if (q == -1) continue; // skip invalid detectors
			PMDFPtr pmdf = instance.getPMDF(p, q);
			if (!pmdf) continue;
			currentObjective -= pmdf->sum;
		}
		// update deltas for detectors
		for (int k = instance.start[p]; k < instance.start[p+1]; k++)
		{
			PMDFPtr pmdf = instance.pmdfs[k];
			DOMINIQS_ASSERT( pmdf );
			detectorDelta[pmdf->q] -= pmdf->sum;
		}
	}
	void removeDetectorContribution(const BSolution& sol, int q)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		// update objective
		for (int p: sol.sources)
		{
			if (p == -1) continue; // skip invalid sources
			PMDFPtr pmdf = instance.getPMDF(p, q);
			if (!pmdf) continue;
			currentObjective -= pmdf->sum;
		}
		// update deltas for sources
		for (int k = instance.start[q]; k < instance.start[q+1]; k++)
		{
			PMDFPtr pmdf = instance.pmdfs[k];
			DOMINIQS_ASSERT( pmdf );
			sourceDelta[pmdf->q] -= pmdf->sum; //< yes: pmdf->q and not pmdf->p!
		}
	}
};


static bool compatible(const Instance& inst, int p, int q)
{
	if (p == q) return false;
	if (getDistance(inst.gdata.positions[p], inst.gdata.positions[q]) < std::max(inst.minDist, inst.minRho)) return false;
	return true;
}


typedef std::pair<int,double> PosObj;
static int pickCandidate(std::vector<PosObj>& candlist, size_t maxCand, std::mt19937& rnd)
{
	if (candlist.empty()) return -1;
	std::sort(candlist.begin(), candlist.end(), [](const PosObj& a, const PosObj& b) { return a.second > b.second; });
	if (candlist.size() > maxCand) candlist.resize(maxCand);
	std::shuffle(candlist.begin(), candlist.end(), rnd);
	return candlist[0].first;
}


static void graspIteration(const Instance& inst, BSolution& sol, size_t maxCand, size_t mult, std::mt19937& rnd)
{
	// setup
	sol.reset();
	LocalDataWeightedObj d(inst);

	// pick initial pair of (source,detectors)
	std::vector<OptoPairObj> candidates;
	for (int p: inst.positions)
	{
		for (int q: inst.positions)
		{
			if ((p >= q) || !compatible(inst,p,q)) continue;
			PMDFPtr pmdf = inst.getPMDF(p, q);
			if (!pmdf) continue;
			candidates.emplace_back(p, q, d.evalPMDF(pmdf));
		}
	}
	DOMINIQS_ASSERT(!candidates.empty());
	std::sort(candidates.begin(), candidates.end(), [](const OptoPairObj& a, const OptoPairObj& b) { return std::get<2>(a) > std::get<2>(b); });
	if (candidates.size() > mult * maxCand) candidates.resize(maxCand);
	std::shuffle(candidates.begin(), candidates.end(), rnd);
	OptoPairObj best = candidates[0];
	int p = std::get<0>(best);
	int q = std::get<1>(best);
	d.appendSource(sol, p);
	d.appendDetector(sol, q);
	DOMINIQS_ASSERT(sol.checkSignal(d.cum));
	DOMINIQS_ASSERT(sol.checkCompatibilities());
	// now we need to add one source and one detector at the time
	std::vector<PosObj> candlist;
	while(true)
	{
		// add a source
		if (sol.sources.size() < inst.nSources)
		{
			candlist.clear();
			for (int p: d.sourcePositions) candlist.emplace_back(p, d.evalAddSource(sol, p));
			int bestp = pickCandidate(candlist, maxCand, rnd);
			if (bestp == -1) break;
			// bestp is now the new source to add
			d.appendSource(sol, bestp);
			DOMINIQS_ASSERT(sol.checkSignal(d.cum));
			DOMINIQS_ASSERT(sol.checkCompatibilities());
		}
		// add a detector
		if (sol.detectors.size() < inst.nDetectors)
		{
			candlist.clear();
			for (int q: d.detectorPositions)  candlist.emplace_back(q, d.evalAddDetector(sol, q));
			int bestq = pickCandidate(candlist, maxCand, rnd);
			if (bestq == -1) break;
			// bestq is now the new detector to add
			d.appendDetector(sol, bestq);
			DOMINIQS_ASSERT(sol.checkSignal(d.cum));
			DOMINIQS_ASSERT(sol.checkCompatibilities());
		}
		// termination condition
		if ((sol.sources.size() >= inst.nSources) && (sol.detectors.size() >= inst.nDetectors)) break;
	}
	sol.eval();
}


template <typename UserData>
inline double getSolutionObjValue(const BSolution& sol, const UserData& ud)
{
	return ud.maxSignal ? -sol.signal : -sol.objective;
}


void grasp(const Instance& inst, BSolution& sol, TerminationCriteria& toStop, bool optTotSignal, size_t maxCand, int seed)
{
	std::mt19937 rnd(seed);
	int iter = 0;
	int noImproveCount = 0;
	UserData ud;
	ud.maxSignal = optTotSignal;

	while (true)
	{
		if (toStop(iter, noImproveCount)) break;
		BSolution currsol(inst);
		graspIteration(inst, currsol, iter ? maxCand : 1, iter ? 5 : 1, rnd);
		localSearch(inst, currsol, optTotSignal);
		if (getSolutionObjValue(currsol, ud) < getSolutionObjValue(sol, ud) - OBJEPS)
		{
			sol = currsol;
			noImproveCount = 0;
		}
		else noImproveCount++;
		std::cout << "Iteration: " << iter << " -> ";
		currsol.print(false);
		iter++;
	}
}


static void computeCandidateSources(const Instance& inst, const BSolution& sol, std::vector<int>& sourcePositions)
{
	sourcePositions = inst.positions;
	// cannot be too close to other sources
	for (int p: sol.sources)
	{
		if (p == -1) continue;
		filterPositions(inst, p, inst.minDist, sourcePositions);
	}
	// cannot be too close to other detectors
	for (int q: sol.detectors)
	{
		if (q == -1) continue;
		filterPositions(inst, q, std::max(inst.minDist, inst.minRho), sourcePositions);
	}
}


static void computeCandidateDetectors(const Instance& inst, const BSolution& sol, std::vector<int>& detectorPositions)
{
	detectorPositions = inst.positions;
	// cannot be too close to other detectors
	for (int q: sol.detectors)
	{
		if (q == -1) continue;
		filterPositions(inst, q, inst.minDist, detectorPositions);
	}
	// cannot be too close to other sources
	for (int p: sol.sources)
	{
		if (p == -1) continue;
		filterPositions(inst, p, std::max(inst.minDist, inst.minRho), detectorPositions);
	}
}


template<typename Observer, typename LocalData>
bool oneOpt(BSolution& sol, Observer& rec)
{
	const Instance& inst = sol.instance;
	LocalData d(inst);
	d.initWithBSolution(sol);
	double oldobj = d.getCurrentObjective();
	bool improved = false;
	// look for improving p
	size_t nsources = sol.sources.size();
	for (size_t i = 0; i < nsources; i++)
	{
		int p = sol.sources[i];
		d.removeSource(sol, i);
		// compute set of candidate positions
		std::vector<int> candidates;
		computeCandidateSources(inst, sol, candidates);
		// try each of the candidates
		int imprp = p;
		for (int newp: candidates)
		{
			if (newp == p) continue;
			// eval contribution of source newp
			double obj = d.evalAddSource(sol, newp);
			if (obj > oldobj + OBJEPS)
			{
				// we have found an improving move
				imprp = newp;
				improved = true;
				break;
			}
		}
		// apply move (could be the no move at all if imprp == p)
		d.setSource(sol, i, imprp);
		if (improved)
		{
			sol.eval();
			break;
		}
	}
	DOMINIQS_ASSERT(sol.checkCompatibilities());

	if (improved) return improved;

	size_t ndetectors = sol.detectors.size();
	for (size_t j = 0; j < ndetectors; j++)
	{
		int q = sol.detectors[j];
		d.removeDetector(sol, j);
		// compute set of candidate positions
		std::vector<int> candidates(inst.positions);
		computeCandidateDetectors(inst, sol, candidates);
		// try each of the candidates
		int imprq = q;
		for (int newq: candidates)
		{
			if (newq == q) continue;
			// eval contribution of detector newq
			double obj = d.evalAddDetector(sol, newq);
			if (obj > oldobj + OBJEPS)
			{
				// we have found an improving move
				imprq = newq;
				improved = true;
				break;
			}
		}
		// apply move (could be the no move at all if imprq == q)
		d.setDetector(sol, j, imprq);
		if (improved)
		{
			sol.eval();
			break;
		}
	}
	DOMINIQS_ASSERT(sol.checkCompatibilities());

	return improved;
}


template<typename Observer, typename LocalData>
bool twoOpt(BSolution& sol, Observer& rec)
{
	const Instance& inst = sol.instance;
	LocalData d(inst);
	d.initWithBSolution(sol);
	double oldobj = d.getCurrentObjective();
	bool improved = false;
	std::vector<OptoPairObj> toremove;
	// evaluate the change in objective for removing each pair
	size_t nsources = sol.sources.size();
	size_t ndetectors = sol.detectors.size();
	for (size_t i = 0; i < nsources; i++)
	{
		for (size_t j = 0; j < ndetectors; j++)
		{
			toremove.emplace_back(i,j,d.evalRemovePair(sol,i,j)-oldobj);
		}
	}
	// sort candidates by decreasing change in objective (better solutions first)
	std::sort(toremove.begin(), toremove.end(), [](const OptoPairObj& x, const OptoPairObj& y)
		{
			return std::get<2>(x) > std::get<2>(y);
		}
	);

	// try removing each pair (p,q)
	for (const OptoPairObj& x: toremove)
	{
		int i = std::get<0>(x);
		int j = std::get<1>(x);
		int p = sol.sources[i];
		int q = sol.detectors[j];
		d.removePair(sol, i, j);
		// now we have removed the pair (p,q)
		// we can start looking for a replacement

		// compute set of candidate positions for sources and detectors
		std::vector<int> candsources(inst.positions);
		computeCandidateSources(inst, sol, candsources);
		std::vector<int> canddetectors(inst.positions);
		computeCandidateDetectors(inst, sol, canddetectors);

		int imprp = p;
		int imprq = q;
		for (int newp: candsources)
		{
			if (newp == p) continue;
			for (int newq: canddetectors)
			{
				if (newq == q) continue;
				if (!compatible(inst,newp,newq)) continue;
				// eval contribution of (newp,newq) using trail
				double obj = d.evalSetPair(sol, i, j, newp, newq);
				if (obj > oldobj + OBJEPS)
				{
					// we have found an improving move
					improved = true;
					imprp = newp;
					imprq = newq;
					break;
				}
			}
			if (improved) break;
		}

		// apply move (again could be no move at all)
		d.setPair(sol, i, j, imprp, imprq);
		if (improved)
		{
			sol.eval();
			break;
		}
	}
	DOMINIQS_ASSERT(sol.checkCompatibilities());
	return improved;
}


static void filterCandidates(const Instance& inst, int point, double thr, std::vector<PosObj>& candidates)
{
	auto func = [&] (const PosObj& item)
	{
		return (inst.getDist(item.first, point) < thr);
	};
	auto newend = std::remove_if(candidates.begin(), candidates.end(), func);
	candidates.erase(newend, candidates.end());
}


void floatOpt(const Instance& instance, std::vector<int>& candlist, const std::vector<double>& allobj, size_t target)
{
	std::vector<PosObj> candidates;
	for (int k: candlist) candidates.emplace_back(k, allobj[k]);
	std::sort(candidates.begin(), candidates.end(), [](const PosObj& a, const PosObj& b) { return a.second > b.second; });
	candlist.clear();
	while (candlist.size() < target && !candidates.empty())
	{
		// pick first
		int newitem = candidates[0].first;
		candlist.push_back(newitem);
		// filter candidates
		filterCandidates(instance, newitem, instance.minDist, candidates);
	}
}


template<typename Observer>
bool twoOptFloat(BSolution& sol, Observer& rec)
{
	const Instance& inst = sol.instance;
	LocalDataMaxSignal d(inst);
	d.initWithBSolution(sol);
	double oldobj = d.getCurrentObjective();
	DOMINIQS_ASSERT(fabs(oldobj - sol.signal) < 1e-9);
	bool improved = false;

	std::vector<int> bestdetectors = sol.detectors;
	d.clearDetectors(sol);
	sol.eval();
	DOMINIQS_ASSERT(fabs(d.getCurrentObjective() - sol.signal) < 1e-9);

	size_t nsources = sol.sources.size();
	// look for improving p
	for (size_t i = 0; i < nsources; i++)
	{
		int p = sol.sources[i];
		d.removeSource(sol, i);
		// compute set of candidate positions
		std::vector<int> candsources;
		computeCandidateSources(inst, sol, candsources);
		// try each of the candidates
		int imprp = p;
		for (int newp: candsources)
		{
			// temporarily set the new source
			d.setSource(sol, i, newp);
			// find set of candidates detectors
			std::vector<int> canddetectors;
			computeCandidateDetectors(inst, sol, canddetectors);
			// compute best subset
			floatOpt(inst, canddetectors, d.detectorDelta, inst.nDetectors);
			// add best subset to solution to evaluate objective
			for (int q: canddetectors) d.appendDetector(sol, q);
			double obj = d.getCurrentObjective();
			// undo for next iteration
			d.removeSource(sol, i);
			d.clearDetectors(sol);
			if (obj > oldobj + OBJEPS)
			{
				// we have found an improving move
				imprp = newp;
				bestdetectors = canddetectors;
				improved = true;
				break;
			}
		}
		// apply move
		d.setSource(sol, i, imprp);
		if (improved)
		{
			break;
		}
	}
	// we need to add the detectors back
	for (int q: bestdetectors) d.appendDetector(sol, q);
	sol.eval();
	DOMINIQS_ASSERT(sol.checkCompatibilities());

	if (improved) return improved;

	std::vector<int> bestsources = sol.sources;
	d.clearSources(sol);
	sol.eval();
	DOMINIQS_ASSERT(fabs(d.getCurrentObjective() - sol.signal) < 1e-9);

	size_t ndetectors = sol.detectors.size();
	for (size_t j = 0; j < ndetectors; j++)
	{
		int q = sol.detectors[j];
		d.removeDetector(sol, j);
		// compute set of candidate positions
		std::vector<int> canddetectors(inst.positions);
		computeCandidateDetectors(inst, sol, canddetectors);
		// try each of the candidates
		int imprq = q;
		for (int newq: canddetectors)
		{
			// temporarily set the new source
			d.setDetector(sol, j, newq);
			// find set of candidates sources
			std::vector<int> candsources;
			computeCandidateSources(inst, sol, candsources);
			// compute best subset
			floatOpt(inst, candsources, d.sourceDelta, inst.nSources);
			// add best subset to solution to evaluate objective
			for (int p: candsources) d.appendSource(sol, p);
			double obj = d.getCurrentObjective();
			// undo for next iteration
			d.removeDetector(sol, j);
			d.clearSources(sol);
			if (obj > oldobj + OBJEPS)
			{
				// we have found an improving move
				imprq = newq;
				bestsources = candsources;
				improved = true;
				break;
			}
		}
		// apply move
		d.setDetector(sol, j, imprq);
		if (improved)
		{
			break;
		}
	}
	// we need to add the sources back
	for (int p: bestsources) d.appendSource(sol, p);
	sol.eval();
	DOMINIQS_ASSERT(sol.checkCompatibilities());

	return improved;
}


template<typename Observer, typename UserData=EmptyUserData>
bool exploreNeighborhood(BSolution& sol, Observer& rec, const UserData& ud)
{
	bool improved = false;
	if (ud.maxSignal)
	{
		improved = oneOpt<Observer,LocalDataMaxSignal>(sol, rec);
		if (!improved) improved = twoOptFloat(sol, rec);
		//if (!improved) improved = twoOpt<Observer,LocalDataMaxSignal>(sol, rec);
	}
	else
	{
		improved = oneOpt<Observer,LocalDataWeightedObj>(sol, rec);
		if (!improved) improved = twoOpt<Observer,LocalDataWeightedObj>(sol, rec);
	}
	return improved;
}


void localSearch(const Instance& inst, BSolution& sol, bool optTotSignal)
{
	IObserver<BSolution> rec; //< does nothing
	UserData ud;
	ud.maxSignal = optTotSignal;
	localSearch(sol, rec, ud);
}


template <typename RandomEngine, typename UserData>
void randomGenerateSolution(BSolution& sol, RandomEngine& engine, const UserData& ud)
{
	const Instance& inst = sol.instance;
	LocalDataWeightedObj d(inst);
	sol.reset();
	while(true)
	{
		// add a source
		if (sol.sources.size() < inst.nSources)
		{
			if (d.sourcePositions.empty()) break;
			std::uniform_int_distribution<> rndint(0, d.sourcePositions.size()-1);
			int elem = rndint(engine);
			int p = d.sourcePositions[elem];
			d.appendSource(sol, p);
			DOMINIQS_ASSERT(sol.checkCompatibilities());
		}
		// add a detector
		if (sol.detectors.size() < inst.nDetectors)
		{
			if (d.detectorPositions.empty()) break;
			std::uniform_int_distribution<> rndint(0, d.detectorPositions.size()-1);
			int elem = rndint(engine);
			int q = d.detectorPositions[elem];
			d.appendDetector(sol, q);
			DOMINIQS_ASSERT(sol.checkCompatibilities());
		}
		// termination condition
		if ((sol.sources.size() >= inst.nSources) && (sol.detectors.size() >= inst.nDetectors)) break;
	}
	sol.eval();
}


void rrls(const Instance& inst, BSolution& sol, TerminationCriteria& toStop, bool optTotSignal, int seed)
{
	IObserver<BSolution> rec; //< does nothing
	std::mt19937 rnd(seed);
	UserData ud;
	ud.maxSignal = optTotSignal;
	randomRestartLocalSearch(sol, rnd, rec, toStop, ud);
}


template<typename RandomEngine, typename UserData>
void perturbeSolution(BSolution& sol, RandomEngine& engine, int K, const UserData& ud)
{
	const Instance& inst = sol.instance;
	std::bernoulli_distribution randbool(0.5);
	for (int k = 0; k < K; k++)
	{
		bool perturbeSources = randbool(engine);
		if (perturbeSources)
		{
			// pick element to change
			DOMINIQS_ASSERT(sol.sources.size());
			std::uniform_int_distribution<> rndint(0, sol.sources.size()-1);
			int i = rndint(engine);
			DOMINIQS_ASSERT(i >= 0 && i < sol.sources.size());
			int p = sol.sources[i];
			std::vector<int> candidates(inst.positions);
			// cannot be too close to remaining sources
			for (int ap: sol.sources)
			{
				if (ap != p)
				{
					filterPositions(inst, ap, inst.minDist, candidates);
				}
			}
			// cannot be too close to detectors
			for (int q: sol.detectors)
			{
				filterPositions(inst, q, std::max(inst.minDist, inst.minRho), candidates);
			}
			// pick new value at random
			std::uniform_int_distribution<> rndint2(0, candidates.size()-1);
			int newp = candidates[rndint2(engine)];
			sol.sources[i] = newp;
			DOMINIQS_ASSERT(sol.checkCompatibilities());
		}
		else
		{
			// pick element to change
			DOMINIQS_ASSERT(sol.detectors.size());
			std::uniform_int_distribution<> rndint(0, sol.detectors.size()-1);
			int i = rndint(engine);
			DOMINIQS_ASSERT(i >= 0 && i < sol.detectors.size());
			int q = sol.detectors[i];
			std::vector<int> candidates(inst.positions);
			// cannot be too close to remaining detectors
			for (int aq: sol.detectors)
			{
				if (aq != q)
				{
					filterPositions(inst, aq, inst.minDist, candidates);
				}
			}
			// cannot be too close to sources
			for (int p: sol.sources)
			{
				filterPositions(inst, p, std::max(inst.minDist, inst.minRho), candidates);
			}
			// pick new value at random
			std::uniform_int_distribution<> rndint2(0, candidates.size()-1);
			int newq = candidates[rndint2(engine)];
			sol.detectors[i] = newq;
			DOMINIQS_ASSERT(sol.checkCompatibilities());
		}
	}
	sol.eval();
}


void ils(const Instance& inst, BSolution& sol, TerminationCriteria& toStop, bool optTotSignal, int seed)
{
	IObserver<BSolution> rec; //< does nothing
	std::mt19937 rnd(seed);
	// HillClimbingCriterion<BSolution> accept;
	RandomWalkCriterion<BSolution> accept;
	UserData ud;
	ud.maxSignal = optTotSignal;
	iteratedLocalSearch(sol, accept, rnd, rec, toStop,
								std::max(3, (int)inst.nSources),
								std::max(5, (int)inst.nSources),
								ud);
}


void multiILS(const Instance& inst, BSolution& sol, TerminationCriteria& toStop, bool optTotSignal, int maxTries, int seed)
{
	IObserver<BSolution> rec; //< does nothing
	std::mt19937 rnd(seed);
	BSolution currSol(sol.instance);
	UserData ud;
	ud.maxSignal = optTotSignal;
	for (int t = 0; t < maxTries; t++)
	{
		RandomWalkCriterion<BSolution> accept;
		iteratedLocalSearch(currSol, accept, rnd, rec, toStop,
									std::max(3, (int)inst.nSources),
									std::max(5, (int)inst.nSources),
									ud);
		if (currSol.objective > sol.objective + OBJEPS)
		{
			sol = currSol;
		}
	}
}
