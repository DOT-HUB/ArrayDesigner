/**
 * @file nirsheur.cpp
 * @brief Heuristics for NIRS
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2018 Domenico Salvagnin
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
		sensorPositions(inst.positions) {}
	// data
	const Instance& instance;
	std::vector<double> cum;
	trail<double> tr;
	std::vector<int> sensorPositions;
	// API
	// eval* methods use the trail
	// append/set/remove method do not
	void initWithSolution(const Solution& sol)
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
	void appendSensor(Solution& sol, int p)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		sol.sensors.push_back(p);
		addSensorContribution(sol, p, false);
		filterPositions(instance, p, instance.minDist, sensorPositions); //< or std::max(instance.minDist, instance.minRho)?
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
	// evaluate objective if we added a sensor in position p to sol (but without actually adding it)
	double evalAddSensor(Solution& sol, int p)
	{
		double cov = 0.0;
		double sig = 0.0;
		double obj = 0.0;
		addSensorContribution(sol, p, true);
		instance.evalObjectives(cum, cov, sig, obj);
		tr.restore();
		return obj;
	}
	// Remove the sensor in position of index i form sol
	void removeSensor(Solution& sol, int i)
	{
		int p = sol.sensors[i];
		// temporarily mark p as invalid
		sol.sensors[i] = -1;
		// remove contribution of source p (without trail!)
		removeSensorContribution(sol, p, false);
	}
	// Set sensor in position of index i to p
	void setSensor(Solution& sol, int i, int p)
	{
		DOMINIQS_ASSERT( sol.sensors[i] == -1 );
		sol.sensors[i] = p;
		addSensorContribution(sol, p, false);
	}
protected:
	// helpers
	void addSensorContribution(const Solution& sol, int p, bool useTrail)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		for (int q: sol.sensors)
		{
			if (q == -1) continue; // skip invalid detectors
			if (q == p) continue;
			PMDFPtr pmdf = instance.getPMDF(p, q);
			if (!pmdf) continue;
			if (useTrail) for (const PMDF::Item& item: pmdf->data) tr.set(item.voxel, cum[item.voxel] + item.value);
			else for (const PMDF::Item& item: pmdf->data) cum[item.voxel] += item.value;
		}
	}
	void removeSensorContribution(const Solution& sol, int p, bool useTrail)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		for (int q: sol.sensors)
		{
			if (q == -1) continue; // skip invalid detectors
			if (q == p) continue;
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
		: instance(inst), deltas(inst.gdata.positions.size(), 0), sensorPositions(inst.positions) {}
	// data
	const Instance& instance;
	std::vector<double> deltas;
	std::vector<int> sensorPositions;
	double currentObjective = 0.0;
	// API
	void initWithSolution(const Solution& sol)
	{
		// set up deltas
		std::fill(deltas.begin(), deltas.end(), 0);
		// first delta for sources
		for (int p: instance.positions)
		{
			double term = 0.0;
			for (int q: sol.sensors)
			{
				if (q == -1) continue; // skip invalid detectors
				if (q == p) continue; // skip invalid detectors
				PMDFPtr pmdf = instance.getPMDF(p, q);
				if (!pmdf) continue;
				term += pmdf->sum;
			}
			deltas[p] = term;
		}
		// compute current objective
		currentObjective = 0.0;
		for (int p: sol.sensors)
		{
			if (p == -1) continue; // skip invalid sources
			for (int q: sol.sensors)
			{
				if (q == -1) continue; // skip invalid detectors
				if (q == p) continue;
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
	void appendSensor(Solution& sol, int p)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		sol.sensors.push_back(p);
		addSensorContribution(sol, p);
		filterPositions(instance, p, instance.minDist, sensorPositions); //< or std::max(instance.minDist, instance.minRho)?
	}
	void clearSensors(Solution& sol)
	{
		sol.sensors.clear();
		currentObjective = 0.0;
		std::fill(deltas.begin(), deltas.end(), 0.0);
	}
	// evaluate objective if we added a given PMDF
	double evalPMDF(PMDFPtr pmdf)
	{
		return pmdf->sum;
	}
	// evaluate objective if we added a sensor in position p to sol (but without actually adding it)
	double evalAddSensor(Solution& sol, int p)
	{
		addSensorContribution(sol, p);
		double obj = currentObjective;
		removeSensorContribution(sol, p);
		return obj;
	}
	// Remove the sensor in position of index i from sol
	void removeSensor(Solution& sol, int i)
	{
		int p = sol.sensors[i];
		// temporarily mark p as invalid
		sol.sensors[i] = -1;
		// remove contribution of sensor p
		removeSensorContribution(sol, p);
	}
	// Set sensor in position of index i to p
	void setSensor(Solution& sol, int i, int p)
	{
		DOMINIQS_ASSERT( sol.sensors[i] == -1 );
		sol.sensors[i] = p;
		addSensorContribution(sol, p);
	}
protected:
	// helpers
	void addSensorContribution(const Solution& sol, int p)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		// update objective
		for (int q: sol.sensors)
		{
			if (q == -1) continue;
			if (q  == p) continue;
			PMDFPtr pmdf = instance.getPMDF(p, q);
			if (!pmdf) continue;
			currentObjective += pmdf->sum;
		}
		// update deltas
		for (int k = instance.start[p]; k < instance.start[p+1]; k++)
		{
			PMDFPtr pmdf = instance.pmdfs[k];
			DOMINIQS_ASSERT( pmdf );
			deltas[pmdf->q] += pmdf->sum;
		}
	}
	void removeSensorContribution(const Solution& sol, int p)
	{
		DOMINIQS_ASSERT( &instance == &(sol.instance) );
		// update objective
		for (int q: sol.sensors)
		{
			if (q == -1) continue;
			if (q == p) continue;
			PMDFPtr pmdf = instance.getPMDF(p, q);
			if (!pmdf) continue;
			currentObjective -= pmdf->sum;
		}
		// update deltas
		for (int k = instance.start[p]; k < instance.start[p+1]; k++)
		{
			PMDFPtr pmdf = instance.pmdfs[k];
			DOMINIQS_ASSERT( pmdf );
			deltas[pmdf->q] -= pmdf->sum;
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


static void graspIteration(const Instance& inst, Solution& sol, size_t maxCand, size_t mult, std::mt19937& rnd)
{
	// setup
	sol.reset();
	LocalDataWeightedObj d(inst);

	// pick initial pair of sensors
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
	d.appendSensor(sol, p);
	d.appendSensor(sol, q);
	DOMINIQS_ASSERT(sol.checkSignal(d.cum));
	DOMINIQS_ASSERT(sol.checkCompatibilities());
	// now we need to add one source and one detector at the time
	std::vector<PosObj> candlist;
	while(sol.sensors.size() < inst.nSources)
	{
		candlist.clear();
		for (int p: d.sensorPositions) candlist.emplace_back(p, d.evalAddSensor(sol, p));
		int bestp = pickCandidate(candlist, maxCand, rnd);
		if (bestp == -1) break;
		// bestp is now the new sensor to add
		d.appendSensor(sol, bestp);
		DOMINIQS_ASSERT(sol.checkSignal(d.cum));
		DOMINIQS_ASSERT(sol.checkCompatibilities());
	}
	sol.eval();
}


template <typename UserData>
inline double getSolutionObjValue(const Solution& sol, const UserData& ud)
{
	return ud.maxSignal ? -sol.signal : -sol.objective;
}


void grasp(const Instance& inst, Solution& sol, TerminationCriteria& toStop, bool optTotSignal, size_t maxCand, int seed)
{
	std::mt19937 rnd(seed);
	int iter = 0;
	int noImproveCount = 0;
	UserData ud;
	ud.maxSignal = optTotSignal;

	while (true)
	{
		if (toStop(iter, noImproveCount)) break;
		Solution currsol(inst);
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


static void computeCandidates(const Instance& inst, const Solution& sol, std::vector<int>& sensorPositions)
{
	sensorPositions = inst.positions;
	// cannot be too close to other sensors
	for (int p: sol.sensors)
	{
		if (p == -1) continue;
		filterPositions(inst, p, inst.minDist, sensorPositions); //< or std::max(inst.minDist, inst.minRho)?
	}
}


template<typename Observer, typename LocalData>
bool oneOpt(Solution& sol, Observer& rec)
{
	const Instance& inst = sol.instance;
	LocalData d(inst);
	d.initWithSolution(sol);
	double oldobj = d.getCurrentObjective();
	bool improved = false;
	// look for improving p
	size_t nsensors = sol.sensors.size();
	for (size_t i = 0; i < nsensors; i++)
	{
		int p = sol.sensors[i];
		d.removeSensor(sol, i);
		// compute set of candidate positions
		std::vector<int> candidates;
		computeCandidates(inst, sol, candidates);
		// try each of the candidates
		int imprp = p;
		for (int newp: candidates)
		{
			if (newp == p) continue;
			// eval contribution of source newp
			double obj = d.evalAddSensor(sol, newp);
			if (obj > oldobj + OBJEPS)
			{
				// we have found an improving move
				imprp = newp;
				improved = true;
				break;
			}
		}
		// apply move (could be the no move at all if imprp == p)
		d.setSensor(sol, i, imprp);
		if (improved)
		{
			sol.eval();
			break;
		}
	}
	DOMINIQS_ASSERT(sol.checkCompatibilities());
	return improved;
}


template<typename Observer, typename UserData=EmptyUserData>
bool exploreNeighborhood(Solution& sol, Observer& rec, const UserData& ud)
{
	bool improved = false;
	if (ud.maxSignal)
	{
		improved = oneOpt<Observer,LocalDataMaxSignal>(sol, rec);
	}
	else
	{
		improved = oneOpt<Observer,LocalDataWeightedObj>(sol, rec);
	}
	return improved;
}


void localSearch(const Instance& inst, Solution& sol, bool optTotSignal)
{
	IObserver<Solution> rec; //< does nothing
	UserData ud;
	ud.maxSignal = optTotSignal;
	localSearch(sol, rec, ud);
}


template <typename RandomEngine, typename UserData>
void randomGenerateSolution(Solution& sol, RandomEngine& engine, const UserData& ud)
{
	const Instance& inst = sol.instance;
	LocalDataWeightedObj d(inst);
	sol.reset();
	while(sol.sensors.size() < inst.nSources)
	{
		if (d.sensorPositions.empty()) break;
		std::uniform_int_distribution<> rndint(0, d.sensorPositions.size()-1);
		int elem = rndint(engine);
		int p = d.sensorPositions[elem];
		d.appendSensor(sol, p);
		DOMINIQS_ASSERT(sol.checkCompatibilities());
	}
	sol.eval();
}


void rrls(const Instance& inst, Solution& sol, TerminationCriteria& toStop, bool optTotSignal, int seed)
{
	IObserver<Solution> rec; //< does nothing
	std::mt19937 rnd(seed);
	UserData ud;
	ud.maxSignal = optTotSignal;
	randomRestartLocalSearch(sol, rnd, rec, toStop, ud);
}


template<typename RandomEngine, typename UserData>
void perturbeSolution(Solution& sol, RandomEngine& engine, int K, const UserData& ud)
{
	const Instance& inst = sol.instance;
	for (int k = 0; k < K; k++)
	{
		// pick element to change
		DOMINIQS_ASSERT(sol.sensors.size());
		std::uniform_int_distribution<> rndint(0, sol.sensors.size()-1);
		int i = rndint(engine);
		DOMINIQS_ASSERT(i >= 0 && i < sol.sensors.size());
		int p = sol.sensors[i];
		std::vector<int> candidates(inst.positions);
		// cannot be too close to remaining sensors
		for (int q: sol.sensors)
		{
			if (q != p)
			{
				filterPositions(inst, q, inst.minDist, candidates);
			}
		}
		// pick new value at random
		std::uniform_int_distribution<> rndint2(0, candidates.size()-1);
		int newp = candidates[rndint2(engine)];
		sol.sensors[i] = newp;
		DOMINIQS_ASSERT(sol.checkCompatibilities());
	}
	sol.eval();
}


void ils(const Instance& inst, Solution& sol, TerminationCriteria& toStop, bool optTotSignal, int seed)
{
	IObserver<Solution> rec; //< does nothing
	std::mt19937 rnd(seed);
	// HillClimbingCriterion<Solution> accept;
	RandomWalkCriterion<Solution> accept;
	UserData ud;
	ud.maxSignal = optTotSignal;
	iteratedLocalSearch(sol, accept, rnd, rec, toStop,
								std::max(3, (int)inst.nSources),
								std::max(5, (int)inst.nSources),
								ud);
}


void multiILS(const Instance& inst, Solution& sol, TerminationCriteria& toStop, bool optTotSignal, int maxTries, int seed)
{
	IObserver<Solution> rec; //< does nothing
	std::mt19937 rnd(seed);
	Solution currSol(sol.instance);
	UserData ud;
	ud.maxSignal = optTotSignal;
	for (int t = 0; t < maxTries; t++)
	{
		RandomWalkCriterion<Solution> accept;
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
