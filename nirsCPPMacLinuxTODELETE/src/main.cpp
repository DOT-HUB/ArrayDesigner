/**
 * @file main.cpp
 * @brief
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2018 Domenico Salvagnin
 */


#include "nirs/nirs.h"
#include "nirs/nirsheur.h"
#include "nirs/timer.h"
#include "nirs/str_utils.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace dominiqs;

template<class SolType>
void runHeur(const Instance& instance, SolType& sol, const std::string& name, bool maxSignal, int maxIter=100, double timeLimit=100.0)
{
	StopWatch chrono(true);
	TerminationCriteria toStop(chrono);
	toStop.iterLimit = maxIter;
	toStop.noImproveLimit = std::max(maxIter / 10, 5);
	toStop.timeLimit = timeLimit;
	if (name == "GRASP") grasp(instance, sol, toStop, maxSignal);
	else if (name == "RRLS") rrls(instance, sol, toStop, maxSignal);
	else if (name == "ILS") ils(instance, sol, toStop, maxSignal);
	else if (name == "multiILS") multiILS(instance, sol, toStop, maxSignal);
	else throw std::runtime_error("Unknown heuristic!");
	chrono.stop();
	std::cout << name << ":" << std::endl;
	sol.print();
	std::cout << "Time: " << chrono.getTotal() << std::endl;
}

static const unsigned int MIN_NUMBER_OF_INPUTS = 11;

int main (int argc, char const *argv[])
{
	static_assert(sizeof(void*) >= 8, "NIRS needs to be compiled as a 64-bit app (use cmake .. -A x64 on Windows)");
	if (argc <= MIN_NUMBER_OF_INPUTS)
	{
		std::cout << "usage: nirsmain scalpFile GMfile pmdfIDXfile pmdfDATAfile weights_file instance_file nS nD minDist coverageThr coverageWeight output_file" << std::endl;
		return -1;
	}
	try
	{
		StopWatch readSW(true);
		NIRSData nirsdata;
		std::string scalpFile = std::string(argv[1]);
		std::string voxelFile = std::string(argv[2]);
		std::string allPMDFfileIDX = std::string(argv[3]);
		std::string allPMDFfileDATA = std::string(argv[4]);
		std::string wFile = std::string(argv[5]);
		nirsdata.read(scalpFile, voxelFile, allPMDFfileIDX, allPMDFfileDATA, wFile);
		Instance instance(nirsdata);
		int s = from_string<int>(argv[7]);
		int d = from_string<int>(argv[8]);
		instance.read(s, d, argv[6]);
		instance.minDist = from_string<double>(argv[9]);
		instance.coverageThr = from_string<double>(argv[10]);
		instance.coverageWeight = from_string<double>(argv[11]);
		readSW.stop();
		std::cout << "Read time: " << readSW.getTotal() << std::endl;
#ifdef NONBIPARTITE_CASE
		Solution sol(instance);
#else
		BSolution sol(instance);
#endif /* NONBIPARTITE_CASE */
		std::cout << "First run optimizing signal only:" << std::endl;
		runHeur(instance, sol, "GRASP", true);
		// use the maximum signal as signal normalization
		instance.sigNorm = sol.signal;
		sol.eval();
		if (sol.coverage >= (double(instance.voxels.size()) - 1e-9))
		{
			std::cout << "No need for second run: optimizing signal gives full coverage :-)" << std::endl;
		}
		else
		{
			std::cout << "Second run optimizing weighted average of signal and coverage:" << std::endl;
			runHeur(instance, sol, "GRASP", false);
		}
		// Output solution
		std::ofstream out(argv[12]);
		sol.save(out);
	}
	catch(const std::exception& e)
	{
		std::cout << "Exception: " << e.what() << std::endl;
		return -1;
	}
	return 0;
}
