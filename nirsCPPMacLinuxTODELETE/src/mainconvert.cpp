/**
 * @file main.cpp
 * @brief
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2018 Domenico Salvagnin
 */


#include "nirs/nirs.h"
#include "nirs/nirsheur.h"
#include "nirs/str_utils.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace dominiqs;

static const unsigned int MIN_NUMBER_OF_INPUTS = 2;

static const std::string scalpFile("../data/down3/scalpPos.txt");
static const std::string voxelFile("../data/down3/GMPos.txt");
static const std::string allPMDFfile("../data/down3/all_pmdfs.txt");
static const std::string wFile("../data/PMDF_MASK_30-60mmWeight_ToastCorrect.txt");


void exportData(const NIRSData& data, const char* fileidx, const char* filedata)
{
	std::ofstream outdata(filedata, std::ios::binary);
	size_t numpositions = data.positions.size();
	size_t numitems = 0;
	// write data file
	for (size_t p = 0; p < numpositions; p++)
	{
		for (size_t q = p+1; q < numpositions; q++)
		{
			PMDFPtr pmdf = data.allpmdfs[p*numpositions+q];
			if (!pmdf)  continue;
			for (const auto& item: pmdf->data)
			{
				int32_t vox = item.voxel;
				outdata.write((const char*)(&vox), sizeof(vox));
				outdata.write((const char*)(&item.value), sizeof(item.value));
				numitems++;
			}
		}
	}
	outdata.close();
	// write index file
	std::ofstream outidx(fileidx);
	for (size_t p = 0; p < numpositions; p++)
	{
		for (size_t q = p+1; q < numpositions; q++)
		{
			size_t counts = 0;
			PMDFPtr pmdf = data.allpmdfs[p*numpositions+q];
			if (pmdf)  counts = pmdf->data.size();
			outidx << counts << " ";
		}
		outidx << std::endl;
	}
	outidx.close();
	std::cout << "Wrote " << numitems << " items to disk" << std::endl;
}


int main (int argc, char const *argv[])
{
	if (argc <= MIN_NUMBER_OF_INPUTS)
	{
		std::cout << "usage: nirsconvert fileidx filedata" << std::endl;
		return -1;
	}
	NIRSData nirsdata;
	nirsdata.read(scalpFile, voxelFile, allPMDFfile, wFile);
	exportData(nirsdata, argv[1], argv[2]);
	return 0;
}
