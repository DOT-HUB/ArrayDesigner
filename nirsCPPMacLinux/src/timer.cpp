/**
 * @file timer.cpp
 * @brief CPU Stopwatch for benchmarks
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2019
 */

#include <nirs/timer.h>

namespace dominiqs {

StopWatch& gStopWatch()
{
	static StopWatch theStopWatch;
	return theStopWatch;
}

} // namespace dominiqs
