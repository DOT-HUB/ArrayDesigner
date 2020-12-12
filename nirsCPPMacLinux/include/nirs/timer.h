/**
 * @file timer.h
 * @brief CPU Stopwatch for benchmarks
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2019
 */

#ifndef TIMER_H
#define TIMER_H

#include "asserter.h"
#include <string>
#include <chrono>


namespace dominiqs {

/**
 * Benchmarking class
 * Simulates a simple stopwatch for wallclock time
 */

class StopWatch
{
public:
	/** default constructor */
	StopWatch(bool autoStart = false);
	/** start stopwatch */
	void start();
	/** stop stopwatch */
	void stop();
	/** reset stopwatch */
	void reset();

	/** @return partial elapsed time in seconds (according to default type)*/
	double getPartial() const;
	/** @return total elapsed time in seconds (according to default type)*/
	double getTotal() const;
	/** @return elapsed time in seconds (according to default type)*/
	double getElapsed() const;
private:
	using HRClock = std::chrono::high_resolution_clock;
	using Seconds = std::chrono::duration<double>;
	HRClock::time_point w_begin;
	HRClock::time_point w_end;
	double w_total;
};

inline StopWatch::StopWatch(bool autoStart) : w_total(0.0)
{
	if (autoStart) start();
}

inline void StopWatch::start()
{
	w_begin = HRClock::now();
}

inline void StopWatch::stop()
{
	w_end = HRClock::now();
	w_total += std::chrono::duration_cast<Seconds>(w_end - w_begin).count();
}

inline void StopWatch::reset()
{
	w_total = 0.0;
}

inline double StopWatch::getPartial() const
{
	return std::chrono::duration_cast<Seconds>(w_end - w_begin).count();
}

inline double StopWatch::getTotal() const
{
	return w_total;
}

inline double StopWatch::getElapsed() const
{
	auto w_now = HRClock::now();
	return std::chrono::duration_cast<Seconds>(w_now - w_begin).count();
}

/**
 * Global Stopwatch
 */

StopWatch& gStopWatch();

/**
 * Automatic stopwatch stopper (RAII principle)
 */

class AutoStopWatch
{
public:
	AutoStopWatch(StopWatch* c) : theStopWatch(c), pending(true)
	{
		DOMINIQS_ASSERT( theStopWatch );
		theStopWatch->start();
	}
	~AutoStopWatch() { stop(); }
	void stop()
	{
		if (pending)
		{
			theStopWatch->stop();
			pending = false;
		}
	}
private:
	StopWatch* theStopWatch;
	bool pending;
};

/**
 * Automatic global stopwatch stopper (RAII principle)
 */

class GlobalAutoStopWatch
{
public:
	GlobalAutoStopWatch() : pending(true) { gStopWatch().start(); }
	~GlobalAutoStopWatch() { stop(); }
	void stop()
	{
		if (pending)
		{
			gStopWatch().stop();
			pending = false;
		}
	}
private:
	bool pending;
};

/**
 * Get current date/time as std::string, format is YYYY-MM-DD.HH:mm:ss
 */
std::string currentDateTime();

} // namespace dominiqs

#endif /* TIMER_H */
