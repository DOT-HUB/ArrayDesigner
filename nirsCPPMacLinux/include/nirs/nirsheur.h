/**
 * @file nirsheur.h
 * @brief Heuristics for NIRS
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2017 Domenico Salvagnin
 */

#ifndef NIRSHEUR_H
#define NIRSHEUR_H

#include "nirs.h"
#include "metaheur_basics.h"

/**
 * Heuristics for NIRS - bipartite case
 */

/* GRASP heuristic */
void grasp(const Instance& inst, BSolution& sol, dominiqs::TerminationCriteria& toStop, bool optTotSignal=false, size_t maxCand=2, int seed=2017);

/* Local Search */
void localSearch(const Instance& inst, BSolution& sol, bool optTotSignal=false);

/* Random-restart LS */
void rrls(const Instance& inst, BSolution& sol, dominiqs::TerminationCriteria& toStop, bool optTotSignal=false, int seed=2017);

/* ILS */
void ils(const Instance& inst, BSolution& sol, dominiqs::TerminationCriteria& toStop, bool optTotSignal=false, int seed=2017);
void multiILS(const Instance& inst, BSolution& sol, dominiqs::TerminationCriteria& toStop, bool optTotSignal=false, int maxTries=10, int seed=2017);

/**
 * Heuristics for NIRS - nonbipartite case
 */

/* GRASP heuristic */
void grasp(const Instance& inst, Solution& sol, dominiqs::TerminationCriteria& toStop, bool optTotSignal=false, size_t maxCand=2, int seed=2017);

/* Local Search */
void localSearch(const Instance& inst, Solution& sol, bool optTotSignal=false);

/* Random-restart LS */
void rrls(const Instance& inst, Solution& sol, dominiqs::TerminationCriteria& toStop, bool optTotSignal=false, int seed=2017);

/* ILS */
void ils(const Instance& inst, Solution& sol, dominiqs::TerminationCriteria& toStop, bool optTotSignal=false, int seed=2017);
void multiILS(const Instance& inst, Solution& sol, dominiqs::TerminationCriteria& toStop, bool optTotSignal=false, int maxTries=10, int seed=2017);

#endif /* end of include guard: NIRSHEUR_H */
