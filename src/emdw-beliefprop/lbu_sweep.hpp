#ifndef LBU_SWEEP_HPP
#define LBU_SWEEP_HPP

/*
 * Author     :  (DSP Group, E&E Eng, US)
 * Created on :
 * Copyright  : University of Stellenbosch, all rights retained
 */

#include <vector>
#include <map>
#include <set>
#include <limits>

#include "emdw.hpp"
#include "clustergraph.hpp"

/**
 * @brief Calculates a message passing order that never revisits a
 * previously encountered source node.
 *
 * @param links Which set of nodes are coupled to each node.
 *
 * @param initNodes The initial nodes from where message propagation
 * should commence. Typically they will be the most informative
 * nodes.
 *
 * @return A vector of {fromCluster,toCluster} pairs to be used for
 * the forward sweep. The backward sweep simply uses this vector in
 * reversed order.
 */

std::vector< std::pair<int, int> > msgOrder(
  const std::vector< std::set<unsigned> >& links,
  const std::set<int>& initNodes);

/**
 * @brief Does full cluster sweep Lauritzen Spiegelhalter message
 * passing on an RIP clustergraph.
 *
 * @param cg: The clustergraph, basically consisting of factors, links
 * between them and a sepset associated with each link.
 *
 * @param ssBeliefs: On entry this specifies any known sepset
 * beliefs. Any others not specified (which may be all of them), will
 * be inferred automatically. On completion the (hopefully) calibrated
 * sepset beliefs are returned in this variable.
 *
 * @param forwardSweep The forward direction messages in the order
 * that they should be executed. Configuration options:
 * (a) If its length is exactly one, the specified source node will be
 * used as the originating node.
 * (b) If empty the system will first try to find all leaf nodes. If
 * that fails it will pick the cluster with greatest distance from its
 * vacuous state.
 * (c) Otherwise the full forward message sequence should be
 * specified. NOTE: Transitions to earlier source nodes are
 * prohibited. The msgOrder function can be used to calculate such
 * valid sweeps.
 *
 * @param dampingFactor: In range [0..1). The higher, the heavier the
 * smoothing of messages. Typically between 0.33 and 0.5, for trees
 * (i.e. no loops) this should be 0.
 *
 * @param deltaThresh: This controls the convergence criterium. When
 * the distance between two subsequent messages is less than this,
 * that message is considered converged. If not converged it will be
 * marked for further processing. A small positive fraction with the
 * given default probably the smallest you would want to go.
 *
 * @param maxNoOfSweeps: When the system has sent this many message
 * sweeps, it will terminate regardless of whether it has
 * converged. At its default value a very large number is
 * pre-specified for this. Also useful if you want to terminate
 * pre-maturely to inspect progress.
 *
 * @return The number of absorb/cancel operations that was done.
 */

unsigned lbu_sweep(
  ClusterGraph& cg,
  std::map< Idx2, rcptr<Factor> >& ssBeliefs,
  std::vector< std::pair<int, int> >& forwardSweep,
  double dampingFactor = 0.0,
  double distThresh = 1E-11,
  unsigned maxNoOfSweeps = std::numeric_limits<unsigned>::max());

#endif // LBU_SWEEP_HPP
