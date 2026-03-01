#ifndef LBU_CG_HPP
#define LBU_CG_HPP

/*******************************************************************************

  AUTHOR:  JA du Preez
  DATE:    November 2013
  PURPOSE: PGM cluster graph inference via residual Lauritzen-Spiegelhalter.

  (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// patrec headers
#include "prlite_vector.hpp"

// emdw headers
#include "factor.hpp"
#include "clustergraph.hpp"
#include "messagequeue.hpp"

// standard headers
#include <string>  // string
#include <iosfwd>  // istream, ostream
#include <vector>
#include <map>
#include <set>
#include <utility> // pair

/**
 * @brief Cluster graph inference via residual
 * Lauritzen-Spiegelhalter.
 *
 * @param cg: The clustergraph, basically consisting of factors, links
 * between them and a sepset associated with each link.
 *
 * @param ssBeliefs: On entry this specifies any known sepset
 * beliefs. Any others not specified (which may be all of them), will
 * be inferred automatically. On completion the (hopefully) calibrated
 * sepset beliefs are returned in this variable. NOTE: the key always
 * has sourceNr < destNr.
 *
 * @param msgQueue: On entry specifies a priority queue of *existing*
 * beliefs to *propagate* next. This will be automatically
 * supplemented with beliefs arising out of changes in the graph,
 * whether from observing variables or inferring new beliefs. On exit
 * it lists any unconverged beliefs, i.e. on convergence it should be
 * empty.
 *
 * @param dampingFactor: In range [0..1). The higher, the heavier the
 * smoothing of messages. Typically between 0.33 and 0.5, for trees
 * (i.e. no loops) this should be 0.
 *
 * @param deltaThresh: This controls the convergence criterium. When
 * the distance between two subsequent messages is less than this,
 * that message is considered converged. If not converged it will be
 * added to the message queue for further processing. A small positive
 * fraction with the given default probably the smallest you would
 * want to go.
 *
 * @param maxNoOfAbsorbs: When the system has sent this many messages
 * per cluster, it will terminate regardless of whether the message
 * queue is empty. At its default value a very large number is
 * pre-specified for this. Also useful if you want to terminate
 * pre-maturely to inspect progress.
 *
 * @param observed: Lists variables for which, in addition to what is
 * already known in the cluster graph, further (oracle) observations
 * are available. NOTE that observations actually change the
 * underlying factor graph. Heisenberg.
 *
 * NOTE: Either use find or consult linkedTo in the cluster graph to
 * see which valid beliefs exist, if you directly poke around in the
 * beliefs map asking about non-existent beliefs, the map inserts
 * all sorts of unwanted stuff in there.
 */
unsigned
loopyBU_CG(
  ClusterGraph& cg,
  std::map< Idx2, rcptr<Factor> >& ssBeliefs,
  MessageQueue& msgQueue ,
  double dampingFactor = 0.0,
  double deltaThresh = 1E-11,
  unsigned maxNoOfAbsorbs = 0,
  const std::map<emdw::RVIdType, AnyType>& observed = std::map<emdw::RVIdType, AnyType>());

/**
 * No nonsense version of the above
 *
 */
unsigned
loopyBU_CG(
    ClusterGraph& cg,
    double dampingFactor = 0.0,
    double deltaThresh = 1E-11,
    unsigned maxNoOfAbsorbs = 0,
    const std::map<emdw::RVIdType, AnyType>& observed = std::map<emdw::RVIdType, AnyType>()
);


/**
 * The order of the variables in queryVars is unimportant
 *
 */
uniqptr<Factor>
queryLBU_CG(
  ClusterGraph& cg,
  std::map< Idx2, rcptr<Factor> >& ssBeliefs,
  const emdw::RVIds& queryVars);


#endif // LBU_CG_HPP
