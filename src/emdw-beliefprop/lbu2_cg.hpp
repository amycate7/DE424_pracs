#ifndef LBU2_CG_HPP
#define LBU2_CG_HPP

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
 * Lauritzen-Spiegelhalter.  Results probably will differ from
 * that op BP_CG - if an absorb reduces a prob to approx zero you can
 * not recover it with a cancel. It therefore is fast, but at the same
 * time also less forgiving of modelling inconsistencies than its
 * Shafer-Shenoy counterpart.
 *
 * @param cg: The clustergraph, basically consisting of factors, links
 * between them and a sepset associated with each link.
 *
 * @param messages: On entry this specifies any known sepset
 * messages. Any others not specified (which may be all of them), will
 * be inferred automatically. On completion the (hopefully) calibrated
 * sepset messages are returned in this variable.
 *
 * @param msgQueue: On entry specifies a priority queue of *existing*
 * messages to *propagate* next. This will be automatically
 * supplemented with messages arising out of changes in the graph,
 * whether from observing variables or inferring new messages. On exit
 * it lists any unconverged messages, i.e. on convergence it should be
 * empty.
 *
 * @param dampingFactor: In range [0..1). The higher, the heavier the
 * smoothing of messages. Typically between 0.33 and 0.5, for trees
 * (i.e. no loops) this should be 0.
 *
 * @param deltaThresh: When the distance between two subsequent
 * messages is less than this, that message is considered converged. A
 * small positive fraction with the given default probably the
 * smallest you would want to go.
 *
 * @param observed: Lists variables for which, in addition to what is
 * already known in the cluster graph, further (oracle) observations
 * are available. NOTE that observations actually change the
 * underlying factor graph. Heisenberg.
 *
 * @param singlMsg: A configuration parameter, factors that do
 * approximation might want to set it to false.
 *
 * NOTE: Either use find or consult linkedTo in the cluster graph to
 * see which valid messages exist, if you directly poke around in the
 * messages map asking about non-existent messages, the map inserts
 * all sorts of unwanted stuff in there.
 */
unsigned
loopyBU2_CG(
  ClusterGraph& cg,
  std::map< Idx2, rcptr<Factor> >& messages,
  MessageQueue& msgQueue,
  double dampingFactor = 0.0,
  double deltaThresh = 1E-11,
  const std::map<emdw::RVIdType, AnyType>& observed = std::map<emdw::RVIdType, AnyType>(),
  bool singleMsg=false);

/**
 * The order of the variables in queryVars is unimportant
 *
 */
uniqptr<Factor>
queryLBU2_CG(
  ClusterGraph& cg,
  std::map< Idx2, rcptr<Factor> >& messages,
  const emdw::RVIds& queryVars);


#endif // LBU2_CG_HPP
