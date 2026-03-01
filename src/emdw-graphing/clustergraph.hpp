#ifndef CLUSTERGRAPH_HPP
#define CLUSTERGRAPH_HPP

/*******************************************************************************

          AUTHOR:    JA du Preez, SF Streicher
          COPYRIGHT: University of Stellenbosch, all rights reserved.
          DATE:      March 2013 - 2018
          PURPOSE:   ClusterGraph PGM cluster graph

*******************************************************************************/

// patrec headers
#include "prlite_vector.hpp"

// emdw headers
#include "factor.hpp"
#include "vecset.hpp"

// standard headers
#include <string>  // string
#include <iosfwd>  // istream, ostream
#include <vector>
#include <map>
#include <set>
#include <utility> // pair


//
/******************************************************************************/

/**
 *  Code to build an RIP cluster graph.
 *
 */

class ClusterGraph {

  //============================ Traits and Typedefs ============================
public:

  class _LTRIP_{}; static _LTRIP_ LTRIP;
  class _BETHE_{}; static _BETHE_ BETHE;
  class _JTREE_{}; static _JTREE_ JTREE;
  class _LPRODS_{}; static _LPRODS_ LPRODS;
  class _LBAGS_{}; static _LBAGS_ LBAGS;
  class _LJTREES_{}; static _LJTREES_ LJTREES;
  class _RLJTREES_{}; static _RLJTREES_ RLJTREES;

  //======================== Constructors and Destructor ========================
  public:

  /**
   * @brief Assembles an RIP cluster graph via the LTRIP algorithm.
   *
   * NOTE: for details on the general automatic construction of RIP
   * cluster graphs, see S.F. Streicher and J.A. du Preez. Graph
   * coloring: Comparing cluster graphs to factor graphs. SAWACMMM
   * '17- Proceedings of the ACM Multimedia 2017 Workshop on South
   * African Academic Participation, Mountain View CA, October 2017.
   *
   * @param theFactorPtrs The PGM factors/potentials.
   *
   * @param observed: If some observations are available, it will be
   * used to simplify the graph. Typically these are things like pixel
   * values in an image, feature vectors in speech processing etc.
   *
   * @param deepCopy: If true the factors internal to the PGM are
   * distinct from theFactorPtrs (i.e. deep copies). If false, shallow
   * copies might be used for efficiency. NOTE: With shallow copies
   * the original factors may be changed in unpredictable ways.
   *
   * @param theMassOfCluster Gives the informedness of each
   * factor/cluster. If empty it is automatically determined via
   * "distanceFromVacuous".
   */
  ClusterGraph(
    _LTRIP_,
    const std::vector< rcptr<Factor> >& theFactorPtrs,
    const std::map<emdw::RVIdType, AnyType>& observed = std::map<emdw::RVIdType, AnyType>(),
    const bool deepCopy = true,
    std::map<size_t,double> theMassOfCluster = std::map<size_t,double>() );

  /**
   * @brief Assembles a Bethe cluster graph (i.e. classical factor graph).
   *
   * @param theFactorPtrs The PGM factors/potentials.
   *
   * @param observed: If some observations are available, it will be
   * used to simplify the graph. Typically these are things like pixel
   * values in an image, feature vectors in speech processing etc.
   *
   * @param deepCopy: If true the factors internal to the PGM are
   * distinct from theFactorPtrs (i.e. deep copies). If false, shallow
   * copies might be used for efficiency. NOTE: With shallow copies
   * the original factors may be changed in unpredictable ways.
   */
  ClusterGraph(
    _BETHE_,
    const std::vector< rcptr<Factor> >& theFactorPtrs,
    const std::map<emdw::RVIdType, AnyType>& observed = std::map<emdw::RVIdType, AnyType>(),
    const bool deepCopy = true);

  /**
   * Assembles a Junction Tree structured cluster graph (also known as
   * a Clique Tree).
   *
   * @param theFactorPtrs: The PGM factors/potentials.
   *
   * @param observed: If some observations are available, it will be
   * used to simplify the graph. Typically these are things like pixel
   * values in an image, feature vectors in speech processing etc.
   *
   * @param deepCopy: If true the factors internal to the PGM are
   * distinct from theFactorPtrs (i.e. deep copies). If false, shallow
   * copies might be used for efficiency. NOTE: With shallow copies
   * the original factors may be changed in unpredictable ways.
   */
  ClusterGraph(
    _JTREE_,
    const std::vector< rcptr<Factor> >& theFactorPtrs,
    const std::map<emdw::RVIdType, AnyType>& observed = std::map<emdw::RVIdType, AnyType>(),
    const bool deepCopy = true);

  ClusterGraph(
    _LPRODS_,
    const std::vector< rcptr<Factor> >& theFactorPtrs,
    const std::map<emdw::RVIdType, AnyType>& observed = std::map<emdw::RVIdType, AnyType>(),
    const bool deepCopy = true,
    size_t maxSz = 16);

  ClusterGraph(
    _LJTREES_,
    const std::vector< rcptr<Factor> >& theFactorPtrs,
    const std::map<emdw::RVIdType, AnyType>& observed = std::map<emdw::RVIdType, AnyType>(),
    const bool deepCopy = true,
    size_t maxSz = 16);

  // en dan wag ons nog vir _RLJTREES_

/**
 * @brief Class Specific constructor. Makes shallow copies of the Factors.
 *
 * @param theFactorPtrs: The PGM factors/potentials.
 *
 * @param theLinkedTo: To which factors is each factor linked
 *
 * @param theSepvecs: The sepset between each pair of linked factors.
 *
 * @param theFactorsWithRV Which factors contain a particular RV. If
 * empty this will be recalculated.
 *
 * @param deepCopy: If true the factors internal to the PGM are
 * distinct from theFactorPtrs (i.e. deep copies). If false, shallow
 * copies might be used for efficiency. NOTE: With shallow copies the
 * original factors may be changed in unpredictable ways.
 */
  ClusterGraph(
    const std::vector< rcptr<Factor> >& theFactorPtrs,
    const std::vector< std::set<unsigned> >& theLinkedTo,
    const std::map<Idx2, emdw::RVIds>& theSepvecs,
    const std::map< emdw::RVIdType, std::vector<unsigned> >& theFactorsWithRV
    = std::map< emdw::RVIdType, std::vector<unsigned> >(),
    const bool deepCopy = true );

/**
   * Kept for legacy reasons. It defers to the LTRIP version above.
   */
  ClusterGraph(
    const std::vector< rcptr<Factor> >& theFactorPtrs,
    const std::map<emdw::RVIdType, AnyType>& observed = std::map<emdw::RVIdType, AnyType>(),
    const bool deepCopy = true);

  ///
  ClusterGraph();

  ///
  ClusterGraph(const ClusterGraph& cg) = default;
  ClusterGraph(ClusterGraph&& cg) = default;

  ///
  virtual ~ClusterGraph();

  //========================= Operators and Conversions =========================
public:

  ///
  ClusterGraph& operator=(const ClusterGraph& d) = default;
  ClusterGraph& operator=(ClusterGraph&& d) = default;

  //================================= Iterators =================================

  //============================= Exemplar Support ==============================
public:

  ///
  virtual ClusterGraph* copy() const;

  ///
  virtual const std::string& isA() const;

  //=========================== Inplace Configuration ===========================
public:

  // class specific
  unsigned classSpecificConfigure(
    const std::vector< rcptr<Factor> >& theFactorPtrs,
    const std::map<emdw::RVIdType, AnyType>& observed,
    const bool deepCopy);

  ///
  virtual std::istream& txtRead(std::istream& file);
  ///
  virtual std::ostream& txtWrite(std::ostream& file) const;

  //===================== Required Virtual Member Functions =====================

  //====================== Other Virtual Member Functions =======================

  //======================= Non-virtual Member Functions ========================
  unsigned noOfFactors() const;

  /**
   * @brief Exports the cluster graph to a graphviz dot file - this is
   *   useful for visually inspecting your graph structure.
   *
   * @param graphname: Name of the graph. The corresponding graphviz
   *   filename will be <graphname>.dot
   *
   * @param idsToShow: Only show cluster nodes that also contains at
   * least one of these variable ids. If empty, the full graph will be
   * shown.
   */
  void
  exportToGraphViz(
    const std::string& graphname,
    const std::set<emdw::RVIdType>& idsToShow = std::set<emdw::RVIdType>() ) const;

  void
  plotNodeStrengths(
    const std::string& graphname) const;

  /**
   * A static helper member that applies the observations to all the
   * applicable factors while also removing RV's that become redundant.
   * If a factor is fully observed (i.e. it has become empty and
   * therefore redundant), it will be removed altogether from the list
   * of factors.
   *
   * @param factorsWithRV maps all the RV's against the factors
   * containing that RV. If it is empty it will be calculated by this
   * routine. On exit it will return a map for the remaining RVs.
   *
   * NOTE: this routine may change the numbering and order of factors.
   *
   * ISSUE: in this code we are assuming, but not enforcing that the
   * supplied FactorOperator will observe and there-after remove the
   * observed variables from the Factor scopes.
   */
  static void observeAndReduce(
    std::vector< rcptr<Factor> >& theFactorPtrs,
    const std::map<emdw::RVIdType, AnyType>& observed,
    std::map< emdw::RVIdType, std::vector<unsigned> >& theFactorsWithRV);

  static void swapFactors(
    std::vector< rcptr<Factor> >& factors,
    size_t idx1,
    size_t idx2,
    std::vector<size_t>& toNew,
    std::vector<size_t>& toOld);

/**
 * This routine typically is the first inspection / sanitising for a
 * set of factors before building them into a PGM. It gets rid of
 * subset factors, while recording other useful bits of information
 * that it had to calculate for this in any case.
 *
 * Any factor that is a subset of another factor, and also is not
 * marked as "keepSeparate", gets absorbed into this other factor and
 * thereafter removed from the list. We also record the factor
 * intersections (forward looking i.e. c > r) and and which clusters
 * has variables in common with which others (all c != r). Lastly we
 * list all the factors containing a particular RV.
 *
 * @param factorsWithRV If empty on entry, determines in which factors
 * a particular RV can be found. If non-empty on entry, accept and use
 * the supplied values (i.e. no determination).
 */
static void absorbSubsetFactors(
  std::vector< rcptr<Factor> >& theFactorPtrs,
  std::vector<std::map<size_t, emdw::RVIds> >& theIntersections,
  std::map<size_t, std::set<size_t> >& theReachableFrom,
  std::map< emdw::RVIdType, std::vector<unsigned> >& theFactorsWithRV);


  /**
   * @theFactorPtrs On entry contains pointers to the collection of
   * Factors, some of which might be null. On exit returns by reference
   * a compacted version where all the null Factors are removed.
   *
   * @param theFactorsWithRV A map that lists all the factors containing
   * a particular RV. May be uninitialized on entry.
   *
   * @param old2new can be uninitialized on entry. It returns by
   * reference to where a specific factor has relocated.
   */
  static void repackFactors(
    std::vector< rcptr<Factor> >& theFactorPtrs,
    std::map< emdw::RVIdType, std::vector<unsigned> >& theFactorsWithRV,
    std::vector<unsigned>& old2new);

  static void repackGraph(
    std::vector< rcptr<Factor> >& theFactorPtrs,
    std::map< emdw::RVIdType, std::vector<unsigned> >& theFactorsWithRV,
    std::vector< std::set<unsigned> >& theLinkedTo,
    std::map<Idx2, emdw::RVIds>& theSepvecs);

private:

  /**
   * A helper member that assembles a cluster graph from factors. We
   * define a clustergraph as a set of a) factors, with b) links
   * between them and c) a sepvec associated with each link.
   *
   * @param factorsWithRV maps all the RV's against the factors
   * containing that RV. If it is empty it will be calculated by this
   * routine. On exit it will return a map showing in which factor each
   * RV occurs.
   *
   * @param theMassOfCluster is a measure specifying how informative
   * the particular cluster is. The LTRIP procedure attempts to focus
   * the factor trees around the more informative clusters. If empty,
   * it will be estimated via the distanceFromVacuous function.
   */
  void buildJunctionTree(
          std::vector<rcptr<Factor>>& theFactorPtrs,
          std::map<emdw::RVIdType, std::vector<unsigned>>& theFactorsWithRV);

  //================================== Friends ==================================
public:

  /// output
  friend
  std::ostream&
  operator<<( std::ostream& file,
              const ClusterGraph& cg );

  /// input
  //friend std::istream& operator>>( std::istream& file,
  //                                 ClusterGraph& cg );

  //=============================== Data Members ================================
public:
  // ---------- definitive members
  std::vector< rcptr<Factor> > factorPtrs;

//private:
  // --------- derivative members

  // to which factors is each factor linked (RIP via a sepset)
  std::vector< std::set<unsigned> > linkedTo;

  // the sepset between all connected factor pairs (both ways)
  std::map<Idx2, emdw::RVIds> sepvecs;

  // in which factors does an RV occur (ascending order)
  std::map< emdw::RVIdType, std::vector<unsigned> > factorsWithRV;

  //============================ DEPRECATED Members =============================

}; // ClusterGraph


//=============================== Non-Class helpers =============================

std::ostream& writeMessagesToText(
  std::ostream& file,
  const std::map< Idx2, rcptr<Factor> >& messages);

emdw::RVIds varElemOrder(const std::vector<rcptr<Factor>> &theFactorPtrs,
                         const emdw::RVIds &theVars = {});

std::vector<std::pair<emdw::RVIds, std::vector<unsigned>>>
  varMarginalizationOrder(const std::vector<rcptr<Factor>> &theFactorPtrs);

#endif // CLUSTERGRAPH_HPP
