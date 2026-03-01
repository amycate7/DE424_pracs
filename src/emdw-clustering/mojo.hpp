#ifndef MOJO_HPP
#define MOJO_HPP

/*******************************************************************************

  AUTHOR:  JA du Preez
  DATE:    November 2013
  PURPOSE: A priority queue for PGM beliefs/distributions/messages.

  (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// emdw headers
#include "emdw.hpp"
#include "factor.hpp"

// standard headers
#include <vector>
#include <map>
#include <set>
#include <cstdlib>

/******************************************************************************/

/******************************************************************************/

/**
 * This class records and maintains the "mojo" pull an (external)
 * cluster (of factors) will have towards various destination clusters
 * recorded here. In particular it gives fast access to the
 * destination cluster with the strongest pull, while also efficiently
 * managing the other linked clusters.
 *
 * Access to toC[n] in log2(N). Direct access to max mojo
 * i.e. O(const)
 */
class Mojo{
  friend class MojoQueue;
  using Itr = std::map<size_t, double>::iterator;
  using ConstItr = std::map<size_t, double>::const_iterator;

public:
  // Need the initializer to stick this into a map. But the Mojo must
  // always be properly initialised with a valid "from" value.
  Mojo( size_t theFromC = std::numeric_limits<size_t>::max() );
  Mojo(const Mojo& gr);
  Mojo(Mojo&& gr);
  Mojo& operator=(const Mojo& d);
  Mojo& operator=(Mojo&& d);

  void insert(size_t toCluster, double sc);
  void erase(size_t key);
  // ONLY if there is an existing toOld, it is replaced with a
  // (toNew,scNew), even if there was an existing toNew.
  void replace(size_t toOld, size_t toNew, double scNew);
  void popMax();
  size_t idx1() const;
  size_t idx2() const;
  double strength() const;
  size_t size() const;
  void show(unsigned level = 1) const;

  ConstItr begin() const;
  ConstItr end() const;
  ConstItr find(size_t key) const;

  //const Itr find(size_t key) {return toC.find(key);}
  //Itr begin() {return toC.begin();}
  //Itr end() {return toC.end();}

private:
  void setMaxItr();
  size_t fromC;
  std::map<size_t, double> toC;
  std::map<size_t, double>::iterator maxItr;
}; // class Mojo

/******************************************************************************/

/**
 * Access to mojo[n] in log2(N). Direct access to max mojo
 * i.e. O(const).
 */
class MojoQueue{

public:

  size_t size() const {return mojoQ.size();} // size

  bool exists(size_t idx1, size_t idx2) const {
    auto itrQ = locate.find(idx1);
    if ( itrQ == locate.end() ) {return false;} // if
    auto itrG = itrQ->second->find(idx2);
    if ( itrG == itrQ->second->end() ) {return false;} // if
    return true;
  } // exists

  size_t maxIdx1() const {return mojoQ.begin()->idx1();} // maxIdx1
  size_t maxIdx2() const {return mojoQ.begin()->idx2();} // maxIdx2

  void insert(const Mojo& gr) {
    locate[gr.idx1()] = mojoQ.insert(gr);
  } // insert

  Mojo popMax() {
    Mojo ret( std::move( *mojoQ.begin() ) );
    //Mojo ret( ( *mojoQ.begin() ) );
    locate.erase( ret.idx1() );
    mojoQ.erase( mojoQ.begin() );
    return ret;
  } // popMax

  Mojo popIdx1(size_t idx1) {
    auto idx1Itr = locate[idx1];
    Mojo ret( std::move( *idx1Itr ) );
    //Mojo ret( ( *idx1Itr ) );
    locate.erase(idx1);
    mojoQ.erase(idx1Itr);
    return ret;
  } // popIdx1

  void show(unsigned level) const {
    for (auto el: mojoQ) {
      el.show(level);
    } // for
  } // show

  /// @brief helper functor to compare Mojo strengths.
  class MojoCompare { // compares two messages to find smaller
  public:
    // returns true if g1.strength() > g2.strength()
    bool
    operator() (
      const Mojo& g1,
      const Mojo& g2) const {
      return g1.strength() > g2.strength();
    } // operator()
  }; // MojoCompare

private:
  std::multiset<Mojo,MojoCompare> mojoQ;
  std::map<size_t, std::multiset<Mojo>::iterator> locate;

}; // MojoQueue

/******************************************************************************/

/**
 * @brief Calculates mojo between factors based on a gravity equations
 * analogy, while also returing several fields that will be useful for
 * other purposes such as clustering factors.
 *
 * @param factors The input factors
 *
 * @param theIntersections[m][n] (with n > m) supplies the RV intersection
 * between factors m and n. NOTE: one-sided.
 *
 * @param theMassOfCluster Supplies the mass for each factor
 * (typically calculated via distanceFromVacuous).
 *
 * @return The mojos between factors. Metaphorically map<size_t, Mojo>
 * gravityMap[r][c] is the accelaration pulling the factor c towards
 * factor r. The accelaration is calculated as the sum of a primary
 * and a secondary accelaration.  NOTE that gravityMap[r][c] !=
 * gravityMap[c][r]. The primary accelaration depends only on the
 * "Jaccard-distance" between the factors, and their individual
 * "masses". The secondary accelaration describes the gravitational
 * pull of other factors that lies "beyond" factor r (i.e. factors
 * reachable from factor r, with some RV overlap with both factor r
 * and factor c).
 */
std::map<size_t, Mojo>
calcGravityMojos(
  const std::vector< rcptr<Factor> >& factors,
  const std::vector< std::map<size_t, emdw::RVIds> >& theIntersections,
  const std::map<size_t, std::set<size_t> >& theReachableFrom,
  const std::map<size_t,double>& theMassOfCluster,
  bool includeSecondary = true);

/******************************************************************************/

void clusterGravity(const std::map<size_t, Mojo>& mojoMap,
                    std::map<size_t, std::set<size_t> >& theReachableFrom,
                    std::map<size_t,double>& massOfCluster,
                    std::map< size_t, std::set<size_t> >& factorsInCluster,
                    std::map< size_t, std::set<emdw::RVIdType> >& varsInCluster,
                    size_t maxRVs);

/******************************************************************************/

/**
 * @brief Uses Prim's algorithm to build a maximum spanning tree.
 *
 * @param available ON ENTRY lists the connection strengths between
 * pairs of nodes (there must be at least two nodes). ON EXIT: The
 * nodes NOT included in the MST (no positive connection strength to
 * any node in the MST). These can be processed with subsequent calls
 * to this routine.
 *
 * @param mst ON ENTRY may contain a pre-existing tree, but typically
 * is empty. ON EXIT contains the MST (pairs of linked nodes).
 */
void primsAlgo(std::map<size_t, Mojo>& available, std::vector< Idx2 >& mst);

/******************************************************************************/

double jaccardIndex(size_t sz1, size_t sz2, size_t iSz);
double jaccardDistance(size_t sz1, size_t sz2, size_t iSz);

template<typename Container>
double jaccardIndex(const Container& s1, const Container& s2){
  Container intersect;
  set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),
                  std::inserter(intersect,intersect.begin()));
  return jaccardIndex( s1.size(), s2.size(), intersect.size() );
} // jaccardIndex

template<typename Container>
double jaccardDistance(const Container& s1, const Container& s2){
  return -log2(jaccardIndex(s1,s2));
} // jaccardDistance

/**
 * Calculates a distance measure (Jaccard distance on RV ids) of how
 * far (on average) a set of factors are from each other. The default
 * percentiles (perc1, perc2) below are set up to return the distance
 * between a distant pair selected from all the close-by factor pairs.
 *
 * @param varsInFactor The key value of the map serves as a unique
 * index to the particular factor, this is ignored. The value field of
 * the map (Container) lists the RVs in that particular factor.
 *
 * @param perc1 For each factor we find its Jaccard distance to all
 * the others. perc1 is a percentile value to select one of these
 * distances to represent the typical distance of this factor to the
 * others. It determines whether we focus on the close-by neighbours,
 * or further removed neighbours. The low default value (0.01) chooses
 * one of the close-by factors. (Non-overlapping factor pairs are
 * common and will have an inf distance).
 *
 * @param perc2 After applying perc1 we have a representive distance
 * between pairs of factors. perc2 is a percentile point to select one
 * of these distances as the final distance value - it determines
 * whether we are interested in close pairs, the median distance
 * etc. The default value of 0.9 chooses one of the furthest of these
 * distances; i.e. a weak link in the chain of close-by factors.
 *
 * @return A value representing an averaged Jaccard distance between
 * the factors in varsInFactor.
 */
template<typename Container>
double betweenFactorDistortion(
  const std::map< size_t, Container >& varsInFactor,
  double perc1 = 0.01,
  double perc2 = 0.9){
  std::vector<double> distances;
  for (auto itr1 = varsInFactor.begin(); itr1 != varsInFactor.end(); itr1++) {

    std::vector<double> tempDist;
    for (auto itr2 = varsInFactor.begin(); itr2 != varsInFactor.end(); itr2++) {
      if (itr2 == itr1) {continue;}
      tempDist.push_back( jaccardDistance(itr1->second, itr2->second) );
    } // for
    std::nth_element(tempDist.begin(), tempDist.begin() + size_t( perc1*tempDist.size() ), tempDist.end());
    distances.push_back(tempDist[size_t( perc1*tempDist.size() )]);

  } // for

  std::nth_element(distances.begin(), distances.begin() + size_t( perc2*distances.size() ), distances.end());
  return(distances[size_t( perc2*distances.size() )]);
} // betweenFactorDistortion


void plotStructure(const std::map<size_t, Mojo>& accel);
void plotStructure(const std::vector<Idx2>& accel);

#endif // MOJO_HPP
