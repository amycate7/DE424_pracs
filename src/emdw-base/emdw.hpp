#ifndef EMDW_HPP
#define EMDW_HPP

/*******************************************************************************
 *                                                                             *
 * AUTHOR:  JA du Preez                                                        *
 * Created on : /02/2013                                                       *
 *                                                                             *
 * PURPOSE: Contains some top level typedefs that are used in all emdw code    *
 * COPYRIGHT: University of Stellenbosch, all rights reserved.                 *
 *                                                                             *
 *******************************************************************************/

// standard headers
#include <string>  // string
#include <iostream> // istream, ostream
#include <fstream>      // std::ofstream
#include <memory>
#include <vector>
#include <set>
#include <random>

#include "oddsandsods.hpp"
#include "prlite_stlvecs.hpp"
#include "anytype.hpp"

template <typename T>
using rcptr = std::shared_ptr<T>;

template <typename T>
using uniqptr = std::unique_ptr<T>;

namespace emdw {

extern bool debug;

  /**
   * RVIdType and RVIds respectively identifies a particular random
   * variable or sequence of random variables. Think of it as the
   * subscript when we talk of X_1, X_2 etc.
   *
   * In our processing we will operate on a set of Random Variables,
   * the members of which can be of different kinds. To identify them
   * each will be associated with an id. Without loss of generality (i
   * believe) we will use an 'unsigned' to identify each RV. Since we
   * are going to use several, we contain their ids in a vector.Eg, if
   * for instance the id vector is {7,9,11} it indicates RVs X7, X9
   * and X11.
   *
   * WHY NOT A SET? Alternatively we can collect the RV ids in a set,
   * after all they are sorted and unique. However, for ease of use we
   * allow the user to specify RV's in either an arbitrary order
   * (i.e. a vector), or a presorted version which is computationally
   * more efficient. If we use a set, it will internally (and outside
   * of our control) also sort the RV's. Externally we will have to
   * find a map between the arbitrary order RV's and the sorted
   * internal set, hopefully without sorting them again. It is not
   * immediately clear to me how to do this more efficiently than with
   * a second sort - therefore I am avoiding the issue by sticking to
   * a vector as internal representation. This has the additional
   * small benefit of making random access slightly more efficient.
   * !!!JADP: It might be worth while to investigate implementing a
   * custom set class (of course based on a redblack tree or similar)
   * with sufficient support for the operations we are interested in.
   */
  typedef unsigned RVIdType;
  typedef std::vector<RVIdType> RVIds;

  /**
   * Random Variable assignments:
   * Now for the values that the RV can take on. Because we will allow
   * multiple types such as bool, int, double, string etc, we will use
   * the AnyType container for this.
   */
  typedef std::vector<AnyType> RVVals;

  /**
   * Graphs: Factor graph, cluster graph, region graph, decision graph
   * log vs linear
   */

  /**
   * Inference algos: sumprod, maxsum, loopy, JT, parentchild
   */

  /**
   * Learning:
   */

  /**
   * Structure learning:
   */

  /// Here we declare random number initialisation constructs to use
  /// over the whole system so that we can do repeatable experiments.
  class RandomEngine{
  public:
    RandomEngine(unsigned theSeedVal = 0);
    std::default_random_engine& operator()();
    void setSeedVal(unsigned theSeedVal);
    unsigned getSeedVal() const;
  private:
    std::random_device randDev;
    unsigned seedVal;
    std::default_random_engine randEng;
  }; // RandomEngine

  extern emdw::RandomEngine randomEngine;

} // namespace

// used to contain pairs of indices to eg Factors
using Idx2 = std::pair<unsigned,unsigned>;
std::ostream& operator<<( std::ostream& file, Idx2 pr ) ;
std::istream& operator>>( std::istream& file, Idx2& pr );

template <typename T>
std::ostream& operator<<( std::ostream& file, std::set<T>& el);
template <typename T>
std::istream& operator>>( std::istream& file, std::set<T>& theSet );

#include "emdw.tcc"

#endif // EMDW_HPP
