#ifndef MVG_HPP
#define MVG_HPP

/*******************************************************************************
 *                                                                             *
 * AUTHOR:  JA du Preez                                                        *
 * Created on : 19/09/2014                                                     *
 *                                                                             *
 * PURPOSE: The abstract base class for all general discrete factors           *
 * COPYRIGHT: University of Stellenbosch, all rights reserved.                 *
 *                                                                             *
 *******************************************************************************/

// standard headers
#include <string>  // string
#include <iosfwd>  // istream, ostream

#include "emdw.hpp"
#include "factor.hpp"

/**
 * An abstract base class to group all general gaussian factors.
 */
class MVG : public Factor {

  //======================== Constructors and Destructor ========================
  public:

  virtual ~MVG() {}

  //========================= Operators and Conversions =========================
public:

  //=========================== Inplace Configuration ===========================
public:
  ///
  virtual std::istream& txtRead(std::istream& file) = 0;
  ///
  virtual std::ostream& txtWrite(std::ostream& file) const = 0;

  //===================== Required Virtual Member Functions =====================
public:

  /**
   * With no params it simply returns a copy of the existing Factor,
   * with params it returns an identical Factor, but with variables
   * replaced by another set.
   *
   * @param newVars These variables will replace, one-for-one AND IN
   * THIS ORDER, the original SORTED set.
   *
   * @param presorted Set this to true if theVars are sorted from
   * smallest to biggest and avoid sorting and repacking costs.
   */
  virtual MVG* copy(const emdw::RVIds& newVars = {},
                       bool presorted = false) const = 0;

  /// note rhs is a raw pointer
  virtual bool isEqual(const Factor* rhsPtr) const = 0;

  /// For the while being, use the default implementation of this member
  double distanceFromVacuous() const = 0; //{return Factor::distanceFromVacuous();}

  ///
  virtual unsigned noOfVars() const = 0;

  /// returns the ids of the variables the factor is defined on, in
  /// *sorted* order
  virtual emdw::RVIds getVars() const = 0;

  /// in *sorted* order
  virtual emdw::RVIdType getVar(unsigned varNo) const = 0;

  /**
   * returns the factor potential at a given variable
   * assignment. *NOTE* the potential is supplied directly as is, no
   * normalization etc is done.
   *
   * @param theVars The RV ids in the order that the values in
   * theirVals assumes them to be.
   *
   * @param theirVals The values assigned to this random vector.
   *
   * @param presorted Set this to true if theVars are sorted from
   * smallest to biggest and avoid sorting and repacking costs.
   */
  virtual double potentialAt(const emdw::RVIds& theVars,
                             const emdw::RVVals& theirVals,
                             bool presorted = false) const = 0;

  //====================== Other Virtual Member Functions =======================
public:

  // /**
  //  * make a copy of this factor, but set the potentials to the default
  //  * unity factor values. The client has to manage the lifetime of
  //  * this pointer - putting it in a rcptr will be a good idea.
  //  */
  // virtual MVG* vacuousCopy() const = 0;

  // /**
  //  * measures the distance between two factors. zero if equal, else
  //  * positive and growing with degree of difference.
  //  */
  // virtual double distance(const Factor* rhsPtr) const = 0;

  //======================= Non-virtual Member Functions ========================
public:

}; // MVG


#endif // MVG_HPP
