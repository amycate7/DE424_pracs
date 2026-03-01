#ifndef CAT1EXPF_HPP
#define CAT1EXPF_HPP

/*******************************************************************************

AUTHOR: JA du Preez (Copyright University of Stellenbosch, all rights reserved.)

DATE: September 2020

PURPOSE: This is a sort of nursary to hatch categorical classes
that can inter-operate with both BP and VMP. These will at a later
stage be moved to fully fledged emdw factors and operators.

For BP see: D. Koller and N. Friedman, Probabilistic Graphical Models:
Principles and Techniques, MIT Press, 2009.

For VMP see: John Winn and Christopher M. Bishop, Variational message
passing, JMLR-2005, vol 6, p 661-694.

A note about VMP:

1. A VMP pdf always remains in conditional form. Multiplying in the
message from the priors simply provide the conditional parameters (via
the expected moments of the prior). In contrast, multiplying in a BP
message changes the distribution to a joint, thereby increasing its
dimension.

2. In standard VMP a message from the prior/parent to the child
undergoes no message cancelation (i.e. at any given stage all messages
from the prior are identical), although in our experience it is indeed
helpful to cancel the previous message from the parent.

3. Marginalisation towards a prior (i.e. child to parent) is done by
re-arranging the exponential form to match that of the prior (easy to
do because of the multi-linear form). The resultant natural parameters
defines the pdf that is the marginalized message destined for the
prior. A very nice trick to avoid intractable integrals and
non-conjugate forms.

4. This form implicitly does message cancellation and therefore no
further message cancellation is involved (or required). However, if
one wants the system to, from the outside, conform to standard BP
rules, one can add the extra offset that message cancelation will then
just reverse.

*******************************************************************************/

#include <map>
#include "emdw.hpp"

// ##################
// ### Cat1EXPF ###
// ##################

/**
 * A very basic unconditional scalar exponential form Categorical. It
 * is used for messages, and also as a building block for Cat1VMP.
 *
 * NOTE: The 1D version of DiscreteTable.
 */

// forward declarations
template <typename AssignType> class Cat1VMP;
template <typename AssignType> class DiscreteTable;

/**
 * A scalar canonical form Categorical.
 *
 * TODO: Currently we do not track the normalization term. This makes
 * it implicitly normalized.
 *
 * NOTE: The 1D version of DiscreteTable.
 */
template <typename AssignType>
class Cat1EXPF {

  friend class Cat1VMP<AssignType>;
  friend class DiscreteTable<AssignType>;

public:
  /// Although we save linear probs, the corresponding natural params
  /// are log probs.
  ///
  /// NOTE: In the vacuous case we take all probabilities to be 1.0
  /// (i.e. it is unnormalized). This makes all natural parameters to
  /// be zero, which is a nice initial for vmp.
  Cat1EXPF( emdw::RVIdType theVarId,
            const rcptr< std::vector<AssignType> >& theDomain,
            double theDefProb = 1.0,
            const std::map<AssignType, double>& theSparseProbs = std::map<AssignType, double>() )
      : varId(theVarId),
        domain(theDomain),
        defProb(theDefProb),
        sparseProbs(theSparseProbs) {} // ctor

  Cat1EXPF* copy () const {
    return new Cat1EXPF(*this);
  } // copy

  Cat1EXPF* vacuousCopy () {
    return new Cat1EXPF(varId, domain);
  } // vacuousCopy

  ///
  void inplaceMpy(const rcptr< Cat1EXPF<AssignType> >& rhs) {
    std::map<AssignType, double>& mapL(sparseProbs);         // convenient alias
    double& defL(defProb);                                 // convenient alias

    const std::map<AssignType, double>& mapR(rhs->sparseProbs); // convenient alias
    const double& defR(rhs->defProb);                         // convenient alias

    // Ensures the number of parameters of the two categorical potential functions are equal
    PRLITE_ASSERT(domain->size() == rhs->domain->size(),
                  "Cat1EXPF::inplaceMpy : lhs size " << domain->size() <<
                  " not equal to rhs size " << rhs->domain->size() <<
                  "\nlhs var " << varId << " rhs var " << rhs->varId );

    // Case 1: with a uniform rhs we can simply return the lhs.
    if ( !mapR.size() ) {
      //std::cout << "Case 1/\n";

      // one can also leave the following block out altogether, in that
      // case the result may be unnormalized
      defL *= defR;
      auto iL = mapL.begin();
      while ( iL != mapL.end() ){
        iL->second *= defR;
        iL++;
      } // while

      return;
    } // if

    // Case 2: efficient handling of relatively sparse rhs with defProb == 0)
    if (defR == 0.0 and mapR.size() <= mapL.size()/10){
      for (auto r : mapR) {
        auto itr = mapL.find(r.first);
        if ( itr != end(mapL) ){
          itr->second *= r.second;
        } // if
        else {
          mapL.insert({r.first, defL*r.second});
        } // else
      } // for r
      return;
    } // if

    // Case 3: For all others, do the double ladder

    // to calc the number of defaults occurring in neither of the maps
    // to calc the number of defaults occurring in neither of the maps
    size_t nDefsL = domain->size()-mapL.size();
    size_t nDefsR = rhs->domain->size() - mapR.size();
    double newDef = defL*defR;

    auto iL = mapL.begin();
    auto iR = mapR.begin();
    while ( iL != mapL.end() && iR != mapR.end() ){

      if (iL->first < iR->first) {      // only in mapL
        nDefsR--;
        iL->second *= defR;
        iL++;
      } // if

      else if (iL->first > iR->first) { // only in mapR
        nDefsL--;
        double prob = defL * iR->second;
        if (std::abs(prob - newDef) > 0.0) {
          mapL.emplace_hint(iL, iR->first, prob);
        } // if
        iR++;
      } // else if

      else {                            //common to both maps
        iL->second *= iR->second;
        iL++; iR++;
      } // else

    } // while

    while ( iL != mapL.end() ){         // some extra spec's in mapL
      nDefsR--;
      iL->second *= defR;
      iL++;
    } // while

    while ( iR != mapR.end() ){         // some extra spec's in mapR
      nDefsL--;
      double prob = defL * iR->second;
      if (std::abs(prob - newDef) > 0.0) {
        mapL.emplace_hint(iL, iR->first, prob);
      } // if
      iR++;
    } // while

    PRLITE_ASSERT(
      nDefsL == nDefsR,
      "Cat1EXPF::inplaceMpy : The number of cases in neither map should match, no matter from which side you look at it.");

    // there are some unused default values
    defL = newDef;
    return;

  } // inplaceMpy

  ///
  Cat1EXPF* mpy(const rcptr<Cat1EXPF>& other) const {
    Cat1EXPF* answ = new Cat1EXPF(*this);
    answ->inplaceMpy(other);
    return answ;
  } // mpy

  ///
  void inplaceDiv(const rcptr< Cat1EXPF<AssignType> >& rhs) {
    std::map<AssignType, double>& mapL(sparseProbs);         // convenient alias
    double& defL(defProb);                                 // convenient alias

    const std::map<AssignType, double>& mapR(rhs->sparseProbs); // convenient alias
    const double& defR(rhs->defProb);                         // convenient alias

    // Ensures the number of parameters of the two categorical potential functions are equal
    PRLITE_ASSERT(domain->size() == rhs->domain->size(),
                  "Cat1EXPF::inplaceDiv : lhs size " << domain->size() <<
                  " not equal to rhs size " << rhs->domain->size() <<
                  "\nlhs var " << varId << " rhs var " << rhs->varId );

    // Case 1: with a uniform rhs we can simply return the lhs.
    if ( !mapR.size() ) {
      //std::cout << "Case 1/\n";

      // one can also leave the following block out altogether, in that
      // case the result may be unnormalized
      defL /= defR;
      auto iL = mapL.begin();
      while ( iL != mapL.end() ){
        iL->second /= defR;
        iL++;
      } // while

      return;
    } // if

    // Case 2: efficient handling of relatively sparse rhs with defProb == 0)
    if (defR == 0.0 and mapR.size() <= mapL.size()/10){
      for (auto r : mapR) {
        auto itr = mapL.find(r.first);
        if ( itr != end(mapL) ){
          itr->second /= r.second;
        } // if
        else {
          mapL.insert({r.first, defL/r.second});
        } // else
      } // for r
      return;
    } // if

    // Case 3: For all others, do the double ladder

    // to calc the number of defaults occurring in neither of the maps
    size_t nDefsL = domain->size()-mapL.size();
    size_t nDefsR = rhs->domain->size() - mapR.size();
    double newDef = defR > 0.0 ? defL/defR : defL;

    auto iL = mapL.begin();
    auto iR = mapR.begin();
    while ( iL != mapL.end() && iR != mapR.end() ){

      if (iL->first < iR->first) {      // only in mapL
        nDefsR--;
        if(defR > 0) {
          iL->second /= defR;
        } // if
        iL++;
      } // if

      else if (iL->first > iR->first) { // only in mapR
        nDefsL--;
        double prob = iR->second > 0.0 ? defL / iR->second : 0.0;
        if (std::abs(prob - newDef) > 0.0) {
          mapL.emplace_hint(iL, iR->first, prob);
        } // if
        iR++;
      } // else if

      else {                            //common to both maps
        if (iR->second > 0.0) {
          iL->second /= iR->second;
        } // if
        iL++; iR++;
      } // else

    } // while

    while ( iL != mapL.end() ){         // some extra spec's in mapL
      nDefsR--;
      if (defR > 0.0) {
        iL->second /= defR;
      } // if
      iL++;
    } // while

    while ( iR != mapR.end() ){         // some extra spec's in mapR
      nDefsL--;
      double prob = iR->second > 0.0 ? defL / iR->second : 0.0;
      if (std::abs(prob - newDef) > 0.0) {
        mapL.emplace_hint(iL, iR->first, prob);
      } // if
      iR++;
    } // while

    PRLITE_ASSERT(
      nDefsL == nDefsR,
      "Cat1EXPF::inplaceDiv : The number of cases in neither map should match, no matter from which side you look at it.");

    // there are some unused default values
    defL = newDef;
    return;

  } // inplaceDiv

  ///
  Cat1EXPF* div(const rcptr<Cat1EXPF>& other) const {
    Cat1EXPF* answ = new Cat1EXPF(*this);
    answ->inplaceDiv(other);
    return answ;
  } // div

  void inplaceNormalize() {

    size_t cnt = 0;
    double sum = 0.0;
    for (auto& el: sparseProbs) {
      cnt++;
      sum += el.second;
    } // for
    sum += (domain->size()-cnt)*defProb;

    defProb /= sum;
    for (auto& el: sparseProbs) {
      el.second /= sum;
    } // for

  } // inplaceNormalize

  ///
  Cat1EXPF* normalize() const {
    Cat1EXPF* answ = new Cat1EXPF(*this);
    answ->inplaceNormalize();
    return answ;
  } // normalize

  /**
   * A more direct version of potentialAt
   *
   * @param theAss The specific assignment to the random variables for
   *   which the probability must be returned. It is assumed that the
   *   corresponding random variables are sorted according to
   *   ascending RVId.
   *
   * @return The (unnormalised) probabilty at the RV assignment given
   *   by theAss.
   */
  double getProb(const AssignType& theAss) const{
    auto itr = sparseProbs.find(theAss);
    if ( itr == sparseProbs.end() ) {return defProb;} // if
    else {return itr->second;} // else
  } // getProb

  emdw::RVIdType getVarId() const {return varId;} // getVarId

  /// Output to file. Currently mostly for information
  std::ostream& txtWrite(std::ostream& file) const {
    file << "Cat1EXPF: "
         << std::endl;
    file << "defProb: " << defProb << std::endl;
    file << "sparseProbs: " << sparseProbs.size() << std::endl;
    for (auto el : sparseProbs) {
      file << el.first << " " << el.second << "  ";
    } // for
    file << std::endl;
    return file;
  } // txtWrite

  friend std::ostream& operator<<( std::ostream& file,
                                   const Cat1EXPF& pdf ) {
    pdf.txtWrite(file);
    return file;
  } // operator<<

private:

  // *******************************
  // *** definitive data members ***
  // *******************************

  emdw::RVIdType varId;

  // domain of catVar
  rcptr< std::vector<AssignType> > domain;

  // natural parameters are log(prob). However, how we represent them
  // internally is up to us. We opt for linear here, it makes
  // interaction with Dirichlet variables much cheaper.
  double defProb;
  std::map<AssignType, double> sparseProbs;

  // ******************************
  // *** secondary data members ***
  // ******************************


}; // Cat1EXPF


#endif // CAT1EXPF_HPP
