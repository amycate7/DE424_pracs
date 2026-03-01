#ifndef DTW_HPP
#define DTW_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez
          COPYRIGHT: University of Stellenbosch, all rights reserved.
          DATE:    December 2020
          PURPOSE: DTW

*******************************************************************************/

// A categorical distance measure that can be used with dtw. Ref and
// Tst must support != operator.
template<typename Ref, typename Tst>
class HitMissDistance{
public:
  double operator()(const Ref& ref, const Tst& tst ) const{
    return (ref != tst);
  } //operator()
}; // HitMissDistance

// utility function to get minimum of three values
double min3(double a, double b, double c) {
  return a < b ? std::min<double>(a,c) : std::min<double>(b,c);
} // min3


/// @brief DTW match continuous sequences. Ref and Tst must both
/// support indexing via [], as well as .size(). Measure is a functor
/// for determining distances (such as euclidean).
///
/// @param ref, tst: Two sequences to be aligned with each other.
///
/// @measure: A functor (or lambda) that can calculate the distance
/// between ref[i] and tst[i].
///
/// @param subset: If true we match tst to only the best fitting
/// subsection of ref.

template<typename Ref, typename Tst, typename Measure>
double dtw(const Ref& ref, const Tst& tst, const Measure& measure, bool subset = false) {

  size_t szRef = ref.size(), szTst = tst.size();
  double dist = szRef+szTst;
  std::vector<double> buff1(szTst); std::vector<double>* ptrCurr = &buff1;
  std::vector<double> buff2(szTst); std::vector<double>* ptrPrev = &buff2;

  (*ptrCurr)[0] = measure(ref[0], tst[0]);
  for (size_t iTst = 1; iTst < szTst; iTst++) {
    (*ptrCurr)[iTst] = (*ptrCurr)[iTst-1] + measure(ref[0], tst[iTst]);
  } // for

  for (size_t iRef = 1; iRef < szRef; iRef++) {
    swap(ptrCurr, ptrPrev);
    if (subset) {
      (*ptrCurr)[0] = measure(ref[iRef], tst[0]);
    } // if
    else {
      (*ptrCurr)[0] = (*ptrPrev)[0] + measure(ref[iRef], tst[0]);
    } // else

    for (size_t iTst = 1; iTst < szTst; iTst++) {
      (*ptrCurr)[iTst] =
        min3(
          (*ptrPrev)[iTst],     // deletion
          (*ptrPrev)[iTst-1],   // substitution or correct
          (*ptrCurr)[iTst-1] )  // insertion
        + measure(ref[iRef], tst[iTst]);
    } // for
    if ( (*ptrCurr)[szTst-1] < dist) {dist = (*ptrCurr)[szTst-1];} // if
  } // for iRef

  return subset ? dist : (*ptrCurr)[szTst-1];
} // dtw

/// @brief SID (Substitutions - Deletions - Insertions) DTW. Ref and
/// Tst must both support indexing via [], as well as .size(). Measure
/// is a functor for determining distances. Differs subtly from
/// standard DTW in how insertions/deletions are handled with
/// consecutive repeating symbols.
///
/// @param ref, tst: Two sequences to be aligned with each other.
///
/// @measure: A functor (or lambda) that can calculate the distance
/// between ref[i] and tst[i].
///
/// @param subset: If true we match tst to only the best fitting
/// subsection of ref.

template<typename Ref, typename Tst, typename Measure>
double sid(const Ref& ref, const Tst& tst, const Measure& measure, bool subset = false) {

  size_t szRef = ref.size(), szTst = tst.size();
  double dist = szRef+szTst;
  std::vector<double> buff1(szTst); std::vector<double>* ptrCurr = &buff1;
  std::vector<double> buff2(szTst); std::vector<double>* ptrPrev = &buff2;

  (*ptrCurr)[0] = measure(ref[0], tst[0]);
  for (size_t iTst = 1; iTst < szTst; iTst++) {
    (*ptrCurr)[iTst] = (*ptrCurr)[iTst-1]+1.0;
  } // for
  //for (auto el : *ptrCurr) {cout << el << " ";} cout << endl;

  for (size_t iRef = 1; iRef < szRef; iRef++) {
    swap(ptrCurr, ptrPrev);
    if (subset) {
      (*ptrCurr)[0] = measure(ref[iRef], tst[0]);
    } // if
    else {
      (*ptrCurr)[0] = (*ptrPrev)[0]+1;
    } // else

    for (size_t iTst = 1; iTst < szTst; iTst++) {
      (*ptrCurr)[iTst] =
        min3(
          (*ptrPrev)[iTst]+1.0,                                         // deletion
          (*ptrPrev)[iTst-1]+ measure(ref[iRef], tst[iTst]),   // substitution or correct
          (*ptrCurr)[iTst-1]+1.0 );                                     // insertion
    } // for
    if ( (*ptrCurr)[szTst-1] < dist) {dist = (*ptrCurr)[szTst-1];} // if
    //for (auto el : *ptrCurr) {cout << el << " ";} cout << endl;
  } // for iRef

  return subset ? dist : (*ptrCurr)[szTst-1];
} // dtw

#endif // DTW_HPP
