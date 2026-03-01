#ifndef VECSET_HPP
#define VECSET_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez (Copyright University of Stellenbosch, all rights reserved.)
          DATE:    March 2013
          PURPOSE: std::vector union and intersection routines

*******************************************************************************/

#include <vector>

/******************************************************************************/

/// Does sorted vector vecL contain sorted vector vecR as a subset? light weight and fast
template<typename T>
bool contains(const std::vector<T>& vecL,
              const std::vector<T>& vecR);

/**
 * Simultaneously finds the union/intersection/partitions of two (sorted) vectors.
 *
 * The cross-referencing indices used in the following code, with example values:
 *
 * vecL = [2,4,5,6]      left side vec, must be pre-sorted
 * vecR = [1,2,3,5]      right side vec, must be pre-sorted
 *
 * vecU = [1,2,3,4,5,6]  the union of vecL and vecR
 * vecI = [2,5]          the intersection of vecL and vecR
 *
 * l2i  = [0,2]          the indices of vecL that are in the intersection
 * lni  = [1,3]          the indices of unique vals in L (i.e. not in the intersection)
 * r2i  = [1,3]          the indices of vecR that are in the intersection
 * rni  = [0,2]          the indices of unique vals in R (i.e. not in the intersection)
 *
 * i2l =  [1,4]          where in the union is l2i[k] located
 * ni2l = [3,5]          where in the union is lni[k] located
 * ni2r = [0,2]          where in the union is rni[k] located
 *
 * Note:
 * vecI = vecL[l2i]
 * vecU[i2l]  = vecL[l2i]
 * vecU[ni2l] = vecL[lni]
 * vecU[ni2r] = vecR[rni]
 *
 *                      vecU = vecL \/ vecR
 *         ______________/\____________
 *        /  ________________           \
 *          /          ______\_________
 *         |  vecL    /.......|        \
 *         |         |..vecI..|  vecR   |
 *          /________|/....../          |
 *         /         /\____________/___/
 *        /         /             /
 *  _____/_________/_____________/__________
 *  | vecL-vecR  | vecL/\vecR | vecR-vecL  |
 *  +------------+------------+------------+
 *  | vecL[lni ] |  vecL[l2i] |    ---     |
 *  |    ---     |  vecR[r2i] | vecR[rni ] |
 *  | vecU[ni2l] |  vecU[i2l] | vecU[ni2r] |
 *  +------------+------------+------------+
 *
 *
 *  Test code (for the above example):

  vector<int> vecL = {2,4,5,6};
  vector<int> vecR = {1,2,3,5};
  vector<size_t> l2i,lni,r2i,rni,i2l,ni2l,ni2r;
  vector<int> vecU =
    cupcap(vecL,vecR,l2i,lni,r2i,rni,i2l,ni2l,ni2r);

  cout << endl << vecL << endl << vecR << endl;
  cout << endl << vecU << endl << extract<int>(vecL, l2i) << endl;
  cout << endl << l2i << endl << lni << endl;
  cout << endl << r2i << endl << rni << endl;
  cout << endl << i2l << endl << ni2l << endl << ni2r << endl;

*/
template<typename T, typename TI>
std::vector<T> cupcap(const std::vector<T>& vecL,
                      const std::vector<T>& vecR,
                      std::vector<TI>& l2i,
                      std::vector<TI>& lni,
                      std::vector<TI>& r2i,
                      std::vector<TI>& rni,
                      std::vector<TI>& i2l,   // same as i2r, so we skip that
                      std::vector<TI>& ni2l,
                      std::vector<TI>& ni2r);

/**
 * Finds the union of two (sorted) integer vectors.
 *
 * The cross-referencing indices used in the following code, with example values:
 *
 * vecL = [2,4,5,6]      left side vec, must be pre-sorted
 * vecR = [1,2,3,5]      right side vec, must be pre-sorted
 * vecU = [1,2,3,4,5,6]  the union of vecL and vecR
 *
 * u2l  = [1,3,4,5]      where in the union does L go to
 * u2r  = [0,1,2,4]      where in the union does R go to
 *
 * Note: union indices 1 and 4 are common to both L and R
 *
 *                   vecU = vecL\/vecR
 *         ______________/\____________
 *        /  _______________            \
 *          /            ___\___________
 *         |   vecL=    /    |          \
 *         | vecU[u2l] |     |  vecR=    |
 *          \__________|____/ vecU[u2r]  |
 *                      \_______________/
 *
 *
 *  Test code:

  vector<int> vecL = {2,4,5,6};
  vector<int> vecR = {1,2,3,5};
  vector<int> u2l;
  vector<int> u2r;
  vector<int> vecU = sortedUnion(vecL, vecR, u2l, u2r);
  cout << "\nvecL: " << vecL
       << "\nvecR: " << vecR
       << "\nvecU: " << vecU
       << "\nu2l: " << u2l
       << "\nu2r: " << u2r
       << endl;
*/

template<typename T, typename TI>
std::vector<T> sortedUnion(const std::vector<T>& vecL,
                           const std::vector<T>& vecR,
                           std::vector<TI>& u2l,
                           std::vector<TI>& u2r);

/**
 * Finds the union of two (sorted) integer vectors. The older version
 * of the code where r2u contains only the elements unique to
 * R. Uncertain whether there are solid reasons for not rather using
 * the other version.
 *
 * The cross-referencing indices used in the following code, with example values:
 *
 * vecL = [2,4,5,6]      left side vec, must be pre-sorted
 * vecR = [1,2,3,5]      right side vec, must be pre-sorted
 * vecU = [1,2,3,4,5,6]  the union of vecL and vecR
 *
 * l2u  = [0,1,2,3]      which indices of L will be used (should be all)
 * u2l  = [1,3,4,5]      to where in the union will L[k] go
 * r2u = [0,2]           which *additional* positions of the R will be used,
 *                       i.e. the indices that are unique to R
 * u2r  = [0,2]          to where in the union will those *additional* unique
 *                       indices of R, i.e. R[r2u[k]], go.
 *
 * Note:
 * vecU[u2l] = vecL[l2u]
 * vecU[u2r] = vecR[r2u]
 *
 *                   vecU = vecL\/vecR
 *         ______________/\____________
 *        /  _______________            \
 *          /            ___\___________
 *         |   vecL=    /    |          \
 *         | vecU[u2l] |     |  vecR=    |
 *          \__________|____/ vecU[u2r]  |
 *                      \_______________/
 *
 *
 *  Test code:

 vector<int> vecL = {2,4,5,6};
 vector<int> vecR = {1,2,3,5};


 vector<int> l2u;
 vector<int> u2l;
 vector<int> r2u;
 vector<int> u2r;
 vector<int> vecU = sortedUnion(vecL, vecR, l2u, u2l, r2u, u2r);
 cout << "\nvecL" << vecL
 << "\nvecR" << vecR
 << "\nvecU" << vecU
 << "\nl2u" << l2u
 << "\nu2l" << u2l
 << "\nr2u" << r2u
 << "\nu2r" << u2r
 << endl;

*/

template<typename T, typename TI>
std::vector<T> sortedUnion(const std::vector<T>& vecL,
                           const std::vector<T>& vecR,
                           std::vector<TI>& l2u,
                           std::vector<TI>& u2l,
                           std::vector<TI>& r2u,
                           std::vector<TI>& u2r);

/******************************************************************************/

/**
 * Returns the intersection of two (sorted) integer vectors, also
 * returning by reference the indices of the unique/difference vals.
 *
 * The cross-referencing indices used in the following code, with example values:
 *
 * vecL = [2,4,5,6]  left side vec, must be pre-sorted
 * vecR = [1,2,3,5]  right side vec, must be pre-sorted
 * vecI = [2,5]      the intersection of vecL and vecR
 *
 * l2i  = [0,2]      the indices in L containing the intersection
 * r2i  = [1,3]      the indices in R containing the intersection
 * lni  = [1,3]      the indices of unique vals in L (not intersection)
 * rni  = [0,2]      the indices of unique vals in R (not intersection)
 *
 * Note:
 * vecI = vecL[l2i] = vecR[r2i]
 * l2i union lni = 0:sz(vecL)-1
 * r2i union rni = 0:sz(vecR)-1
 *
 *                      vecU = vecL\/vecR
 *            ______________/\____________
 *           /  ________________           \
 *             /          ______\_________
 *            |  vecL    /.......|        \
 *            |         |..vecI..|  vecR   |
 *             /________|../..../          |
 *            /          \/___________/___/
 *           /           /           /
 *          /     vecI=vecL/\vecR   /
 *  _______/__________ /___________/________
 *  | vecL-vecR: |    vecI:   | vecR-vecL: |
 *  +------------+------------+------------+
 *  | vecL[lni ] |  vecL[l2i] |    ---     |
 *  |    ---     |  vecR[r2i] | vecR[rni ] |
 *  +------------+------------+------------+
 *
 *  Test code:

 vector<int> vecL = {2,4,5,6};
 vector<int> vecR = {1,2,3,5};

 vector<int> l2i;
 vector<int> r2i;
 vector<int> vecI = sortedIntersection(vecL, vecR, l2i, r2i);
 cout << "\nvecL" << vecL
 << "\nvecR" << vecR
 << "\nvecI" << vecI
 << "\nl2i" << l2i
 << "\nr2i" << r2i
 << endl;

*/

template<typename T, typename TI>
std::vector<T> sortedIntersection(const std::vector<T>& vecL,
                                  const std::vector<T>& vecR,
                                  std::vector<TI>& l2i,
                                  std::vector<TI>& r2i,
                                  std::vector<TI>& lni,
                                  std::vector<TI>& rni);

// slightly faster version without the unique/difference indices
template<typename T, typename TI>
std::vector<T> sortedIntersection(const std::vector<T>& vecL,
                                  const std::vector<T>& vecR,
                                  std::vector<TI>& l2i,
                                  std::vector<TI>& r2i);

/******************************************************************************/

/**
 * Insert an element into the correct position of a presumingly sorted set
 * of values if no duplicate exists (otherwise do nothing).
 *
 * vecL: A vector containing a set of values
 * elem: A value to be added to the set of vecL
 *
 *  Test code:

 vector<int> vecL = {};
 for(const int &i: {2,4,5,6,2,3,5}){
    sortedSetInsert(vecL, i);
 }
 cout << "vecL: " << vecL << endl;
 //5  2 3 4 5 6
*/

template<typename T, typename TI>
void sortedSetInsert(std::vector<T> &vecL,
                     const TI &elem);


/******************************************************************************
 * Create a vector from the values passed to this function.
 */

template <class T, class... Ts>
std::vector<T> make_vector(T&& first, Ts&&... args);


#include "vecset.tcc"

#endif // VECSET_HPP
