/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    October 2014
          PURPOSE: Information theory routines

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/
#include "anytype.hpp"
#include "prlite_stlvecs.hpp"

template <typename T>
double entropy(
  const std::map< std::vector<T>, double>& p) {

  double sumP = 0.0;
  for (auto elem : p) {
    sumP += elem.second;
  } // for

  double h = 0;
  for (auto& elem : p) {
    double pn = elem.second/sumP;
    if (pn > 0) { // p log(p) == 0 with p == 0
      h -= pn*log2(pn);
    } // if
  } // for

  return h;
} // entropy

/**
 * @brief returns the possible value assignments that each
 * dimension/column of the data matrix can take on. Typically a helper
 * function required for Dirichlet smoothing.
 *
 * @param data is a matrix containing in its rows a sequence of (equal
 * length) data-vectors.
 *
 * @return An array of vectors where the i'th vector returns the
 * sorted sequence of values found in the i'th column of the 'data'
 * input matrix.
 */
template <typename T>
std::vector< std::vector<T> > getDomains(const std::vector< std::vector<T> >& data) {

  unsigned nVecs = data.size();
  unsigned nDim = data[0].size();
  std::vector< std::set<unsigned> > domainsets(nDim);

  for (unsigned r = 0; r < nVecs; r++) {
    for (unsigned c = 0; c < nDim; c++) {
      domainsets[c].insert(data[r][c]);
    } // for
  } // for

  std::vector< std::vector<T> > domainvecs(nDim);
  for (unsigned n = 0; n < nDim; n++) {
    domainvecs[n].insert( domainvecs[n].end(), domainsets[n].begin(), domainsets[n].end() );
  } // for

  return domainvecs;
} // getDomains

/**
 * @brief Calculates the mutual information between a set of RVs X and another set of RVs y
 *
 * @param datInt the full discrete set of observed values, including also somewhere the target values
 * @param domains the range of values each RV can take on
 * @param varsX the indices into datInt containing the X RVs
 * @param varsY the indices into datInt containing the y RVs
 * @param entY returns by reference the estimated joint entropy for the y RVs
 * @param smoothing The type of Dirichlet smoothing to apply 0: none, 1: partial, 2: full
 */
template<typename T>
double estMI(
  const std::vector< std::vector<T> >& datInt,
  std::vector< std::vector<T> > domains,
  const std::vector<unsigned>& varsX,
  const std::vector<unsigned>& varsY,
  double& entY,
  unsigned smoothing) {

  unsigned nVecs = datInt.size();
  //unsigned nVars = datInt[0].size();
  std::vector<unsigned> varsXY = varsX; varsXY.insert( varsXY.end(), varsY.begin(), varsY.end() );
  std::vector< std::vector<T> > datSel;

  // first extract the relevant subset of rvs
  for (unsigned n = 0; n < nVecs; n++) {
    datSel.push_back(extract<T>(datInt[n], varsXY));
  } // for n
  //cout << datSel << endl; exit(-1);


  // set the indexing into datSel
  std::vector<unsigned> idxX( varsX.size() );
  std::vector<unsigned> idxY( varsY.size() );
  std::vector<unsigned> idxXY( varsXY.size() );
  prlite::setVals(idxX,0u,1u);
  prlite::setVals(idxY,static_cast<unsigned>(idxX.size()),1u);
  prlite::setVals(idxXY,0u,1u);

  // the counts go in here
  std::map<std::vector<T>,double> pX;
  std::map<std::vector<T>,double> pY;
  std::map<std::vector<T>,double> pXY;

  // and count them
  for (unsigned n = 0; n < datSel.size(); n++) {
    pX[extract<T>(datSel[n], idxX)]++;
    pY[extract<T>(datSel[n], idxY)]++;
    pXY[extract<T>(datSel[n], idxXY)]++;
  } // for n
  //for (auto cnt : pXY) {cout << cnt.first << ": " << cnt.second << endl;}  cout << endl; exit(-1);
  datSel.clear();

  if (smoothing==1) {
  //*

  // this is a partial dirichlet phantom set, we only extend the
  // feature std::vectors that were actually observed.

    std::vector< std::vector<T> > outcomesY =
      getAllCombinations( extract< std::vector<T> >(domains, varsY) );

    for (auto& outcomeX : pX) {
      //cout << outcomeX.first << endl;
      for (auto& outcomeY : outcomesY) {
        std::vector<T> outcome = outcomeX.first;
        outcome.insert( outcome.end(), outcomeY.begin(), outcomeY.end() );
        datSel.push_back(outcome);
        //cout << datSel[datSel.size()-1] << endl;
      } // for
    } // for

    // and count them
    //double weight = 1.0;
    double weight = 1.0/pow(datSel.size(), 0.25);
    // cout << "Adding another " << datSel.size()
    //      << " vectors, each with weight " << weight << endl;
    for (unsigned n = 0; n < datSel.size(); n++) {
      pX[extract<T>(datSel[n], idxX)]+=weight;
      pY[extract<T>(datSel[n], idxY)]+=weight;
      pXY[extract<T>(datSel[n], idxXY)]+=weight;
    } // for n

    // */
  } // if

  if (smoothing==2) {

    //*
    // this is a full dirichlet phantom observation set, we extend all
    // possible feature outcomes. Because there are so many we reduce
    // the weight of each outcome

    datSel = getAllCombinations( extract< std::vector<T> >(domains, varsXY) );
    double weight = 1.0; // /sqrt(datSel.size());
    //cout << varsXY << endl;
    //cout << "Adding another " << datSel.size()
    //      << " vectors, each with weight " << weight << endl;

    // and count them
    for (unsigned n = 0; n < datSel.size(); n++) {
      pX[extract<T>(datSel[n], idxX)]+=weight;
      pY[extract<T>(datSel[n], idxY)]+=weight;
      pXY[extract<T>(datSel[n], idxXY)]+=weight;
    } // for n

    // */
  } // if

  //for (auto elem : pX) {cout << "X: " << elem.first << " " << elem.second << endl;}
  //for (auto elem : pY) {cout << "Y: " << elem.first << " " << elem.second << endl;}
  //for (auto elem : pXY) {cout << "XY: " << elem.first << " " << elem.second << endl;}
  // cout << "H(X):    " << entropy(pX)
  //      << ", H(Y):   " << entropy(pY)
  //      << ", H(X,Y): " << entropy(pXY) << endl;

  // calculate the entropies and MI
  entY = entropy(pY);
  return entropy(pX) +  entY - entropy(pXY);
} // estMI
