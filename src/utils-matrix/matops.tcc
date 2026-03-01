/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    December 2013
          PURPOSE: some useful matrix routines
          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/


// lower triangle of m -= v*v' with v a column vector
template<typename TOut, typename TInVec>
void
subtractOuterProdSym(
  prlite::RowMatrix<TOut>& m,
  const TInVec& v){
  auto vPtr = v.begin();
  auto mrPtr = m.rowBegin();
  while ( vPtr != v.end() ) {
    auto wPtr = v.begin();
    auto v2End = vPtr+1;
    auto mPtr = (*mrPtr).begin();
    while (wPtr != v2End) {
      *mPtr++ -= *vPtr * *wPtr++;
    } // while
    vPtr++;
    mrPtr++;
  } // while
} // subtractOuterProd

// lower triangle of  m += v*v' with v a column vector
template<typename TOut, typename TInVec>
void
addOuterProdSym(
  prlite::RowMatrix<TOut>& m,
  const TInVec& v){
  auto vPtr = v.begin();
  auto mrPtr = m.rowBegin();
  while ( vPtr != v.end() ) {
    auto wPtr = v.begin();
    auto v2End = vPtr+1;
    auto mPtr = (*mrPtr).begin();
    while (wPtr != v2End) {
      *mPtr++ += *vPtr * *wPtr++;
    } // while
    vPtr++;
    mrPtr++;
  } // while
} // addOuterProd

// lower triangle of m = v*w' with v a column vector and w = v*scalar
template<typename TOut, typename TInVec>
prlite::RowMatrix<TOut>
outerProdSym(
  const TInVec& v,
  const TInVec& w) {

  prlite::RowMatrix<TOut> m(v.size(),w.size());

  auto vPtr = v.begin();
  auto mrPtr = m.rowBegin();
  unsigned szv = 0;
  while ( vPtr != v.end() ) {
    unsigned szw = ++szv;
    auto wPtr = w.begin();
    auto mPtr = (*mrPtr).begin();
    while (szw--) {
      *mPtr++ = *vPtr * *wPtr++;
    } // while
    vPtr++;
    mrPtr++;
  } // while

  return m;
} // outerProdSym

// lower triangle of m = v*w' with v a column vector, w = v*scalar, m
// symmetrical
template<typename TOut, typename TInVec>
void
addOuterProdSym(
  prlite::RowMatrix<TOut>& m,
  const TInVec& v,
  const TInVec& w){

  auto vPtr = v.begin();
  auto mrPtr = m.rowBegin();
  unsigned szv = 0;
  while ( vPtr != v.end() ) {
    unsigned szw = ++szv;
    auto wPtr = w.begin();
    auto mPtr = (*mrPtr).begin();
    while (szw--) {
      *mPtr++ += *vPtr * *wPtr++;
    } // while
    vPtr++;
    mrPtr++;
  } // while

} // addOuterProdSym


// inplace adds m += v*w'
template<typename TOut, typename TInVec>
void
addOuterProd(
  prlite::RowMatrix<TOut>& m,
  const TInVec& v,
  const TInVec& w){

  //m += v*w.transpose();
  auto vPtr = v.begin();
  auto mrPtr = m.rowBegin();
  while ( vPtr != v.end() ) {
    auto wPtr = w.begin();
    auto mPtr = (*mrPtr).begin();
    while ( wPtr != w.end() ) {
      *mPtr++ += *vPtr * *wPtr++;
    } // while
    vPtr++;
    mrPtr++;
  } // while

} // addOuterProd

template<typename T>
prlite::ColVector<T> packLowerDiagsToVec(const prlite::RowMatrix<T>& mat) {

  unsigned dim = mat.rows();
  prlite::ColVector<T> vec(dim*(dim+1)/2);
  unsigned n = 0;
  for (unsigned r = 0; r < dim; r++) {
    for (unsigned d = r; d < dim; d++) {
      //cout << d << "," << d-r << " ";
      vec[n++] = mat(d,d-r);
    } // for d
    //cout << endl;
  } // for r

  return vec;
} // packLowerDiagsToVec

template<typename T>
prlite::RowMatrix<T> packVecToLowerDiags(const prlite::ColVector<T>& vec, unsigned dim){
  prlite::RowMatrix<T> mat(dim,dim);
  unsigned n = 0;
  for (unsigned r = 0; r < dim; r++) {
    for (unsigned d = r; d < dim; d++) {
      //cout << d << "," << d-r << " ";
      mat(d,d-r) = vec[n++];
    } // for d
    //cout << endl;
  } // for r

  return mat;
} // packVecToLowerDiags

//This function returns a subset of a Matrix consisting of the rows and vectors specified
template<typename T>
prlite::RowMatrix<T> subsetMat(prlite::RowMatrix<T> mat, std::vector<int> rows, std::vector<int> cols){
  prlite::RowMatrix<T> submat(rows.size(),cols.size());

  //Just loop through the matrix and keep what you need
  for (unsigned i = 0; i < rows.size(); i++){
    for (unsigned j = 0; j < cols.size(); j++){
      submat(i,j) = mat(rows[i], cols[j]);
    }
  }
  return submat;
}

//This function returns a subset of a RowVector
template<typename T>
prlite::ColVector<T> subsetVec(prlite::ColVector<T> vec, std::vector<int> rows){
  prlite::ColVector<T> subvec(rows.size());

  //Just loop through the vector and keep what you need
  for (unsigned i = 0; i < rows.size(); i++){
      subvec[i] = vec[rows[i]];
  }
  return subvec;
}
