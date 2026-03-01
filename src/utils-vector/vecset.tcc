/*******************************************************************************

          AUTHOR:  JA du Preez (Copyright University of Stellenbosch, all rights reserved.)
          DATE:    March 2013
          PURPOSE: std::vector union and intersection routines

*******************************************************************************/

// light weight and fast
template<typename T>
bool contains(const std::vector<T>& vecL,
              const std::vector<T>& vecR) {

  auto iL = vecL.begin();
  auto iR = vecR.begin();
  while ( iL != vecL.end() && iR != vecR.end() ){

    if (*iL > *iR) { // only in vecR
      return false;
    } // else if

    else if (*iL < *iR) {      // only in vecL
      iL++;
    } // if

    else {                            //common to both vecs
      iL++; iR++;
    } // else

  } // while

  if ( iR != vecR.end() ){         // some extra spec's in vecR
    return false;
  } // while

  return true;
} // contains


template<typename T, typename TI>
std::vector<T> cupcap(const std::vector<T>& vecL,
                      const std::vector<T>& vecR,
                      std::vector<TI>& l2i,
                      std::vector<TI>& lni,
                      std::vector<TI>& r2i,
                      std::vector<TI>& rni,
                      std::vector<TI>& i2l,
                      std::vector<TI>& ni2l,
                      std::vector<TI>& ni2r) {


  size_t szL = vecL.size(); size_t szR = vecR.size();
  size_t maxI = std::max(szL,szR);
  std::vector<T> vecU; vecU.reserve(szL+szR);
  size_t iL = 0; size_t iR = 0;

  l2i.clear(); l2i.reserve(maxI);
  lni.clear(); lni.reserve(vecL.size());
  r2i.clear(); r2i.reserve(maxI);
  rni.clear(); rni.reserve(vecR.size());

  i2l.clear(); i2l.reserve(maxI);
  ni2l.clear(); ni2l.reserve(vecL.size());
  ni2r.clear(); ni2r.reserve(vecR.size());

  while (iL < szL && iR < szR){

    if (vecL[iL] < vecR[iR]) { // unique to vecL
      lni.push_back(iL);
      ni2l.push_back(vecU.size());
      vecU.push_back(vecL[iL++]);
    } // if

    else if (vecL[iL] > vecR[iR]) { // unique to vecR
      rni.push_back(iR);
      ni2r.push_back(vecU.size());
      vecU.push_back(vecR[iR++]);
    } // else if

    else { //common to both
      l2i.push_back(iL);
      i2l.push_back(vecU.size());
      vecU.push_back(vecL[iL++]);
      r2i.push_back(iR++);
    } // else

  } // while

  // fill in the surplus
  while (iL < szL){
    lni.push_back(iL);
    ni2l.push_back(vecU.size());
    vecU.push_back(vecL[iL++]);
  } // while

  while (iR < szR){
    rni.push_back(iR);
    ni2r.push_back(vecU.size());
    vecU.push_back(vecR[iR++]);
  } // while

  return vecU;
} // cupcap

template<typename T, typename TI>
std::vector<T> sortedUnion(const std::vector<T>& vecL,
                           const std::vector<T>& vecR,
                           std::vector<TI>& u2l,
                           std::vector<TI>& u2r) {
  size_t szL = vecL.size(); size_t szR = vecR.size();
  size_t iL = 0; size_t iR = 0;
  std::vector<T> vecU; vecU.reserve(szL+szR);
  u2l.clear(); u2l.reserve( vecL.size() );
  u2r.clear(); u2r.reserve( vecR.size() );

  while (iL < szL && iR < szR){

    if (vecL[iL] < vecR[iR]) { // unique to vecL
      u2l.push_back( vecU.size() );
      vecU.push_back(vecL[iL++]);
    } // if

    else if (vecL[iL] > vecR[iR]) { // unique to vecR
      u2r.push_back( vecU.size() );
      vecU.push_back(vecR[iR++]);
    } // else if

    else { //common to both
      u2l.push_back( vecU.size() );
      u2r.push_back( vecU.size() );
      vecU.push_back(vecL[iL++]);
      iR++;
    } // else

  } // while

  // fill in the surplus
  while (iL < szL){
    u2l.push_back( vecU.size() );
    vecU.push_back(vecL[iL++]);
  } // while

  while (iR < szR){
    u2r.push_back( vecU.size() );
    vecU.push_back(vecR[iR++]);
  } // while

  return vecU;
} // sortedUnion

template<typename T, typename TI>
std::vector<T> sortedUnion(const std::vector<T>& vecL,
                           const std::vector<T>& vecR,
                           std::vector<TI>& l2u,
                           std::vector<TI>& u2l,
                           std::vector<TI>& r2u,
                           std::vector<TI>& u2r) {
  size_t szL = vecL.size(); size_t szR = vecR.size();
  size_t iL = 0; size_t iR = 0;
  std::vector<T> vecU; vecU.reserve(szL+szR);
  l2u.clear(); l2u.reserve( vecL.size() );
  u2l.clear(); u2l.reserve( vecL.size() );
  r2u.clear(); r2u.reserve( vecR.size() );
  u2r.clear(); u2r.reserve( vecR.size() );

  while (iL < szL && iR < szR){

    if (vecL[iL] < vecR[iR]) { // unique to vecL
      l2u.push_back(iL);
      u2l.push_back( vecU.size() );
      vecU.push_back(vecL[iL++]);
    } // if

    else if (vecL[iL] > vecR[iR]) { // unique to vecR
      r2u.push_back(iR);
      u2r.push_back( vecU.size() );
      vecU.push_back(vecR[iR++]);
    } // else if

    else { //common to both
      l2u.push_back(iL);
      u2l.push_back( vecU.size() );
      //r2u.push_back(iR);
      //u2r.push_back( vecU.size() );
      vecU.push_back(vecL[iL++]);
      iR++;
    } // else

  } // while

  // fill in the surplus
  while (iL < szL){
    l2u.push_back(iL);
    u2l.push_back( vecU.size() );
    vecU.push_back(vecL[iL++]);
  } // while

  while (iR < szR){
    r2u.push_back(iR);
    u2r.push_back( vecU.size() );
    vecU.push_back(vecR[iR++]);
  } // while

  return vecU;
} // sortedUnion

// also check out version without unique/difference indices (below)
template<typename T, typename TI>
std::vector<T> sortedIntersection(const std::vector<T>& vecL,
                                  const std::vector<T>& vecR,
                                  std::vector<TI>& l2i,
                                  std::vector<TI>& r2i,
                                  std::vector<TI>& lni,
                                  std::vector<TI>& rni) {
  size_t szL = vecL.size(); size_t szR = vecR.size();
  size_t iL = 0; size_t iR = 0;
  std::vector<T> vecI; vecI.reserve( std::max(szL,szR) );

  l2i.clear(); l2i.reserve(vecI.size());
  r2i.clear(); r2i.reserve(vecI.size());
  lni.clear(); lni.reserve(szL);
  rni.clear(); rni.reserve(szR);

  while (iL < szL && iR < szR){

    if (vecL[iL] < vecR[iR]) { // unique to vecL
      lni.push_back(iL);
      iL++;
    } // if

    else if (vecL[iL] > vecR[iR]) { // unique to vecR
      rni.push_back(iR);
      iR++;
    } // if

    else {
      vecI.push_back(vecL[iL]); // common to both
      l2i.push_back(iL++);
      r2i.push_back(iR++);
    } // else

  } // while

  // fill in the surplus
  while (iL < szL){
    lni.push_back(iL++);
  } // while

  while (iR < szR){
    rni.push_back(iR++);
  } // while

  return vecI;
} // sortedIntersection

// also check out version with unique/difference indices (above)
template<typename T, typename TI>
std::vector<T> sortedIntersection(const std::vector<T>& vecL,
                                  const std::vector<T>& vecR,
                                  std::vector<TI>& l2i,
                                  std::vector<TI>& r2i) {
  size_t szL = vecL.size(); size_t szR = vecR.size();
  size_t iL = 0; size_t iR = 0;
  std::vector<T> vecI; vecI.reserve( std::max(szL,szR) );

  l2i.clear(); l2i.reserve(vecI.size());
  r2i.clear(); r2i.reserve(vecI.size());

  while (iL < szL && iR < szR){

    if (vecL[iL] < vecR[iR]) { // unique to vecL
      iL++;
    } // if

    else if (vecL[iL] > vecR[iR]) { // unique to vecR
      iR++;
    } // if

    else {
      vecI.push_back(vecL[iL]); // common to both
      l2i.push_back(iL++);
      r2i.push_back(iR++);
    } // else

  } // while

  return vecI;
} // sortedIntersection



template<typename T, typename TI>
void sortedSetInsert(std::vector<T> &vecL,
                     const TI &elem){
  auto i = std::lower_bound(vecL.begin(), vecL.end(), elem);
  if (i == vecL.end() || elem < *i){
    vecL.insert(i, elem);
  }
}

template <class T, class... Ts>
std::vector<T> make_vector(T&& first, Ts&&... args) {
  using first_type = std::decay_t<T>;
  return std::vector<first_type>{
      std::forward<T>(first),
      std::forward<Ts>(args)...
  };
}
