/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    January 2015
          PURPOSE: Sorting via the indices of a container

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/


template <typename T, typename Lambda>
std::vector<size_t>
sortIndices(
  const T& v,
  Lambda&& comp) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort( idx.begin(), idx.end(),
        [&v,&comp](size_t i1, size_t i2) {return comp(v[i1],v[i2]);} );
  return idx;
} // sortIndices

template <typename T, typename Lambda>
std::vector<size_t>
sortIndices(
  const T& iBegin,
  const T& iFrom,
  const T& iPastTo,
  Lambda&& comp) {

  // initialize original index locations
  std::vector<size_t> idx(iPastTo-iFrom);
  size_t offset = iFrom-iBegin;
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i+offset;

  // sort indexes based on comparing values in v
  sort( idx.begin(), idx.end(),
        [&iBegin,&comp](size_t i1, size_t i2) {return comp( *(iBegin+i1),*(iBegin+i2) );} );
  return idx;
} // sortIndices
