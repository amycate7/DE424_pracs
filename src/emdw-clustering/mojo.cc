#include "mojo.hpp"

#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <functional>
#include <random>

using namespace std;
using namespace emdw;

/******************************************************************************/

double jaccardIndex(size_t sz1, size_t sz2, size_t iSz) {
  return double(iSz)/(sz1+sz2-iSz);
} // jaccardIndex

double jaccardDistance(size_t sz1, size_t sz2, size_t iSz) {
  return -log2( jaccardIndex(sz1,sz2,iSz) );
} // jaccardDistance

/******************************************************************************/

/******************************************************************************/

void Mojo::setMaxItr() {

  maxItr = toC.begin();
  if ( maxItr == toC.end() ){
    //cout << __FILE__ << __LINE__ << " empty list" << endl;
    return;
  } // if

  auto mapItr = next(maxItr);
  while ( mapItr != toC.end() ) {
    if ( mapItr->second > maxItr->second ) {
      maxItr = mapItr;
    } // if
    mapItr++;
  } // while

} // setMaxItr

Mojo::Mojo(size_t theFromC): fromC(theFromC), toC(), maxItr( toC.end() ) {
  //cout << __FILE__ << __LINE__ << " Mojo default constructor\n";
} // def ctor

Mojo::Mojo(const Mojo& gr)
    : fromC(gr.fromC), toC(gr.toC), maxItr() {
  //cout << __FILE__ << __LINE__ << " Mojo copy constructor\n";
  setMaxItr();
} // cpy ctor

//TODO: Clang-Tidy: Move constructor initializes class member by calling a copy constructor
//TODO: Clang-Tidy: Move constructors should be marked noexcept
Mojo::Mojo(Mojo&& gr)
    : fromC(gr.fromC), toC(gr.toC), maxItr() {
  //cout << __FILE__ << __LINE__ << " Mojo move constructor\n";
  setMaxItr();
} // mv ctor

Mojo& Mojo::operator=(const Mojo& d){
  //cout << __FILE__ << __LINE__ << " Mojo copy assignment\n";
  fromC = d.fromC;
  toC = d.toC;
  setMaxItr();
  return *this;
} // operator=

Mojo& Mojo::operator=(Mojo&& d){
  //cout << __FILE__ << __LINE__ << " Mojo move assignment\n";
  fromC = d.fromC;
  toC = std::move(d.toC);
  //maxItr = d.maxItr;
  setMaxItr();
  return *this;
} // operator=

void Mojo::insert(size_t toCluster, double sc) {
  //cout << __FILE__ << __LINE__ << "\n";
  bool ok; std::map<size_t, double>::iterator insItr;
  std::tie(insItr, ok) = toC.insert(pair<size_t,double>{toCluster, sc});
  //auto [insItr, ok] = toC.insert(pair<size_t,double>{toCluster, sc});

  if (!ok) { // an existing dest
    insItr->second = sc;

    if (sc > maxItr->second) { // new top scorer
      maxItr = insItr;
    } // if

    else if (toCluster == maxItr->first) { // lost the top scorer
      maxItr = toC.begin();
      auto itr = maxItr;
      while ( ++itr != toC.end() ) {
        if (itr->second > maxItr->second) {
          maxItr = itr;
        } // if
      } // while
    } // else if

  } // if

  else { // a new dest
    if (maxItr == toC.end() or sc > maxItr->second) {maxItr = insItr;} // if
  } // else

} // insert

void Mojo::erase(size_t key) {

  auto itr = toC.find(key);
  if ( itr != toC.end() ) {toC.erase(itr);} // if

  if (itr == maxItr) {
    setMaxItr();
  } // if

} // erase

void Mojo::replace(size_t toOld, size_t toNew, double scNew) {
  //cout << __FILE__ << __LINE__ << "\n";

  auto itr = toC.find(toOld);
  if ( itr == toC.end() ) {return;} // if

  double prevBestSc = maxItr->second;
  size_t prevBestTo = maxItr->first;
  bool ok; std::map<size_t, double>::iterator insItr;
  std::tie(insItr, ok) = toC.insert( make_pair(toNew, scNew) );
  //auto [insItr, ok] = toC.insert( make_pair(toNew, scNew) );
  if (!ok) {insItr->second = scNew;} // if

  if (scNew > prevBestSc) { // we will have to replace maxItr
    maxItr = insItr;
  } // if

  else if (toOld == prevBestTo) { // we will have to search for the new  maxItr
    maxItr = toC.begin();
    auto itr = maxItr;
    while ( ++itr != toC.end() ) {
      if (itr->second > maxItr->second) {
        maxItr = itr;
      } // if
    } // while
  } // else if

  if (toOld != toNew) {  toC.erase(itr); } // if
} // replace

void Mojo::popMax() {
  //cout << __FILE__ << __LINE__ << "\n";
  toC.erase(maxItr);
  maxItr = toC.begin();
  if ( maxItr != toC.end() ) {
    auto itr = maxItr;
    while ( ++itr != toC.end() ) {
      if (itr->second > maxItr->second) {
        maxItr = itr;
      } // if
    } // while
  } // if
} // popMax


// const std::map<size_t, double>::iterator Mojo::getMax() const {
//   return maxItr;
// } // getMax

size_t Mojo::idx1() const {
  //cout << __FILE__ << __LINE__ << "\n";
  return fromC;
} // idx1

size_t Mojo::idx2() const {
  //cout << __FILE__ << __LINE__ << "\n";
  if ( maxItr == toC.end() ) {
    cout << __FILE__ << __LINE__ << ": This should not happen\n";
  } // if
  return maxItr->first;
} // idx2

double Mojo::strength() const {
  //cout << __FILE__ << __LINE__ << "\n";
  if ( maxItr == toC.end() ) {
    return 0.0;
  } // if
  return maxItr->second;
} // strength

// bool Mojo::operator<(const Mojo& cg) const {
//   return maxItr->second > cg.maxItr->second;
// } // operator<

size_t Mojo::size() const {
  //cout << __FILE__ << __LINE__ << "\n";
  return toC.size();
} // size

void Mojo::show(unsigned level) const {
  //cout << __FILE__ << __LINE__ << "\n";
  cout << fromC << "<-" << flush;

  if ( maxItr == toC.end() ) {
    cout << "empty\n";
  } // if
  else {
    cout << maxItr->first << " : " << maxItr->second
         << " (" << size() << ")" << endl;
  } // else

  if (level) {
    for (auto el : toC) {
      cout << "\t" << el.first
                << ": " << el.second << endl;
    } // for
  } // if
} // show

Mojo::ConstItr Mojo::begin() const {
  //cout << __FILE__ << __LINE__ << "\n";
  return toC.begin();
}

Mojo::ConstItr Mojo::end() const {
  //cout << __FILE__ << __LINE__ << "\n";
  return toC.end();
}

Mojo::ConstItr Mojo::find(size_t key) const {
  //cout << __FILE__ << __LINE__ << "\n";
  return toC.find(key);
}

/******************************************************************************/

/*
  // Now we are going to start clustering factors into sets that
  // will later become separate JT's. We are going to record these
  // sets in a map since the clustering process will keep on
  // reducing them, and we don't want to renumber the cluster
  // references in accel.
  map< size_t, set<size_t> > factorsInCluster;
  map< size_t, set<RVIdType> > varsInCluster;
  for (size_t m = 0; m < M; m++) {
    factorsInCluster[m].insert({m});
    const RVIds& theVars( factors[m]->getVars() );
    varsInCluster[m].insert( theVars.begin(), theVars.end() );
  } // for


//*/

std::map<size_t, Mojo>
calcGravityMojos(
  const std::vector< rcptr<Factor> >& factors,
  const vector< map<size_t, emdw::RVIds> >& theIntersections,
  const std::map<size_t, std::set<size_t> >& theReachableFrom,
  const map<size_t,double>& massOfCluster,
  bool includeSecondary) {

  auto M = unsigned(factors.size());
  std::uniform_real_distribution<> addNoise(0, 5E-4);

  // initialize the mojos
  map<size_t, Mojo> gravityMap;
  for (size_t m = 0; m < M; m++) {
    gravityMap.insert( make_pair(m, Mojo(m)) );
  } // for

  for (size_t m = 0; m < M; m++) {
    for (auto toN : theIntersections[m]) {
      auto n = toN.first;
      size_t iSz = toN.second.size();
      double dist2 = jaccardDistance(factors[m]->noOfVars(), factors[n]->noOfVars(), iSz);
      dist2 *= dist2; dist2 += addNoise( emdw::randomEngine() );
      auto itr = massOfCluster.find(m);
      gravityMap[m].insert(n, itr->second/dist2);
      itr = massOfCluster.find(n);
      gravityMap[n].insert(m, itr->second/dist2);
    } // for
  } // for

  if (includeSecondary) {
    // The secondary accelaration. NOTE: could be expensive - triple loop.
    for (auto toN : gravityMap) { // from factor m // UPDATE THIS
      size_t n = toN.first;
      auto fromItr = theReachableFrom.find(n);
      for (size_t m : fromItr->second) { // from m via n
        for (auto k : fromItr->second) { // to k

          if (k != m) {

            auto mnOvlItr =
              ( n > m ? theIntersections[m].find(n) : theIntersections[n].find(m) );
            auto nkOvlItr =
              ( k > n ?  theIntersections[n].find(k) : theIntersections[k].find(n) );
            auto& mnOvl(mnOvlItr->second);
            auto& nkOvl(nkOvlItr->second);
            vector<RVIdType> mnkOvl;
            set_intersection(mnOvl.begin(),mnOvl.end(),nkOvl.begin(),nkOvl.end(),
                             std::inserter(mnkOvl,mnkOvl.begin()));

            size_t iSz = mnkOvl.size();
            if (iSz) {

              // now add the secondary accel
              double mnkDist2 =
                jaccardDistance(factors[m]->noOfVars(), factors[n]->noOfVars(), iSz) +
                jaccardDistance(factors[n]->noOfVars(), factors[k]->noOfVars(), iSz);
              mnkDist2 *= mnkDist2;

              Mojo& mGrav(gravityMap[n]);
              double oldGrav = mGrav.find(m)->second;
              auto itr = massOfCluster.find(k);
              mGrav.insert(m, oldGrav + itr->second/mnkDist2);

              // cout << m << ","
              //      << n << ","
              //      << k << ": "
              //      << mnkOvl << " :: "
              //      << mnkDist << endl;

            } // if


          } // if

        } // for k
      } // for n
    } // for m
  } // if includeSecondary

  return gravityMap;
} // calcGravityMojos

/******************************************************************************/

void clusterGravity(const std::map<size_t, Mojo>& mojoMap,
                    std::map<size_t, std::set<size_t> >& theReachableFrom,
                    map<size_t,double>& massOfCluster,
                    map< size_t, set<size_t> >& factorsInCluster,
                    map< size_t, set<RVIdType> >& varsInCluster,
                    size_t maxRVs) {

  MojoQueue mojoQ;
  for (auto toN : mojoMap) { // from cluster m
    mojoQ.insert(std::move(toN.second)); //TODO: Clang-Tidy: Passing result of std::move() as a const reference argument; no move will actually happen
  } // for fromM

  std::uniform_real_distribution<> addNoise(0, 5E-4);

  /*
   * In the following we use the analogy of gravitation pull to
   * cluster factors together.
   *
   * IMPORTANT STRUCTURES:
   *
   *  map< size_t, set<size_t> > factorsInCluster[m] records the set
   *  of all factors that are allocated to cluster m. Initially each
   *  factorsInCluster[m] consists of only the single factor m.
   *
   * map<size_t,double> massOfCluster[m] records the "mass"
   * (informedness) of cluster m. For single factor clusters this
   * measures how far it deviates from a uniform
   * distribution. Calculating this for multi-factor clusters will
   * be expensive, therefore we approximate it as the weighted sum
   * of the individual masses.
   *
   * map< size_t, set<RVIdType> > varsInCluster[m] records the set
   * of all RVs in cluster m (i.e. the union of the RVs in each of
   * the factors in cluster m).
   *
   * map<size_t, set<size_t> > theReachableFrom[r] mirrors mojoQ,
   * but from the other side. It lists clusters that cluster r can
   * merge with.
   *
   * MojoQueue mojoQ serves as a priority queue from where we
   * can efficiently pick the Mojo set with the highest
   * dest<-srce gravity pull, or alternatively, any other set of
   * Gravities towards that particular dest. The Mojo class
   * itself lists all the possible source clusters that can gravitate
   * towards a particular destination cluster.
   *
   * Initialise the above structures using earlier existing
   * info. Then clear all earlier info.
   */

  size_t cnt = 0;
  while ( mojoQ.size() ) {

    // cout << __FILE__ << __LINE__ << endl;
    // cout << "\nVars in clusters\n";
    // for (auto m : massOfCluster) {
    //   cout << right << setw(3) << m.first << ": "
    //        << left << setw(10) << m.second << " :: "
    //        << varsInCluster[m.first] << " (" << factorsInCluster[m.first].size() << ")\n";
    // } // for
    // cout << "\nNo of clusters: " << massOfCluster.size() << "\n";

    // cout << "Reachable from:\n";
    // for (auto m : theReachableFrom) {
    //   cout << m.first << ": ";
    //   for (auto n : m.second) {
    //     cout << n << " ";
    //   } // for
    //   cout << endl;
    // } // for

    // // print gravity queue
    // cout << __FILE__ << __LINE__ << "\nQ:\n";
    // mojoQ.show(0);

    //if (cnt) exit(-1);

    // 1. Remove the top-scoring Mojo destCluster from
    // mojoQ. If there are none to pick, stop.
    //cout << __FILE__ << __LINE__ << endl;

    Mojo maxGravToDest = mojoQ.popMax();
    if ( !maxGravToDest.size() ) {
      cout << "Not linked\n";
      continue;
    } // if not linked at all
    auto destCluster = maxGravToDest.idx1();
    auto srceCluster = maxGravToDest.idx2();
    //cout << "Considering to merge " << destCluster << "<-" << srceCluster << endl;

    // 2. Check that merging the best srceCluster->destCluster
    // results in a valid RV union size. If not valid, remove this
    // srceCluster->destCluster combination from this Mojo set. If
    // there are remaining gravities in this set, re-insert it into
    // the mojoQ. Restart from step 1.
    //cout << __FILE__ << __LINE__ << endl;

    set<RVIdType> rvUnion;
    set_union(varsInCluster[destCluster].begin(),varsInCluster[destCluster].end(),
              varsInCluster[srceCluster].begin(),varsInCluster[srceCluster].end(),
              std::inserter(rvUnion, rvUnion.begin()));
    auto uSz = rvUnion.size();

    if (uSz <= maxRVs) { // we now merge cluster srceCluster into cluster destCluster
      //cout << "\nMerge " << srceCluster << "->" << destCluster << endl;

      // 3. Add mass of srceCluster to that of destCluster. Remove
      // mass of srceCluster. Merge factors in srceCluster into
      // those of destCluster. Remove factors in srceCluster. Merge
      // vars in srceCluster] into those of destCluster. Remove vars
      // in srceCluster.
      //cout << __FILE__ << __LINE__ << endl;

      massOfCluster[destCluster] = 1.0*massOfCluster[destCluster] + 1.0*massOfCluster[srceCluster];
      massOfCluster.erase(srceCluster);
      factorsInCluster[destCluster].insert( factorsInCluster[srceCluster].begin(), factorsInCluster[srceCluster].end() );
      factorsInCluster.erase(srceCluster);
      varsInCluster[destCluster].insert( varsInCluster[srceCluster].begin(), varsInCluster[srceCluster].end() );
      varsInCluster.erase(srceCluster);

      // 4. Remove the Mojo set toSrceCluster from the
      // mojoQ. Initialise a new toDestCluster Mojo set.
      //cout << __FILE__ << __LINE__ << endl;

      Mojo oldGravToSrce( mojoQ.popIdx1(srceCluster) );
      Mojo newGravToDest(destCluster);

      // 5. Update all the original dest<-k, except k==srce
      //cout << __FILE__ << __LINE__ << endl;

      for (auto itr = maxGravToDest.begin(); itr != maxGravToDest.end(); itr++) {
        size_t k = itr->first;
        if (k != srceCluster) {
          double jDist2 = jaccardDistance(varsInCluster[destCluster], varsInCluster[k]);
          jDist2 *= jDist2; jDist2 += addNoise( emdw::randomEngine() );
          newGravToDest.insert(k, massOfCluster[destCluster]/jDist2);
        } // if
      } // for

        // 6. Update all the original k<-dest, except k == srce
        //cout << __FILE__ << __LINE__ << endl;

      for (auto k : theReachableFrom[destCluster]) {
        if (k != srceCluster) {

          // Remove the Mojo set toK from the mojoQ. Update
          // k<-dest and re-insert into the mojoQ.
          Mojo gravToK = mojoQ.popIdx1(k);
          double jDist2 = jaccardDistance(varsInCluster[destCluster], varsInCluster[k]);
          jDist2 *= jDist2; jDist2 += addNoise( emdw::randomEngine() );
          gravToK.insert(destCluster, massOfCluster[k]/jDist2);
          mojoQ.insert(gravToK);

        } // if
      } // for k

        // 7. Redirect all srce<-k to dest<-k, except k==dest or k
        // already handled in previous step 5.
      //cout << __FILE__ << __LINE__ << endl;

      for (auto itr = oldGravToSrce.begin(); itr != oldGravToSrce.end(); itr++) {
        size_t k = itr->first;
        if (k != destCluster) {
          theReachableFrom[k].erase(srceCluster);
          theReachableFrom[k].insert(destCluster);
          if (newGravToDest.find(k) == newGravToDest.end() ) {
            double jDist2 = jaccardDistance(varsInCluster[destCluster], varsInCluster[k]);
            jDist2 *= jDist2; jDist2 += addNoise( emdw::randomEngine() );
            newGravToDest.insert(k, massOfCluster[destCluster]/jDist2);
          } // if
        } // for
      } // for

        // 8. Redirect all k<-srce to k<-dest, except k==dest or k
        // already handled in previous step 6.
        //cout << __FILE__ << __LINE__ << endl;

      for (auto k : theReachableFrom[srceCluster]) {

        if (k != destCluster) {
          theReachableFrom[destCluster].insert(k);

          // Remove the Mojo set toK from the mojoQ. Calculate
          // the gravity k<-dest and REPLACE k<-srce with
          // this. Re-insert into the mojoQ.
          Mojo gravToK = mojoQ.popIdx1(k);
          gravToK.erase(srceCluster);
          if ( gravToK.find(destCluster) == gravToK.end() ) {
            double jDist2 = jaccardDistance(varsInCluster[destCluster], varsInCluster[k]);
            jDist2 *= jDist2; jDist2 += addNoise( emdw::randomEngine() );
            gravToK.insert(destCluster, massOfCluster[k]/jDist2);
          } // if
          mojoQ.insert(gravToK);

        } // if
      } // for k

        // 9. Now we can remove srceCluster from the reachable clusters
        //cout << __FILE__ << __LINE__ << endl;

      theReachableFrom.erase(srceCluster);
      theReachableFrom[destCluster].erase(srceCluster);

      // 10. If the new dest Gravities are non-empty, put it back in the
      // mojoQ and repeat from step 1.
      //cout << __FILE__ << __LINE__ << endl;

      if ( newGravToDest.size() ) {mojoQ.insert(newGravToDest);} // if

    } // if ok to merge

    else { // could not merge this one, too many variables in destCluster
      //cout << __FILE__ << __LINE__ << endl;
      // cout << "Could not merge cluster " << srceCluster
      //      << " with cluster " << destCluster << endl;
      // cout << "Result would have had " << uSz
      //      << " RVs which exceeds the max of " << maxRVs
      //      << endl;

      maxGravToDest.popMax();
      theReachableFrom[srceCluster].erase(destCluster);
      if ( maxGravToDest.size() ) {mojoQ.insert(maxGravToDest);} // if
    } // else

    cnt++;

  } // while


  } // clusterGravity

/******************************************************************************/

void primsAlgo(std::map<size_t, Mojo>& available, std::vector< Idx2 >& mst) {

  //A single node is already a tree.
  if (available.size() < 2) {return;} // if

  // -------------------- initialization
      // cout << "available: " << endl;
      // for (auto el: available) {
      //    cout << el.first << ": "; el.second.show(1);
      // } // for

  auto maxItr = available.begin();
  if ( maxItr == available.end() ) {
    cout << __FILE__ << __LINE__ << ": Hey come on, that's not funny!\n";
  } // if
  //cout << "1Max at "; maxItr->second.show(1);
  //exit(-1);

  // find the strongest attraction
  auto itr = next(maxItr);
  while ( itr != available.end() ) {
    if ( itr->second.strength() > maxItr->second.strength() ) {
      maxItr = itr;
    } // if
    itr++;
  } //while

  //cout << "2Max at "; maxItr->second.show(1);
  size_t idx1 = maxItr->second.idx1();
  size_t idx2 = maxItr->second.idx2();


  map<size_t, Mojo> inTree;

  // cout << __FILE__ << __LINE__ << "\n";
  // cout << "Processing " << idx1 << "<-" << idx2 << endl;
  // remove idx1 and idx2 from the "key-side" in the available set,
  // and add them to the tree

  // cout << __FILE__ << __LINE__ << "\n";
  // cout << "available: " << endl;
  // for (auto el: available) {
  //   cout << el.first << ": "; el.second.show(1);
  // } // for
  itr = available.find(idx1);
  if ( itr == available.end() ){cout << "WTF\n";}
  inTree.insert( make_pair(idx1,itr->second) );
  //cout << __FILE__ << __LINE__ << "\n";
  available.erase(itr);

  // cout << __FILE__ << __LINE__ << "\n";
  // cout << "available: " << endl;
  // for (auto el: available) {
  //   cout << el.first << ": "; el.second.show(1);
  // } // for
  itr = available.find(idx2);
  if ( itr == available.end() ){cout << "WTF\n";}
  //cout << __FILE__ << __LINE__ << "\n";
  inTree.insert( make_pair(idx2,itr->second) );
  //cout << __FILE__ << __LINE__ << "\n";
  available.erase(itr);
      // cout << __FILE__ << __LINE__ << "\n";
      // cout << "available: " << endl;
      // for (auto el: available) {
      //    cout << el.first << ": "; el.second.show(1);
      // } // for

  // remove idx1 and idx2 from the "value/mojo-side" in the available set
  for (auto itr = available.begin(); itr != available.end(); itr++) {
    //cout << __FILE__ << __LINE__ << "\n";
    itr->second.erase(idx1);
    //cout << __FILE__ << __LINE__ << "\n";
    itr->second.erase(idx2);
    //cout << __FILE__ << __LINE__ << "\n";
  } // for

  //cout << __FILE__ << __LINE__ << "\n";
  // now remove idx1 and idx2 from the "value/mojo-side" in the tree
  // while also noting the key that has the maximum strength.
  maxItr = inTree.begin();
  if ( maxItr == inTree.end() ) {
    cout << __FILE__ << __LINE__ << ": Hey come on, that's not funny!\n";
  } // if
  maxItr->second.erase(idx1);
  maxItr->second.erase(idx2);
  // find the strongest attraction
  itr = next(maxItr);
  //cout << __FILE__ << __LINE__ << "\n";
  while ( itr != inTree.end() ) {
    //cout << __FILE__ << __LINE__ << "\n";
    itr->second.erase(idx1);
    itr->second.erase(idx2);

    if ( itr->second.size() ) {
      if ( itr->second.strength() > maxItr->second.strength() ) {
        maxItr = itr;
      } // if
      itr++;
    } // if

    else {
      itr = inTree.erase(itr);
    } // else

      //cout << __FILE__ << __LINE__ << "\n";
  } //while

  // first two nodes are linked
  if (idx2 > idx1) {mst.push_back( make_pair(idx1, idx2) );} // if
  else {mst.push_back( make_pair(idx2, idx1) );} // else

  // -------------------- now iterate greedily

  //cout << __FILE__ << __LINE__ << "\n";
  while ( available.size() and inTree.size() ){

    size_t idx1 = maxItr->second.idx1();
    size_t idx2 = maxItr->second.idx2();

    // remove idx2 from the "key-side" in the available set, and
    // insert it into the tree
    itr = available.find(idx2);
    inTree[idx2] = itr->second;
    available.erase(itr);

    // remove idx2 from the "value/mojo-side" in the available set
    for (auto itr = available.begin(); itr != available.end(); itr++) {
      itr->second.erase(idx2);
    } // for

    // now remove idx2 from the "value/mojo-side" in the tree while
    // also noting the key that has the maximum strength.
    maxItr = inTree.begin();
    if ( maxItr == inTree.end() ) {
      cout << __FILE__ << __LINE__ << ": Hey come on, that's not funny!\n";
    } // if
    maxItr->second.erase(idx2);
    // find the strongest attraction
    auto itr = next(maxItr);
    while ( itr != inTree.end() ) {
      itr->second.erase(idx2);

      if ( itr->second.size() ) {
        if ( itr->second.strength() > maxItr->second.strength() ) {
          maxItr = itr;
        } // if
        itr++;
      } // if

      else {
        itr = inTree.erase(itr);
      } // else

    } //while

    // two nodes are linked
    if (idx2 > idx1) {mst.push_back( make_pair(idx1, idx2) );} // if
    else {mst.push_back( make_pair(idx2, idx1) );} // else

  } // while

} // primsAlgo

/******************************************************************************/


/******************************************************************************/

/******************************************************************************/

void plotStructure(const map<size_t, Mojo>& accel) {
  string filename = "accel.dot";
  ofstream gvfile( filename.c_str() );

  // top stuff
  gvfile << "digraph cluster_accelarations {" << std::endl;
  gvfile << "\t" << "rankdir = LR;" << std::endl;
  gvfile << "\t" << "orientation = landscape;" << std::endl;
  // gvfile << "\t" << "concentrate = true;" << std::endl;
  gvfile << "\t" << "size = \"7.5,10.0\";" << std::endl;
  gvfile << std::endl;

  gvfile << "\t" << "node [shape=ellipse]" << std::endl;

  for (auto fromM : accel) {
    size_t m = fromM.second.idx1();
    size_t n = fromM.second.idx2();
    gvfile << "\t" << m << " -> " << n << std::endl;
  } // for

  // and terminate things
  gvfile << "}" << std::endl;
  gvfile.close();
} // plotStructure

void plotStructure(const vector<Idx2>& accel) {
  string filename = "mst.dot";
  ofstream gvfile( filename.c_str() );

  // top stuff
  gvfile << "graph prims {" << std::endl;
  gvfile << "\t" << "rankdir = LR;" << std::endl;
  gvfile << "\t" << "orientation = landscape;" << std::endl;
  // gvfile << "\t" << "concentrate = true;" << std::endl;
  gvfile << "\t" << "size = \"7.5,10.0\";" << std::endl;
  gvfile << std::endl;

  gvfile << "\t" << "node [shape=ellipse]" << std::endl;

  for (auto fromM : accel) {
    size_t m = fromM.first;
    size_t n = fromM.second;
    gvfile << "\t" << m << " -- " << n << std::endl;
  } // for

  // and terminate things
  gvfile << "}" << std::endl;
  gvfile.close();
} // plotStructure


std::set<int> abs(const std::set<int>& s1) {
  std::set<int> ret;
  for (auto el : s1) {
    ret.insert( std::abs(el) );
  } // for
  return ret;
} // abs
