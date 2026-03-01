/*
 * Author     : JA du Preez (DSP Group, E&E Eng, US)
 * Created on : /02/2013
 * Copyright  : University of Stellenbosch, all rights retained
 */

// patrec headers
#include "emdw.hpp"  // emdw

// standard headers
#include <string>  // string
#include <iostream>  // cout, endl
//#include <typeinfo>

using namespace std;
using namespace emdw;

bool emdw::debug = false;

emdw::RandomEngine::RandomEngine(unsigned theSeedVal)
    : seedVal( theSeedVal ? theSeedVal : randDev() ),
      randEng(seedVal) {
  std::ofstream ofile("emdw_RandomEngine_seedVal.txt");
  ofile << seedVal << std::endl;
  ofile.close();
} // ctor

std::default_random_engine& emdw::RandomEngine::operator()() {
  return randEng;
} // operator()

void emdw::RandomEngine::setSeedVal(unsigned theSeedVal) {
  this->~RandomEngine();
  new(this)RandomEngine(theSeedVal);
} // setSeedVal

unsigned emdw::RandomEngine::getSeedVal() const {
  return seedVal;
} // getSeedVal

emdw::RandomEngine emdw::randomEngine;


std::ostream& operator<<( std::ostream& file, Idx2 pr ) {
  file << pr.first << " " << pr.second;
  return file;
}

std::istream& operator>>( std::istream& file, Idx2& pr ){
  file >> pr.first >> pr.second;
  return file;
}
