/*
 * Author     :  (DSP Group, E&E Eng, US)
 * Created on :
 * Copyright  : University of Stellenbosch, all rights retained
 *
 * Hash entries of a vector, extended for arbitrary byte size from
 * https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
 */

#ifndef HASHFUNCTIONS_HPP
#define HASHFUNCTIONS_HPP

#include <vector>
#include <cstdint>
#include <cstring>

/**
 * @brief Calculate a 64bit hash for a given vector of any type
 * that can forcefully be byte-accessible.
 *
 * @param vec The vector from which the hash is calculated.
 *
 * @param seed A seed value to start the hash from (useful for
 * combining a previous hash).
 *
 * @return the 64bit hash as uint64_t (usually same as ulong)
 */
template<typename T>
uint64_t vectorToHash(std::vector<T> const& vec,
                      uint64_t seed = 0){
  auto tSize = sizeof(T);

  uint64_t hash = seed; //64bits = 8bytes
  uint64_t val = 0    ; //64bits = 8bytes
  char *valAddr  = (char*)(&val); //access to raw bytes

  //Keeping the tumble code seperate from doing the actual tumbling
  auto tumble = [&hash, &val, &valAddr](char* dataAddr, size_t BSize){
    auto BLong = sizeof(uint64_t);
    for (auto i = 0*BSize; i<BSize; i+=BLong) {
      auto j = (i+BLong > BSize) ? (BSize) : (i+BLong); //(cond)?(true->last byte):(false->i+8)

      val*=0;
      memcpy(valAddr, dataAddr+i, j-i);
      hash ^= val + 0b1000101000110000110110100000101101000101010000101010111110011100 //random
              + (hash << 50) + (hash << 10) + (hash >> 3); //adapted to uint64_t (not crypto secure)
    }
  };

  //Start tumbing with the size of the vec, then with the actual data
  auto vec_seed = (uint64_t)(vec.size());
  tumble((char*)(&vec_seed), sizeof(vec_seed));

  for(auto& i : vec) {
    tumble((char*)(&i), tSize);
  }

  return hash;
}

#endif //HASHFUNCTIONS_HPP
