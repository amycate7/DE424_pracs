#ifndef NESTEDLOOP_HPP
#define NESTEDLOOP_HPP

/*******************************************************************************

  AUTHOR:  SF Streicher
  EMAIL:   sfstreicher@gmail.com
  DATE:  September 2016
  PURPOSE: A helper to loop over multiple ranges in a nested fashion.

  Code example:

    NestedLoop looper;
    std::vector<unsigned> maxes = {2, 3, 2, 2};
    looper.reset(maxes);
    int i = 0;

    std::cout << "First test" << std::endl;
    cout << "Premature indices: " << looper.idxes << " (do not use)" << endl;
    while(looper.next()){
      std::cout << "Indices: " << looper.idxes << ", Last nest incremented: " << looper.nestLevel << std::endl;
      if(i == 5){
        std::cout << "...Jump Second Nest (index 1)..." << std::endl;
        looper.jumpNest(1);
      }
      i++;
    }

    std::cout << "\nSecond test" << std::endl;
    maxes = {2, 2, 3, 2};
    looper.reset(maxes);
    while(looper.next()){
      std::cout << "Indices: " << looper.idxes << ", Last nest incremented: " << looper.nestLevel << std::endl;
    }

  Code output:

    First test
    Premature indices: 4  4294967295 2 1 1  (do not use)
    Indices: 4  0 0 0 0 , Last nest incremented: 0
    Indices: 4  0 0 0 1 , Last nest incremented: 3
    Indices: 4  0 0 1 0 , Last nest incremented: 2
    Indices: 4  0 0 1 1 , Last nest incremented: 3
    Indices: 4  0 1 0 0 , Last nest incremented: 1
    Indices: 4  0 1 0 1 , Last nest incremented: 3
    ...Jump Second Nest (index 1)...
    Indices: 4  0 2 0 0 , Last nest incremented: 1
    Indices: 4  0 2 0 1 , Last nest incremented: 3
    Indices: 4  0 2 1 0 , Last nest incremented: 2
    Indices: 4  0 2 1 1 , Last nest incremented: 3
    Indices: 4  1 0 0 0 , Last nest incremented: 0
    Indices: 4  1 0 0 1 , Last nest incremented: 3
    Indices: 4  1 0 1 0 , Last nest incremented: 2
    Indices: 4  1 0 1 1 , Last nest incremented: 3
    Indices: 4  1 1 0 0 , Last nest incremented: 1
    Indices: 4  1 1 0 1 , Last nest incremented: 3
    Indices: 4  1 1 1 0 , Last nest incremented: 2
    Indices: 4  1 1 1 1 , Last nest incremented: 3
    Indices: 4  1 2 0 0 , Last nest incremented: 1
    Indices: 4  1 2 0 1 , Last nest incremented: 3
    Indices: 4  1 2 1 0 , Last nest incremented: 2
    Indices: 4  1 2 1 1 , Last nest incremented: 3

    Second test
    Indices: 4  0 0 0 0 , Last nest incremented: 0
    Indices: 4  0 0 0 1 , Last nest incremented: 3
    Indices: 4  0 0 1 0 , Last nest incremented: 2
    Indices: 4  0 0 1 1 , Last nest incremented: 3
    Indices: 4  0 0 2 0 , Last nest incremented: 2
    Indices: 4  0 0 2 1 , Last nest incremented: 3
    Indices: 4  0 1 0 0 , Last nest incremented: 1
    Indices: 4  0 1 0 1 , Last nest incremented: 3
    Indices: 4  0 1 1 0 , Last nest incremented: 2
    Indices: 4  0 1 1 1 , Last nest incremented: 3
    Indices: 4  0 1 2 0 , Last nest incremented: 2
    Indices: 4  0 1 2 1 , Last nest incremented: 3
    Indices: 4  1 0 0 0 , Last nest incremented: 0
    Indices: 4  1 0 0 1 , Last nest incremented: 3
    Indices: 4  1 0 1 0 , Last nest incremented: 2
    Indices: 4  1 0 1 1 , Last nest incremented: 3
    Indices: 4  1 0 2 0 , Last nest incremented: 2
    Indices: 4  1 0 2 1 , Last nest incremented: 3
    Indices: 4  1 1 0 0 , Last nest incremented: 1
    Indices: 4  1 1 0 1 , Last nest incremented: 3
    Indices: 4  1 1 1 0 , Last nest incremented: 2
    Indices: 4  1 1 1 1 , Last nest incremented: 3
    Indices: 4  1 1 2 0 , Last nest incremented: 2
    Indices: 4  1 1 2 1 , Last nest incremented: 3


*******************************************************************************/

#include <vector>

//TODO: template <class T=int>
class NestedLoop {
  public:
    //Variables
    std::vector<uint> maxes;
    std::vector<uint> idxes;
    int N=0;
    int nestLevel=0;

    NestedLoop();
    NestedLoop(int numberOfNests, uint Max);
    NestedLoop(std::vector<uint> theMaxes);

    void reset(int numberOfNests, uint Max);
    void reset(std::vector<uint> theMaxes);

    bool next();
    void jumpNest(int theNest);

  private:
    void clear();
  //void skipMostOuter();
};

//Initialisations
NestedLoop::NestedLoop(){}

NestedLoop::NestedLoop(int numberOfNests, uint Max) {
  reset(numberOfNests, Max);
}
NestedLoop::NestedLoop(std::vector<uint> theMaxes) {
  reset(theMaxes);
}

void NestedLoop::clear(){
  maxes.clear();
  idxes.clear();
  N = 0;
  nestLevel = 0;
}

//Reset the scene
void NestedLoop::reset(int numberOfNests, uint Max){
  std::vector<uint> theMaxes(numberOfNests);
  for(int i =0; i < numberOfNests; i++) theMaxes[i] = Max;
  reset(theMaxes);
}

void NestedLoop::reset(std::vector<uint> theMaxes){
  clear();
  N = int(theMaxes.size());

  maxes.resize(N);
  for(int i = 0; i < N; i++) maxes[i] = theMaxes[i];

  idxes.resize(N);
  idxes[0] = -1;
  for(int i=1; i<N; i++) idxes[i]=theMaxes[i]-1;
}

bool NestedLoop::next(){
  idxes[N-1]+=1;

  for(int i=N-1; i>=0; i--){
    if(idxes[i]>=maxes[i]) {
      idxes[i] = 0;

      if( i > 0 ){ //actually, if i > 0 is needed
        idxes[i-1] += 1;
      }else{
        return false;
      }
    }else{
      nestLevel = i;
      break;
    }
  }
  return true;
}

void NestedLoop::jumpNest(int theNest){
  for(int i = N-1; i>theNest; i--) {
    idxes[i] = maxes[i]-1;
  }
}

#endif // NESTEDLOOP_HPP
