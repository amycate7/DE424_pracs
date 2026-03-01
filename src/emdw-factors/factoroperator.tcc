/*
 * Author     : JA du Preez
 * Created on : /02/2013
 * Copyright  : University of Stellenbosch, all rights retained
 */

//------------------Family 0

//TO BE DEPRECATED SOMETIME
template<typename FactorType>
emdw::RVVals Operator1<FactorType>::sample(const FactorType*) {
  std::cout << "sample(FactorType*) not supported\n";
  throw "sample(FactorType*) not supported";
} // sample

//------------------Family 1

//TO BE DEPRECATED SOMETIME
template<typename FactorType>
Factor* Operator1<FactorType>::process(const FactorType*){
  std::cout << "process(FactorType*) not supported\n";
  throw "process(FactorType*) not supported";
} // process

//TO BE DEPRECATED SOMETIME
template<typename FactorType>
void Operator1<FactorType>::inplaceProcess(FactorType*){
  std::cout << "inplaceProcess(FactorType*) not supported\n";
  throw "inplaceProcess(FactorType*) not supported";
} // inplaceProcess

//------------------Family 2

//TO BE DEPRECATED SOMETIME
template<typename FactorType>
Factor* Operator1<FactorType>::process(const FactorType*,
                                       const Factor*){
  std::cout << "process(FactorType*, Factor*) not supported\n";
  throw "process(FactorType*, Factor*) not supported";
} // process

//TO BE DEPRECATED SOMETIME
template<typename FactorType>
void Operator1<FactorType>::inplaceProcess(FactorType*,
                                           const Factor*){
  std::cout << "inplaceProcess(FactorType*, Factor*) not supported\n";
  throw "inplaceProcess(FactorType*, Factor*) not supported";
} // inplaceProcess

//------------------Family 3

//TO BE DEPRECATED SOMETIME
template<typename FactorType>
Factor* Operator1<FactorType>::process(const FactorType*,
                                       const emdw::RVIds&,
                                       bool,
                                       const Factor*){
  std::cout << "marginalize(FactorType*, const emdw::RVIds&, bool, const Factor*) not supported\n";
  throw "marginalize(FactorType*, const emdw::RVIds&, bool, const Factor*) not supported";
} // process

//------------------Family 4

//TO BE DEPRECATED SOMETIME
template<typename FactorType>
Factor* Operator1<FactorType>::process(const FactorType* lhsPtr,
                                       const emdw::RVIds& variables,
                                       const emdw::RVVals& assignedVals,
                                       bool presorted) {
  std::cout << "process(FactorType*, const emdw::RVIds&, const emdw::RVVals&, bool) not supported\n";
  throw "process(FactorType*, const emdw::RVIds&, const emdw::RVVals&, bool) not supported";
} // process

//------------------Family 5

//TO BE DEPRECATED SOMETIME
template<typename FactorType>
double Operator1<FactorType>::inplaceProcess(FactorType* lhsPtr,
                                             const Factor* rhsPtr,
                                             double df) {
  std::cout << "inplaceProcess(FactorType* lhsPtr, const Factor* rhsPtr, double df) not supported\n";
  throw "inplaceProcess(FactorType* lhsPtr, const Factor* rhsPtr, double df) not supported";
} // inplaceProcess


// AND THIS IS THE NEW SET

//------------------Family 0

template<typename FactorType>
emdw::RVVals SampleOperator<FactorType>::sample(const FactorType*) {
  std::cout << "sample(FactorType*) not supported\n";
  throw "sample(FactorType*) not supported";
} // sample

//------------------Family 1

template<typename FactorType>
Factor* NormalizeOperator<FactorType>::process(const FactorType*){
  std::cout << "process(FactorType*) not supported\n";
  throw "process(FactorType*) not supported";
} // process

template<typename FactorType>
void InplaceNormalizeOperator<FactorType>::inplaceProcess(FactorType*){
  std::cout << "inplaceProcess(FactorType*) not supported\n";
  throw "inplaceProcess(FactorType*) not supported";
} // inplaceProcess

//------------------Family 2

template<typename FactorType>
Factor* AbsorbOperator<FactorType>::process(const FactorType*,
                                       const Factor*){
  std::cout << "process(FactorType*, Factor*) not supported\n";
  throw "process(FactorType*, Factor*) not supported";
} // process

template<typename FactorType>
Factor* CancelOperator<FactorType>::process(const FactorType*,
                                       const Factor*){
  std::cout << "process(FactorType*, Factor*) not supported\n";
  throw "process(FactorType*, Factor*) not supported";
} // process

template<typename FactorType>
void InplaceAbsorbOperator<FactorType>::inplaceProcess(FactorType*,
                                           const Factor*){
  std::cout << "inplaceProcess(FactorType*, Factor*) not supported\n";
  throw "inplaceProcess(FactorType*, Factor*) not supported";
} // inplaceProcess

template<typename FactorType>
void InplaceCancelOperator<FactorType>::inplaceProcess(FactorType*,
                                           const Factor*){
  std::cout << "inplaceProcess(FactorType*, Factor*) not supported\n";
  throw "inplaceProcess(FactorType*, Factor*) not supported";
} // inplaceProcess

//------------------Family 3

template<typename FactorType>
Factor* MarginalizeOperator<FactorType>::process(const FactorType*,
                                                 const emdw::RVIds&,
                                                 bool,
                                                 const Factor*){
  std::cout << "marginalize(FactorType*, const emdw::RVIds&, bool, const Factor*) not supported\n";
  throw "marginalize(FactorType*, const emdw::RVIds&, bool, const Factor*) not supported";
} // process

//------------------Family 4

template<typename FactorType>
Factor* ObserveAndReduceOperator<FactorType>::process(const FactorType* lhsPtr,
                                       const emdw::RVIds& variables,
                                       const emdw::RVVals& assignedVals,
                                       bool presorted) {
  std::cout << "process(FactorType*, const emdw::RVIds&, const emdw::RVVals&, bool) not supported\n";
  throw "process(FactorType*, const emdw::RVIds&, const emdw::RVVals&, bool) not supported";
} // process

//------------------Family 5

template<typename FactorType>
double InplaceDampenOperator<FactorType>::inplaceProcess(FactorType* lhsPtr,
                                             const Factor* rhsPtr,
                                             double df) {
  std::cout << "inplaceProcess(FactorType* lhsPtr, const Factor* rhsPtr, double df) not supported\n";
  throw "inplaceProcess(FactorType* lhsPtr, const Factor* rhsPtr, double df) not supported";
} // inplaceProcess
