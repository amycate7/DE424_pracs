/*
 * Author     : JA du Preez
 * Created on : /02/2013
 * Copyright  : University of Stellenbosch, all rights retained
 */

//========================= Operators and Conversions =========================
//=========================== Inplace Configuration ===========================
//===================== Required Virtual Member Functions =====================
//====================== Other Virtual Member Functions =======================
//======================= Non-virtual Member Functions ========================

//------------------Family 0

template<typename FactorType>
inline emdw::RVVals Factor::dynamicSample(
  FactorOperator* procPtr,
  const FactorType* self) const {
  return dynamic_cast<SampleOperator<FactorType>&>(*procPtr).sample(self);
} //dynamicSample

//------------------Family 1

template<typename FactorType>
inline Factor* Factor::dynamicApply(FactorOperator* procPtr,
                                    const FactorType* self) const {
  return dynamic_cast<Operator1<FactorType>&>(*procPtr).process(self);
} //dynamicApply

template<typename FactorType>
inline void Factor::dynamicInplaceApply(FactorOperator* procPtr,
                                        FactorType* self) {
  dynamic_cast<Operator1<FactorType>&>(*procPtr).inplaceProcess(self);
} //dynamicInplaceApply

//------------------Family 2

template<typename FactorType>
inline Factor* Factor::dynamicApply(FactorOperator* procPtr,
                                    const FactorType* self,
                                    const Factor* other) const {
  return dynamic_cast<Operator1<FactorType>&>(*procPtr).process(self, other);
} //dynamicApply

template<typename FactorType>
inline void Factor::dynamicInplaceApply(FactorOperator* procPtr,
                                        FactorType* self,
                                        const Factor* other){

  dynamic_cast<Operator1<FactorType>&>(*procPtr).inplaceProcess(self, other);
} //dynamicInplaceApply

//------------------Family 3

template<typename FactorType>
inline Factor* Factor::dynamicApply(FactorOperator* procPtr,
                                    const FactorType* self,
                                    const emdw::RVIds& variablesToKeep,
                                    bool presorted,
                                    const Factor* peekAheadPtr) const {
  return dynamic_cast<Operator1<FactorType>&>(*procPtr).process(self, variablesToKeep, presorted, peekAheadPtr);
} //dynamicApply

  //------------------Family 4

template<typename FactorType>
inline Factor* Factor::dynamicApply(FactorOperator* procPtr,
                                    const FactorType* self,
                                    const emdw::RVIds& variables,
                                    const emdw::RVVals& assignedVals,
                                    bool presorted) const {
  return dynamic_cast<Operator1<FactorType>&>(*procPtr).process(self, variables, assignedVals, presorted);
} //dynamicApply


template<typename FactorType>
inline void Factor::dynamicInplaceApply(FactorOperator* procPtr,
                                        FactorType* self,
                                        const emdw::RVIds& variables,
                                        const emdw::RVVals& assignedVals,
                                        bool presorted) {
  dynamic_cast<Operator1<FactorType>&>(*procPtr).inplaceProcess(self, variables, assignedVals, presorted);
} //dynamicInplaceApply


  //------------------Family 5

  template<typename FactorType>
  inline double Factor::dynamicInplaceApply(FactorOperator* procPtr,
                                    FactorType* self,
                                    const Factor* other,
                                    double df){
    return dynamic_cast<Operator1<FactorType>&>(*procPtr).inplaceProcess(self, other, df);
} //dynamicInplaceApply


//================================== Friends ==================================
