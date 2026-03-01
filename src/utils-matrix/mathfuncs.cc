/*******************************************************************************

          AUTHORS: DJP Badenhorst, JA du Preez
          COPYRIGHT: University of Stellenbosch, all rights reserved.
          DATE:    Sept 2014
          PURPOSE: Various mathematical and probability-related functions

*******************************************************************************/

#include "mathfuncs.hpp"

// digamma(1.0) = 0.5772156649015328606065120900824024310421;
// The digamma function is called in the calculation of the distance between two Dirichlet potentials
double digamma(double input){
  // Ensure that the argument given to the digamma function is positive

  if (input > 0.0) {

    double output=0;
    // If the given argument is below 6, use recurrent formula to express as a function of digamma(x) where x>6
    double x=input;
    if(x<=1){
      output=output-1/x-1/(x+1)-1/(x+2)-1/(x+3)-1/(x+4)-1/(x+5);
      x=x+6;
    }
    else if(x<=2){
      output=output-1/x-1/(x+1)-1/(x+2)-1/(x+3)-1/(x+4);
      x=x+5;
    }
    else if(x<=3){
      output=output-1/x-1/(x+1)-1/(x+2)-1/(x+3);
      x=x+4;
    }
    else if(x<4){
      output=output-1/x-1/(x+1)-1/(x+2);
      x=x+3;
    }
    else if(x<5){
      output=output-1/x-1/(x+1);
      x=x+2;
    }
    else if(x<6){
      output=output-1/x;
      x=x+1;
    }
    // Once we have ensured that the x is greater than 6, use below equation to approximate digamma(x)
    output = output+log(x)-1/(2*x)-1/(12*pow(x,2))+1/(120*pow(x,4))-1/(252*pow(x,6))+1/(240*pow(x,8))-5/(660*pow(x,10))+691/(32760*pow(x,12))-1/(12*pow(x,14));
    return(output);

  } // if
  else if (input < 0.0) {
    PRLITE_ASSERT(false,"Digamma function can only be applied to positive numbers:" << input);
  } // else if

  // according to WolframAlpha actually complex infinity. But looking
  // at eps->0+ values I think this will work.
  return -std::numeric_limits<double>::infinity();

  } // digamma

double logSumExp(double l1, double l2) {
  if (l1 > l2) {
    return l1 + log( 1.0+exp(l2-l1) );
  } // if
  else {
    return l2 + log( 1.0+exp(l1-l2) );
  } // else
} // logSumExp
