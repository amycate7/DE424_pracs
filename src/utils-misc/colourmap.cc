/*******************************************************************************

          AUTHOR:  JA du Preez
          COPYRIGHT: University of Stellenbosch, all rights reserved.
          DATE:    March 2013
          PURPOSE: Various image processing routines

*******************************************************************************/

// must be first for opencv REAL
#include "colourmap.hpp"

#include "prlite_vector.hpp" // FIXME must be before prlite_stlvecs.hpp
#include "prlite_matrix.hpp"

using namespace std;

/*******************************************************************************************/
/*               Colours in RGB format                                                     */
/*******************************************************************************************/

const std::vector<float> BLACK{0,0,0};
const std::vector<float> RED{1,0,0};
const std::vector<float> GREEN{0,1,0};
const std::vector<float> BLUE{0,0,1};
const std::vector<float> CYAN{0,1,1};
const std::vector<float> MAGENTA{1,0,1};
const std::vector<float> YELLOW{1,1,0};
const std::vector<float> WHITE{1,1,1};

/*******************************************************************************************/
/*               Various mappings of range 0:1 numbers to RGB colours                      */
/*******************************************************************************************/

void
colourRange(
  const vector< vector<float> >& colSeq,
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  unsigned nSections = colSeq.size()-1;

  normVal = pow(normVal, warp);
  unsigned sectionNr = unsigned(normVal*nSections);
  float lower = float(sectionNr)/nSections;
  float upper = float(sectionNr+1)/nSections;

  if (sectionNr == nSections) {
    R = colSeq[nSections][0];
    G = colSeq[nSections][1];
    B = colSeq[nSections][2];
  } // else
  else {
    R = (upper-normVal)/(upper-lower)*colSeq[sectionNr][0] + (normVal-lower)/(upper-lower)*colSeq[sectionNr+1][0];
    G = (upper-normVal)/(upper-lower)*colSeq[sectionNr][1] + (normVal-lower)/(upper-lower)*colSeq[sectionNr+1][1];
    B = (upper-normVal)/(upper-lower)*colSeq[sectionNr][2] + (normVal-lower)/(upper-lower)*colSeq[sectionNr+1][2];
  } // else

} // colourRange

void
colourRangeWK(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{WHITE,BLACK};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeWB

void
colourRangeKW(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{BLACK,WHITE};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeWB

void
colourRangeBCWYR(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{BLUE,CYAN,WHITE,YELLOW,RED};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeBCWYR

void
colourRangeRWB(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{RED,WHITE,BLUE};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeCBM

void
colourRangeBWR(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{BLUE,WHITE,RED};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeCBM

void
colourRangeCBM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{CYAN,BLUE,MAGENTA};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeCBM

void
colourRangeCM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{CYAN,MAGENTA};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeCM

void
colourRangeGCBM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{GREEN,CYAN,BLUE,MAGENTA};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeGCBM

void
colourRangeRYWCB(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{RED,YELLOW,WHITE,CYAN,BLUE};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeRYWCB

void
colourRangeWCB(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{RED,YELLOW,WHITE,CYAN,BLUE};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeWCB

void
colourRangeWCBM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{WHITE,CYAN,BLUE,MAGENTA};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeWCBM

void
colourRangeMBCW(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{MAGENTA,BLUE,CYAN,WHITE};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeMBCW

void
colourRangeWGCBM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{WHITE,GREEN,CYAN,BLUE,MAGENTA};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeWGCBM

void
colourRangeWYGCBM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{WHITE,YELLOW,GREEN,CYAN,BLUE,MAGENTA};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeWYGCBM

void
colourRangeYGCBM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp) {

  vector< vector<float> > seq{YELLOW,GREEN,CYAN,BLUE,MAGENTA};
  colourRange(seq, normVal, R, G, B, warp);
} // colourRangeYGCBM


/*******************************************************************************************/
/*                    Creates a uniformly coloured BGR image.                              */
/*******************************************************************************************/

std::vector< prlite::RowMatrix<float> >
filledImageBGR(
  unsigned nRows,
  unsigned nCols,
  std::vector<double> fillBGR) {

  std::vector< prlite::RowMatrix<float> > bgr(3);
  bgr[0].resize(nRows,nCols);
  bgr[1].resize(nRows,nCols);
  bgr[2].resize(nRows,nCols);
  unsigned nChans = bgr.size();

  for (unsigned r = 0; r < nRows; r++) {
    for (unsigned c = 0; c < nCols; c++) {
      for (unsigned k = 0; k < nChans; k++) {
        bgr[k](r,c) = fillBGR[k];
      } // for k
    } // for c
  } // for r

  return bgr;
} // filledImageBGR
