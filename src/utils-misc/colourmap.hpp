#ifndef COLOURMAP_HPP
#define COLOURMAP_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez
          COPYRIGHT: University of Stellenbosch, all rights reserved.
          DATE:    March 2013
          PURPOSE: Various image processing routines

*******************************************************************************/

// standard headers
#include <iostream> // istream, ostream
#include <vector>
#include <set>
#include <map>

#include "prlite_matrix.hpp"

/*******************************************************************************************/
/*               Colours in RGB format                                                     */
/*******************************************************************************************/

extern const std::vector<float> BLACK;
extern const std::vector<float> RED;
extern const std::vector<float> GREEN;
extern const std::vector<float> BLUE;
extern const std::vector<float> CYAN;
extern const std::vector<float> MAGENTA;
extern const std::vector<float> YELLOW;
extern const std::vector<float> WHITE;

/*******************************************************************************************/
/*               Various mappings of range 0:1 numbers to RGB colours                      */
/*******************************************************************************************/

/**
 * Maps the [0..1] range onto a continuously varying White to Black
 * monochrome scale. Nice for an 'inverted' monochrome represention
 * range [0:white..1:black].
 *
 * @param normVal is the input value to be mapped to a colour. Range
 *   [0..1].
 *
 * @param R,G,B returns by reference the corresponding RGB values.
 *
 * @param warp can be used to concentrate more of the colours into the
 *   lower (0.5 <= warp < 1) or higher (1 < warp <= 2.0) areas. For
 *   most purposes the linear value (warp = 1.0) should be fine.
 */
void
colourRangeWK(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp);

void
colourRangeKW(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp);

/**
 * Maps the [0..1] range onto a continuously varying White Cyan Blue
 * Magenta colour scale. Soft on the eyes and the details at both ends
 * are nicely visible. Particularly nice for representing range [0..1]
 * values such as edges etc.
 *
 * @param normVal is the input value to be mapped to a colour. Range
 *   [0..1].
 *
 * @param R,G,B returns by reference the corresponding RGB values.
 *
 * @param warp can be used to concentrate more of the colours into the
 *   lower (0.5 <= warp < 1) or higher (1 < warp <= 2.0) areas. For
 *   most purposes the linear value (warp = 1.0) should be fine.
 */
void
colourRangeWCBM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp = 1.0);

void
colourRangeMBCW(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp = 1.0);

/**
 * Maps the [0..1] range onto a continuously varying Blue Cyan White
 * Yellow Red colour scale. With Blue on one extreme, White in the
 * middle and Red on the other extreme, this is nice for showing
 * values on the [-1..1] range. However, the values must then it must
 * first be remapped to a [0..1] range via (x+1)/2.
 *
 * @param normVal is the input value to be mapped to a colour. Range
 *   [0..1].
 *
 * @param R,G,B returns by reference the corresponding RGB values.
 *
 * @param warp can be used to concentrate more of the colours into the
 *   lower (0.5 <= warp < 1) or higher (1 < warp <= 2.0) areas. For
 *   most purposes the linear value (warp = 1.0) should be fine.
 */
void
colourRangeBCWYR(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp = 1.0);

void
colourRangeRWB(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp);

void
colourRangeBWR(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp);


void
colourRangeCBM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp);

void
colourRangeCM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp);

void
colourRangeGCBM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp = 1.0);

void
colourRangeRYWCB(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp = 1.0);

void
colourRangeWCB(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp = 1.0);

void
colourRangeWGCBM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp = 1.0);

void
colourRangeWYGCBM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp = 1.0);

void
colourRangeYGCBM(
  float normVal,
  float& R,
  float& G,
  float& B,
  float warp = 1.0);

/*******************************************************************************************/
/*               Generates an RGB colour map for the number range 0:1                      */
/*******************************************************************************************/

/**
 * @brief Generates an RGB colour map for the number range 0:1 via one
 * of the colourRange functions (to use in Matlab etc).
 * NOTE: lr flipped compared to BGR as opencv uses.
 *
 * @param colourRange Use the name of one of the colourRange functions
 * (eg colourRangeWCBM) available from colourmap.hpp here as
 * function parameter.
 *
 * @param sz The number of colours in the colour map.
 *
 * @param warp [0.5:2.0] control where on the 0:1 scale we want to
 * emphasize edge strength differences. 1: Neutral, lower emphasises
 * (i.e. stretches) bottom, higher emphasizes/stretches top.
 */
template<typename Lambda>
prlite::RowMatrix<float>
colourMap(
  Lambda&& colourRange,
  unsigned sz,
  float warp) {

  prlite::RowMatrix<float> cmap(sz,3);
  for (unsigned idx = 0; idx < sz; idx++) {
    float strength = float(idx)/(sz-1);
    colourRange(strength, cmap(idx,0), cmap(idx,1), cmap(idx,2), warp);
  } // for

  return cmap;
} //colourMap


/*******************************************************************************************/
/*                    Creates a uniformly coloured BGR image.                              */
/*******************************************************************************************/

/**
 * @brief creates a uniformly coloured BGR image, usefull for creating a
 * blank slate to draw on. (NOTE not RGB but flipped for OCV usage.)
 *
 * @param fillBGR: BGR colours, but on a normalized [0:1] scale
 */
std::vector< prlite::RowMatrix<float> >
filledImageBGR(
  unsigned nRows,
  unsigned nCols,
  std::vector<double> fillBGR);


/*******************************************************************************************/
/*                    Colour map a matrix.                                                 */
/*******************************************************************************************/

/**
 * @brief Transforms a matrix of doubles to a colour mapped image.
 *
 * @param pMat The information to display on a [0:1] scale.
 *
 * @param colourRange Use the name of one of the colourRange functions
 * (eg colourRangeWCBM) available from colourmap.hpp here as
 * function parameter.
 *
 * @param normalize If set the pMat will be renormalized to utilize
 * the full 0:1 scale.
 *
 * @param warp [0.5:2.0] control where on the 0:1 scale we want to
 * emphasize value differences. 1: Neutral, lower emphasises
 * (i.e. stretches) bottom, higher emphasizes/stretches top.
 */
template<typename Lambda>
std::vector< prlite::RowMatrix<float> >
colourmapImage(
  const prlite::RowMatrix<float>& pMat,
  Lambda&& colourRange,
  bool normalize = false,
  float warp = 1.0) {

  unsigned nRows = pMat.rows();
  unsigned nCols = pMat.cols();

  std::vector< prlite::RowMatrix<float> >  cmapImage(3);
  cmapImage[0].resize(nRows,nCols);
  cmapImage[1].resize(nRows,nCols);
  cmapImage[2].resize(nRows,nCols);

  // normalization parameters
  double minv = 0.0;
  double maxv = 1.0;
  if (normalize) {
    minv = pMat(0,0);
    maxv = pMat(0,0);
    for (unsigned r = 0; r < nRows; r++) {
      for (unsigned c = 0; c < nCols; c++) {
        if (pMat(r,c) > maxv) maxv = pMat(r,c);
        if (pMat(r,c) < minv) minv = pMat(r,c);
      } // for
    } // for
  } // if

  for (unsigned r = 0; r < nRows; r++) {
    for (unsigned c = 0; c < nCols; c++) {
      float strength = (pMat(r,c) - minv) / (maxv-minv);
      colourRange(strength, cmapImage[2](r,c), cmapImage[1](r,c), cmapImage[0](r,c), warp);
    } // for
  } // for

  return cmapImage;
} // colourmapImage


/*******************************************************************************************/
/*                    Overlays information on top of an existing image.                    */
/*******************************************************************************************/

/**
 * @brief Overlays information on top of an existing image.
 *
 * @param backdrop. The backdrop / original image in bgr format (each
 * channel in range [0:1]. This routine will overlay the additional
 * information in ovl on top of this image.
 *
 * @param ovl The information to overlay on a [0:1] scale. The
 * backdrop and ovl must have the same matrix size.
 *
 * @param colourRange Use the name of one of the colourRange functions
 * (eg colourRangeWCBM) available from colourmap.hpp here as
 * function parameter.
 *
 * @param normalize If set the overlay will be renormalized to utilize
 * the full 0:1 scale.
 *
 * @param thresh Only show information exceeding this value.
 *
 * @param warp [0.5:2.0] control where on the 0:1 scale we want to
 * emphasize edge strength differences. 1: Neutral, lower emphasises
 * (i.e. stretches) bottom, higher emphasizes/stretches top.
 */
template<typename Lambda>
std::vector< prlite::RowMatrix<float> >
overlayImage(
  const std::vector< prlite::RowMatrix<float> > & backdrop,
  const prlite::RowMatrix<float>& ovl,
  Lambda&& colourRange,
  bool normalize = false,
  float thresh = 0.0,
  float warp = 1.0) {

  unsigned nChans = backdrop.size();
  std::vector< prlite::RowMatrix<float> > bgr(nChans);
  for (unsigned k = 0; k < nChans; k++) {bgr[k] = backdrop[k];}
  //bgr = backdrop; // deepshallow issues

  unsigned nRows = ovl.rows();
  unsigned nCols = ovl.cols();
  PRLITE_ASSERT(bgr[0].rows() == static_cast<int>(nRows) and bgr[0].cols() == static_cast<int>(nCols),
         "backdrop image must have same dimensions as overlay");

  // normalization parameters
  double minv = 0.0;
  double maxv = 1.0;
  if (normalize) {
    minv = ovl(0,0);
    maxv = ovl(0,0);
    for (unsigned r = 0; r < nRows; r++) {
      for (unsigned c = 0; c < nCols; c++) {
        if (ovl(r,c) > maxv) maxv = ovl(r,c);
        if (ovl(r,c) < minv) minv = ovl(r,c);
      } // for
    } // for
  } // if

  for (unsigned r = 0; r < nRows; r++) {
    for (unsigned c = 0; c < nCols; c++) {
      float strength = (ovl(r,c) - minv) / (maxv-minv);
      if (strength > thresh) {
        colourRange(strength, bgr[2](r,c), bgr[1](r,c), bgr[0](r,c), warp);
      } // if
    } // for
  } // for

  return bgr;
} // overlayImage


#endif // COLOURMAP_HPP
