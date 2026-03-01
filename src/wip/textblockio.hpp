#ifndef TEXTBLOCKIO_HPP
#define TEXTBLOCKIO_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    March 2020
          PURPOSE: Block IO for numeric matrices

*******************************************************************************/

#include <string>
#include "gLinear/Matrix.h"

/******************************************************************************/

/**
 * This loads a blocked/matrix file of numbers
 */
template <typename AssignType>
gLinear::gRowMatrix<AssignType> loadBlock(const std::string& fileName);

/**
 * This saves a blocked/matrix file of numbers
 */
template <typename AssignType>
void saveBlock(const std::string& fileName, const gLinear::gRowMatrix<AssignType>& dataBuffer);

#include "textblockio.tcc"

#endif // TEXTBLOCKIO_HPP
