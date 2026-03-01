/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    March 2020
          PURPOSE: Block IO for numeric matrices

*******************************************************************************/

#include <iostream>
#include <sstream>
#include "prlite_error.hpp"

#include "prlite_zfstream.hpp"  //gzifstream gzofstream

template <typename AssignType>
gLinear::gRowMatrix<AssignType> loadBlock(const std::string& fileName){
  prlite::gzifstream fileIn(fileName.c_str(), std::ios::in);
  PRLITE_ASSERT(fileIn, "TextBlockFeatures::findFileDimensions : cannot open file '" << fileName << "'");

  size_t fileVecs = 0;
  size_t fileDim = 0;

  // ***************************************
  // * Determine the size of the file
  // ***************************************

  while ( !fileIn.eof() ) {

    std::string s;
    getline(fileIn, s);
    std::stringstream ist(s);
    float tmpFloat = 0;
    ist >> tmpFloat;
    // there has to be something on the line
    if (!ist) {
      if ( !fileIn.eof() ) {
        std::cout << "loadBlock: trailing lines in file " << fileName << std::endl;
      } // if
      break;
    } // if

    if (!fileVecs) {  // first line
      fileDim++;
      while (ist >> tmpFloat) {
        fileDim++;
      } // while
    } // if
    fileVecs++;

  } // while
  PRLITE_ASSERT(fileVecs > 0 && fileDim > 0, "TextBlockFeatures::findFileDimensions : the file " << fileName << " seems to be empty");

  // ******************************
  // * now do the actual loading
  // ******************************

  fileIn.close();
  fileIn.open(fileName.c_str(), std::ios::in);

  gLinear::gRowMatrix<AssignType> dataBuffer(fileVecs, fileDim);

  // seek to begin pos
  fileVecs = 0;

  int nrRead = 0;
  while ( !fileIn.eof() ) {
    std::string s;
    getline(fileIn, s);
    std::stringstream ist(s);
    int lineDim = 0;
    while ( ist >> dataBuffer(nrRead,lineDim) ) {
      lineDim++;
    } // while
    nrRead++; fileVecs++;
  } // while
  fileIn.close();

  return dataBuffer;
} // loadBlock

template <typename AssignType>
void saveBlock(const std::string& fileName, const gLinear::gRowMatrix<AssignType>& dataBuffer) {
  prlite::gzofstream file;
  file.open( fileName.c_str(), static_cast<std::ios::openmode>(std::ios::out) );
  PRLITE_ASSERT(file, "Could not open file '" << fileName << "'");
  //file << setw(precis + 2) << setprecision(precis);

  int nRows = dataBuffer.rows();
  int nCols = dataBuffer.cols();
  for (int i = 0; i < nRows; i++) {
    int cnt = 10;
    for (int k = 0; k < nCols; k++) {
      file << dataBuffer(i, k) << ' ';
      if (!--cnt) {
        file << "  ";
        cnt = 10;
      } // if
    } // for k
    file << std::endl;
  } // for i
  file.close();

} // saveBlock
