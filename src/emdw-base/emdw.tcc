template <typename T>
std::ostream& operator<<( std::ostream& file, std::set<T>& theSet){
  file << theSet.size() << "  " << std::flush;
    for (auto el: theSet) {
      file << el << " " ;
    } // for
    return file;
} // operator<<

template <typename T>
std::istream& operator>>( std::istream& file, std::set<T>& theSet ){
  unsigned sz;
  file >> sz;
  for (unsigned i = 0; i < sz; i++) {
    T tmp;
    file >> tmp;
    theSet.insert(tmp);
  } // for
  return file;
} // operator>>

#define FILEPOS std::cout << __FILE__ << __LINE__ << std::endl
