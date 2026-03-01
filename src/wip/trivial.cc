class Trivial {
public:
  Trivial();
  Trivial(const Trivial& triv);
  ~Trivial();
  Trivial& operator=(const Trivial& p);
  Trivial* copy() const;

  static unsigned noOfObjects;
}; // Trivial

unsigned Trivial::noOfObjects = 0u;

Trivial::Trivial() {
  ++noOfObjects;
  std::cout << '<' << noOfObjects << '>' << std::flush;
} // def ctor

Trivial::Trivial(const Trivial&) {
  ++noOfObjects;
  std::cout << '<' << noOfObjects << '>' << std::flush;
} // cpy ctor

Trivial::~Trivial() {
  --noOfObjects;
  std::cout << '<' << noOfObjects << '>' << std::flush;
} // cpy ctor

Trivial& Trivial::operator=(const Trivial& p) {
  if (this != &p) {
    this->~Trivial();              // Destroy existing
    new ( static_cast<void *>(this) ) Trivial(p); // do an inplace construction
  } // if

  return *this;
} // operator=

Trivial* Trivial::copy() const {
  return new Trivial(*this);
} // copy
