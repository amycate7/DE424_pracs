#ifndef CONDITIONALDIRICHLET_HPP
#define CONDITIONALDIRICHLET_HPP

// To support the old deprecated name

#include "dirichletset.hpp"

// deprecated name for a DirichletSet
template <typename CondType, typename DirType>
using ConditionalDirichlet = DirichletSet<CondType,DirType>;

#endif // CONDITIONALDIRICHLET_HPP
