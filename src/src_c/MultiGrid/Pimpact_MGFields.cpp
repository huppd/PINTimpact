#include "Pimpact_MGFields.hpp"

#ifdef COMPILE_ETI
template class Pimpact::MGFields< Pimpact::MGSpaces< Pimpact::Space<double,int,3,4>, Pimpact::Space<double,int,3,2> >, Pimpact::ScalarField >;
template class Pimpact::MGFields< Pimpact::MGSpaces< Pimpact::Space<double,int,4,4>, Pimpact::Space<double,int,4,2> >, Pimpact::ScalarField >;
template class Pimpact::MGFields< Pimpact::MGSpaces< Pimpact::Space<double,int,3,4>, Pimpact::Space<double,int,3,2> >, Pimpact::VectorField >;
template class Pimpact::MGFields< Pimpact::MGSpaces< Pimpact::Space<double,int,4,4>, Pimpact::Space<double,int,4,2> >, Pimpact::VectorField >;
#endif
