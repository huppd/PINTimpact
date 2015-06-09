#ifdef COMPILE_ETI

#include "Pimpact_TimeField.hpp"


template class Pimpact::TimeField<Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > >;
template class Pimpact::TimeField<Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > >;
template class Pimpact::TimeField<Pimpact::VectorField< Pimpact::Space<double,int,4,2> > >;
template class Pimpact::TimeField<Pimpact::VectorField< Pimpact::Space<double,int,4,4> > >;

#endif
