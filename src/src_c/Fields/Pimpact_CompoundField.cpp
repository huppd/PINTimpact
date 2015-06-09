#ifdef COMPILE_ETI

#include "Pimpact_CompoundField.hpp"

template class Pimpact::CompoundField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> >, Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > >;
template class Pimpact::CompoundField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> >, Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > >;
template class Pimpact::CompoundField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > >, Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > > >;
template class Pimpact::CompoundField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > >, Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > > >;
template class Pimpact::CompoundField< Pimpact::ModeField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > >, Pimpact::ModeField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > > >;
template class Pimpact::CompoundField< Pimpact::ModeField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > >, Pimpact::ModeField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > > >;
template class Pimpact::CompoundField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > >, Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > > >;
template class Pimpact::CompoundField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > >, Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > > >;

#endif
