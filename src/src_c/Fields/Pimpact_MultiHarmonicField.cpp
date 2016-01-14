#ifdef COMPILE_ETI


#include "Pimpact_MultiHarmonicField.hpp"




template class Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > >;
template class Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > >;
template class Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > >;
template class Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > >;
template class Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > >;
template class Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > >;
template class Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > >;
template class Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > >;


#endif
