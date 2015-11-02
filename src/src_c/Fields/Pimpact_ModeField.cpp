#ifdef COMPILE_ETI


#include "Pimpact_ModeField.hpp"




template class Pimpact::ModeField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > >;
template class Pimpact::ModeField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > >;
template class Pimpact::ModeField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > >;
template class Pimpact::ModeField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > >;


#endif
