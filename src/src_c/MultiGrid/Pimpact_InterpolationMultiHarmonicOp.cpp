#include "Pimpact_InterpolationMultiHarmonicOp.hpp"



#ifdef COMPILE_ETI
template class Pimpact::InterpolationMultiHarmonicOp< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,2> > >;
template class Pimpact::InterpolationMultiHarmonicOp< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,4> > >;
template class Pimpact::InterpolationMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,2> > > >;
template class Pimpact::InterpolationMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,4> > > >;
#endif
