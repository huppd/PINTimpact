#include "Pimpact_TransferMultiHarmonicOp.hpp"




#ifdef COMPILE_ETI
template class Pimpact::TransferMultiHarmonicOp< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,2> > >;
template class Pimpact::TransferMultiHarmonicOp< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,4> > >;
template class Pimpact::TransferMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,2> > > >;
template class Pimpact::TransferMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp< Pimpact::Space<double,int,3,4> > > >;
#endif
