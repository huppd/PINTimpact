#ifdef COMPILE_ETI

#include "BelosPimpactAdapter.hpp"

// ScalarFields
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > > >;

// VectorFields
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > > >;

// MultiHarmonicFieldsBelos
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > > > >;

// TimeFields 
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > > > >;

// CompoundFields
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::CompoundField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > >, Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > > > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::CompoundField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > >, Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > > > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::CompoundField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > >, Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > > > > >;
template class Belos::MultiVecTraits< double, Pimpact::MultiField< Pimpact::CompoundField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > >, Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > > > > >;


// ScalarFields
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > > > >;

// VectorFields
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > > > >;

// MultiHarmonicFieldsBelos::
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > > > > >;

// TimeFields 
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > > > > >;

//// CompoundFields
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::CompoundField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > >, Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > > > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::CompoundField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,2> > >, Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,2> > > > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::CompoundField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > >, Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > > > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::CompoundField< Pimpact::MultiHarmonicField< Pimpact::VectorField< Pimpact::Space<double,int,3,4> > >, Pimpact::MultiHarmonicField< Pimpact::ScalarField< Pimpact::Space<double,int,3,4> > > > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::CompoundField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > >, Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > > > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::CompoundField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,2> > >, Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,2> > > > > > >;
template class Belos::OperatorTraits< double, Pimpact::MultiField< Pimpact::CompoundField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > >, Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > > > >, Pimpact::OperatorBase< Pimpact::MultiField< Pimpact::CompoundField< Pimpact::TimeField< Pimpact::VectorField< Pimpact::Space<double,int,4,4> > >, Pimpact::TimeField< Pimpact::ScalarField< Pimpact::Space<double,int,4,4> > > > > > >;
#endif
