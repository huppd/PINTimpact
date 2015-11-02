#include "NOX_Pimpact_Vector.hpp"




#ifdef COMPILE_ETI

template class NOX::Pimpact::Vector<
Pimpact::MultiField<Pimpact::VectorField<Pimpact::Space<double,int,3,2>
> > >; 
template class NOX::Pimpact::Vector< Pimpact::CompoundField<
Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<Pimpact::Space<double,int,3,4>
> > >,
	Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<Pimpact::Space<double,int,3,4>
	> > > > >; 
template class NOX::Pimpact::Vector<
		Pimpact::MultiField<
			Pimpact::CompoundField<
				Pimpact::MultiHarmonicField<Pimpact::VectorField<Pimpact::Space<double,int,3,4> > >,
				Pimpact::MultiHarmonicField<Pimpact::ScalarField<Pimpact::Space<double,int,3,4> > > 
			>
		>
	>; 
template class NOX::Pimpact::Vector<
		Pimpact::MultiField<
			Pimpact::CompoundField<
				Pimpact::TimeField<Pimpact::VectorField<Pimpact::Space<double,int,4,4> > >,
				Pimpact::TimeField<Pimpact::ScalarField<Pimpact::Space<double,int,4,4> > > 
			>
		>
	>; 
#endif
