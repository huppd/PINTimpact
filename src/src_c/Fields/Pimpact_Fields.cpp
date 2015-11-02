#include "Pimpact_Fields.hpp"




#ifdef COMPILE_ETI
namespace Pimpact {

	// MultiFields
	// ScalarFields
	template class MultiField< ScalarField< Space<double,int,3,2> > >;
	template class MultiField< ScalarField< Space<double,int,3,4> > >;
	template class MultiField< ScalarField< Space<double,int,4,2> > >;
	template class MultiField< ScalarField< Space<double,int,4,4> > >;

	// VectorFields
	template class MultiField< VectorField< Space<double,int,3,2> > >;
	template class MultiField< VectorField< Space<double,int,3,4> > >;
	template class MultiField< VectorField< Space<double,int,4,2> > >;
	template class MultiField< VectorField< Space<double,int,4,4> > >;

	// TimeFields
	template class MultiField< TimeField<ScalarField< Space<double,int,4,2> > > >;
	template class MultiField< TimeField<ScalarField< Space<double,int,4,4> > > >;
	template class MultiField< TimeField<VectorField< Space<double,int,4,2> > > >;
	template class MultiField< TimeField<VectorField< Space<double,int,4,4> > > >;

	// ModeFields
	template class MultiField< ModeField< ScalarField< Space<double,int,3,2> > > >;
	template class MultiField< ModeField< ScalarField< Space<double,int,3,4> > > >;
	template class MultiField< ModeField< VectorField< Space<double,int,3,2> > > >;
	template class MultiField< ModeField< VectorField< Space<double,int,3,4> > > >;

	// MultiHarmonicFieldsMultiField< 
	template class MultiField< MultiHarmonicField< ScalarField< Space<double,int,3,2> > > >;
	template class MultiField< MultiHarmonicField< ScalarField< Space<double,int,3,4> > > >;
	template class MultiField< MultiHarmonicField< VectorField< Space<double,int,3,2> > > >;
	template class MultiField< MultiHarmonicField< VectorField< Space<double,int,3,4> > > >;

	// CompoundFields
	template class MultiField< CompoundField< VectorField< Space<double,int,3,2> >, ScalarField< Space<double,int,3,2> > > >;
	template class MultiField< CompoundField< VectorField< Space<double,int,3,4> >, ScalarField< Space<double,int,3,4> > > >;
	template class MultiField< CompoundField< TimeField< VectorField< Space<double,int,4,2> > >, TimeField< ScalarField< Space<double,int,4,2> > > > >;
	template class MultiField< CompoundField< TimeField< VectorField< Space<double,int,4,4> > >, TimeField< ScalarField< Space<double,int,4,4> > > > >;
	template class MultiField< CompoundField< ModeField< VectorField< Space<double,int,3,2> > >, ModeField< ScalarField< Space<double,int,3,2> > > > >;
	template class MultiField< CompoundField< ModeField< VectorField< Space<double,int,3,4> > >, ModeField< ScalarField< Space<double,int,3,4> > > > >;
	template class MultiField< CompoundField< MultiHarmonicField< VectorField< Space<double,int,3,2> > >, MultiHarmonicField< ScalarField< Space<double,int,3,2> > > > >;
	template class MultiField< CompoundField< MultiHarmonicField< VectorField< Space<double,int,3,4> > >, MultiHarmonicField< ScalarField< Space<double,int,3,4> > > > >;
} // end of namespace Pimpact


#endif
