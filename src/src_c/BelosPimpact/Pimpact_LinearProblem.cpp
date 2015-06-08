#include "Pimpact_Fields.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_LinearProblem.hpp"


#ifdef COMPILE_ETI
namespace Pimpact {

	template class LinearProblem< MultiField<ScalarField< Space<double,int,3,2> > > >;
	template class LinearProblem< MultiField<ScalarField< Space<double,int,3,4> > > >;
	template class LinearProblem< MultiField<ScalarField< Space<double,int,4,2> > > >;
	template class LinearProblem< MultiField<ScalarField< Space<double,int,4,4> > > >;

	// VectorFields
	template class LinearProblem< MultiField<VectorField< Space<double,int,3,2> > > >;
	template class LinearProblem< MultiField<VectorField< Space<double,int,3,4> > > >;
	template class LinearProblem< MultiField<VectorField< Space<double,int,4,2> > > >;
	template class LinearProblem< MultiField<VectorField< Space<double,int,4,4> > > >;

	// MultiHarmonicFields
	template class LinearProblem< MultiField <MultiHarmonicField< ScalarField< Space<double,int,3,2> > > > >;
	template class LinearProblem< MultiField <MultiHarmonicField< ScalarField< Space<double,int,3,4> > > > >;

	template class LinearProblem< MultiField <MultiHarmonicField< VectorField< Space<double,int,3,2> > > > >;
	template class LinearProblem< MultiField <MultiHarmonicField< VectorField< Space<double,int,3,4> > > > >;

	// TimeFields
	template class LinearProblem< MultiField <TimeField< ScalarField< Space<double,int,4,2> > > > >;
	template class LinearProblem< MultiField <TimeField< ScalarField< Space<double,int,4,4> > > > >;

	template class LinearProblem< MultiField <TimeField< VectorField< Space<double,int,4,2> > > > >;
	template class LinearProblem< MultiField <TimeField< VectorField< Space<double,int,4,4> > > > >;

	// CompoundFields
	template class LinearProblem< MultiField< CompoundField< MultiHarmonicField< VectorField< Space<double,int,3,2> > >, MultiHarmonicField< ScalarField< Space<double,int,3,2> > > > > >;
	template class LinearProblem< MultiField< CompoundField< MultiHarmonicField< VectorField< Space<double,int,3,4> > >, MultiHarmonicField< ScalarField< Space<double,int,3,4> > > > > >;
	template class LinearProblem< MultiField< CompoundField< TimeField< VectorField< Space<double,int,4,2> > >, TimeField< ScalarField< Space<double,int,4,2> > > > > >;
	template class LinearProblem< MultiField< CompoundField< TimeField< VectorField< Space<double,int,4,4> > >, TimeField< ScalarField< Space<double,int,4,4> > > > > >;


} // end of namespace Pimpact
#endif
