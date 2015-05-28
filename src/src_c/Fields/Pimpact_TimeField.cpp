#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_TimeField.hpp"



namespace Pimpact {

	template class TimeField<ScalarField< Space<double,int,4,2> > >;
	template class TimeField<ScalarField< Space<double,int,4,4> > >;
	template class TimeField<VectorField< Space<double,int,4,2> > >;
	template class TimeField<VectorField< Space<double,int,4,4> > >;

} // end of namespace Pimpact
