#pragma once
#ifndef PIMPACT_FIELDS_HPP
#define PIMPACT_FIELDS_HPP


/// \defgroup Field Fields
///
/// back bone of Pimpact,  can be seen as linerar algebra vectors

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_FieldFactory.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_TimeField.hpp"
#include "Pimpact_VectorField.hpp"




#ifdef COMPILE_ETI
namespace Pimpact {

	// MultiFields
	// ScalarFields
	extern template class MultiField< ScalarField< Space<double,int,3,2> > >;
	extern template class MultiField< ScalarField< Space<double,int,3,4> > >;
	extern template class MultiField< ScalarField< Space<double,int,4,2> > >;
	extern template class MultiField< ScalarField< Space<double,int,4,4> > >;

	// VectorFields
	extern template class MultiField< VectorField< Space<double,int,3,2> > >;
	extern template class MultiField< VectorField< Space<double,int,3,4> > >;
	extern template class MultiField< VectorField< Space<double,int,4,2> > >;
	extern template class MultiField< VectorField< Space<double,int,4,4> > >;

	// TimeFields
	extern template class MultiField< TimeField<ScalarField< Space<double,int,4,2> > > >;
	extern template class MultiField< TimeField<ScalarField< Space<double,int,4,4> > > >;
	extern template class MultiField< TimeField<VectorField< Space<double,int,4,2> > > >;
	extern template class MultiField< TimeField<VectorField< Space<double,int,4,4> > > >;

	// ModeFields
	extern template class MultiField< ModeField< ScalarField< Space<double,int,3,2> > > >;
	extern template class MultiField< ModeField< ScalarField< Space<double,int,3,4> > > >;
	extern template class MultiField< ModeField< VectorField< Space<double,int,3,2> > > >;
	extern template class MultiField< ModeField< VectorField< Space<double,int,3,4> > > >;

	// MultiHarmonicFieldsMultiField< 
	extern template class MultiField< MultiHarmonicField< ScalarField< Space<double,int,3,2> > > >;
	extern template class MultiField< MultiHarmonicField< ScalarField< Space<double,int,3,4> > > >;
	extern template class MultiField< MultiHarmonicField< VectorField< Space<double,int,3,2> > > >;
	extern template class MultiField< MultiHarmonicField< VectorField< Space<double,int,3,4> > > >;

	// CompoundFields
	extern template class MultiField< CompoundField< VectorField< Space<double,int,3,2> >, ScalarField< Space<double,int,3,2> > > >;
	extern template class MultiField< CompoundField< VectorField< Space<double,int,3,4> >, ScalarField< Space<double,int,3,4> > > >;
	extern template class MultiField< CompoundField< TimeField< VectorField< Space<double,int,4,2> > >, TimeField< ScalarField< Space<double,int,4,2> > > > >;
	extern template class MultiField< CompoundField< TimeField< VectorField< Space<double,int,4,4> > >, TimeField< ScalarField< Space<double,int,4,4> > > > >;
	extern template class MultiField< CompoundField< ModeField< VectorField< Space<double,int,3,2> > >, ModeField< ScalarField< Space<double,int,3,2> > > > >;
	extern template class MultiField< CompoundField< ModeField< VectorField< Space<double,int,3,4> > >, ModeField< ScalarField< Space<double,int,3,4> > > > >;
	extern template class MultiField< CompoundField< MultiHarmonicField< VectorField< Space<double,int,3,2> > >, MultiHarmonicField< ScalarField< Space<double,int,3,2> > > > >;
	extern template class MultiField< CompoundField< MultiHarmonicField< VectorField< Space<double,int,3,4> > >, MultiHarmonicField< ScalarField< Space<double,int,3,4> > > > >;

} // end of namespace Pimpact
#endif


#endif // end of #ifndef PIMPACT_FIELDS_HPP
