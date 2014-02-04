#pragma once
#ifndef PIMPACT_FIELDFACTORY_HPP
#define PIMPACT_FIELDFACTORY_HPP

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiField.hpp"


namespace Pimpact {


template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<VectorField<Scalar,Ordinal> > > > createMultiModeVectorField( int n=1 ) {
	auto fS = Pimpact::createFieldSpace<Ordinal>();
	auto iIS = Pimpact::createInnerFieldIndexSpaces<Ordinal>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<Ordinal>();

	auto velc = Pimpact::createVectorField<Scalar,Ordinal>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<Scalar,Ordinal>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );

	return Pimpact::createMultiField<Pimpact::ModeField<Pimpact::VectorField<Scalar,Ordinal> >,double,int>( *vel, n );
}


template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<ScalarField<Scalar,Ordinal> > > > createMultiModeScalarField( int n=1 ) {
	auto sVS = Pimpact::createFieldSpace<Ordinal>();

	auto scalc = Pimpact::createScalarField<Scalar,Ordinal>(sVS);
	auto scals = Pimpact::createScalarField<Scalar,Ordinal>(sVS);

	auto scal = Pimpact::createModeField( scalc, scals );

	return Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<Scalar,Ordinal> >,Scalar,Ordinal>( *scal, n );
}


template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<VectorField<Scalar,Ordinal> > > > createInitMVF(
		EFlowType flowType, Teuchos::RCP<FieldSpace<Ordinal> > fS, )

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_FIELDFACTORY_HPP
