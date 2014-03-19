#pragma once
#ifndef PIMPACT_FIELDFACTORY_HPP
#define PIMPACT_FIELDFACTORY_HPP

//#include "Teuchos_RCP.hpp"

#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"


/// \file
/// \todo make general templated create function

/// \defgroup Field Fields
///
/// back bone of Pimpact,  can be seen as linerar algebra vectors

namespace Pimpact {


template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<VectorField<Scalar,Ordinal> > > > createMultiModeVectorField( int n=1 ) {
	auto fS  = createFieldSpace<Ordinal>();
	auto iIS = createInnerFieldIndexSpaces<Ordinal>();
	auto fIS = createFullFieldIndexSpaces<Ordinal>();

	auto velc = createVectorField<Scalar,Ordinal>(fS,iIS,fIS);
	auto vels = createVectorField<Scalar,Ordinal>(fS,iIS,fIS);

	auto vel = createModeField( velc, vels );

	return( createMultiField<ModeField<VectorField<Scalar,Ordinal> > >( *vel, n ) );
}


template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<ScalarField<Scalar,Ordinal> > > > createMultiModeScalarField( int n=1 ) {
	auto sVS = createFieldSpace<Ordinal>();

	auto scalc = createScalarField<Scalar,Ordinal>(sVS);
	auto scals = createScalarField<Scalar,Ordinal>(sVS);

	auto scal = createModeField( scalc, scals );

	return( createMultiField<ModeField<ScalarField<Scalar,Ordinal> > >( *scal, n ) );
}


template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<VectorField<Scalar,Ordinal> > > > createInitMVF(
		EFlowType flowType,
		const Teuchos::RCP<const FieldSpace<Ordinal> >& fS,
		const Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<Ordinal> > >& iIS,
		const Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<Ordinal> > >& fIS,
		Scalar re=1., Scalar omega=1., Scalar px=1.) {

	auto velc = createVectorField<double,int>(fS,iIS,fIS);
	auto vels = createVectorField<double,int>(fS,iIS,fIS);

	auto vel  = createModeField( velc, vels );

	auto u   = Pimpact::createMultiField< ModeField<VectorField<Scalar,Ordinal> > >(*vel,1);

	switch( EFlowType( flowType ) ) {
	case Pimpact::Zero2DFlow :
		u->getField(0).getCFieldPtr()->initField( ZeroProf );
		u->getField(0).getSFieldPtr()->initField( ZeroProf );
		break;
	case Pimpact::Poiseuille_inX :
		u->getField(0).getCFieldPtr()->initField( Poiseuille2D_inX );
		u->getField(0).getSFieldPtr()->initField( ZeroProf );
		break;
	case Poiseuille_inY :
		u->getField(0).getCFieldPtr()->initField( Poiseuille2D_inY );
		u->getField(0).getSFieldPtr()->initField( ZeroProf );
		break;
	case Pulsatile_inX :
		u->getField(0).getCFieldPtr()->initField( Pulsatile2D_inXC, re, omega, px );
		u->getField(0).getSFieldPtr()->initField( Pulsatile2D_inXS, re, omega, px );
		break;
	case Pulsatile_inY :
		u->getField(0).getCFieldPtr()->initField( Pulsatile2D_inYC, re, omega, px );
		u->getField(0).getSFieldPtr()->initField( Pulsatile2D_inYS, re, omega, px );
		break;
	case Streaming2DFlow :
		u->getField(0).getCFieldPtr()->initField( Streaming2D, re, omega, px );
		u->getField(0).getSFieldPtr()->initField( ZeroProf, re, omega, px );
		break;
	}
	return( u );
}


template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<ScalarField<Scalar,Ordinal> > > > createInitMSF(
		const Teuchos::RCP< const FieldSpace<Ordinal> >& fS ) {
	auto scac = Pimpact::createScalarField<Scalar,Ordinal>(fS);
	auto scas = Pimpact::createScalarField<Scalar,Ordinal>(fS);

	auto sca = Pimpact::createModeField( scac, scas );

	return( Pimpact::createMultiField<ModeField<ScalarField<Scalar,Ordinal> > >(*sca,1) );
}



/// \brief creates a multi-harmonic scalar field.
///
/// \param fS scalar Vector Space to which returned vector belongs
/// \param nf amount of modes
/// \return field vector
template<class S, class O>
Teuchos::RCP< MultiHarmonicField< ScalarField<S,O> > > createMultiHarmonicScalarField(
    const Teuchos::RCP<const FieldSpace<O> >& fS, int nf
    ) {
  auto field0 = createScalarField<S,O>( fS );
  auto mfield = createModeField< ScalarField<S,O> >( field0, field0 );
  auto fields = createMultiField< ModeField< ScalarField<S,O> > >( *mfield, nf );
  return Teuchos::rcp(
        new MultiHarmonicField< ScalarField<S,O> >( field0, fields ) );
}



/// \brief creates a multi-harmonic vector field.
///
/// \param fS scalar Vector Space to which returned vector belongs
/// \param nf amount of modes
/// \return field vector
template<class S, class O>
Teuchos::RCP< MultiHarmonicField< VectorField<S,O> > > createMultiHarmonicVectorField(
    const Teuchos::RCP<const FieldSpace<O> >& fieldS,
    typename VectorField<S,O>::IndexSpaces innerIS,
    const typename VectorField<S,O>::IndexSpaces& fullIS,
    int nf ) {
  auto field0 = createVectorField<S,O>( fieldS, innerIS, fullIS );
  auto mfield = createModeField< VectorField<S,O> >( field0, field0 );
  auto fields = createMultiField< ModeField< VectorField<S,O> > >( *mfield, nf );
  return Teuchos::rcp(
        new MultiHarmonicField< VectorField<S,O> >( field0, fields ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_FIELDFACTORY_HPP
