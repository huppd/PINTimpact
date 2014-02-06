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


/**
 * \file
 * \todo make generat templated create function
 */
namespace Pimpact {


template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<VectorField<Scalar,Ordinal> > > > createMultiModeVectorField( int n=1 ) {
	auto fS  = createFieldSpace<Ordinal>();
	auto iIS = createInnerFieldIndexSpaces<Ordinal>();
	auto fIS = createFullFieldIndexSpaces<Ordinal>();

	auto velc = createVectorField<Scalar,Ordinal>(fS,iIS,fIS);
	auto vels = createVectorField<Scalar,Ordinal>(fS,iIS,fIS);

	auto vel = createModeField( velc, vels );

	return( createMultiField<ModeField<VectorField<Scalar,Ordinal> >,double,int>( *vel, n ) );
}


template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<ScalarField<Scalar,Ordinal> > > > createMultiModeScalarField( int n=1 ) {
	auto sVS = createFieldSpace<Ordinal>();

	auto scalc = createScalarField<Scalar,Ordinal>(sVS);
	auto scals = createScalarField<Scalar,Ordinal>(sVS);

	auto scal = createModeField( scalc, scals );

	return( createMultiField<ModeField<ScalarField<Scalar,Ordinal> >,Scalar,Ordinal>( *scal, n ) );
}


template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<VectorField<Scalar,Ordinal> > > > createInitMVF(
		EFlowType flowType,
		Teuchos::RCP<const FieldSpace<Ordinal> > fS,
		Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<Ordinal> > >iIS,
		Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<Ordinal> > >fIS,
		Scalar re=1., Scalar omega=1., Scalar px=1.) {

	auto velc = createVectorField<double,int>(fS,iIS,fIS);
	auto vels = createVectorField<double,int>(fS,iIS,fIS);

	auto vel  = createModeField( velc, vels );

	auto u   = Pimpact::createMultiField< ModeField<VectorField<Scalar,Ordinal> > , double, int >(*vel,1);

	switch( EFlowType( flowType ) ) {
	case Pimpact::ZeroFLow :
		u->GetVec(0).getFieldC()->init_field( ZeroProf );
		u->GetVec(0).getFieldS()->init_field( ZeroProf );
		break;
	case Pimpact::Poiseuille_inX :
		u->GetVec(0).getFieldC()->init_field( Poiseuille2D_inX );
		u->GetVec(0).getFieldS()->init_field( ZeroProf );
		break;
	case Poiseuille_inY :
		u->GetVec(0).getFieldC()->init_field( Poiseuille2D_inY );
		u->GetVec(0).getFieldS()->init_field( ZeroProf );
		break;
	case Pulsatile_inX :
		u->GetVec(0).getFieldC()->init_field( Pulsatile2D_inXC, re, omega, px );
		u->GetVec(0).getFieldS()->init_field( Pulsatile2D_inXS, re, omega, px );
		break;
	case Pulsatile_inY :
		u->GetVec(0).getFieldC()->init_field( Pulsatile2D_inYC, re, omega, px );
		u->GetVec(0).getFieldS()->init_field( Pulsatile2D_inYS, re, omega, px );
		break;
	case Streaming2DFlow :
		u->GetVec(0).getFieldC()->init_field( Streaming2D, re, omega, px );
		u->GetVec(0).getFieldS()->init_field( ZeroProf, re, omega, px );
		break;
	}
	return( u );
}


template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<ScalarField<Scalar,Ordinal> > > > createInitMSF(
		Teuchos::RCP< const FieldSpace<Ordinal> > fS ) {
	auto scac = Pimpact::createScalarField<Scalar,Ordinal>(fS);
	auto scas = Pimpact::createScalarField<Scalar,Ordinal>(fS);

	auto sca = Pimpact::createModeField( scac, scas );

	return( Pimpact::createMultiField<ModeField<ScalarField<Scalar,Ordinal> >, double, int >(*sca,1) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_FIELDFACTORY_HPP
