#pragma once
#ifndef PIMPACT_FIELDFACTORY_HPP
#define PIMPACT_FIELDFACTORY_HPP

//#include "Teuchos_RCP.hpp"

#include "Pimpact_Space.hpp"
//#include "Pimpact_IndexSpace.hpp"

//#include "Pimpact_ScalarField.hpp"
//#include "Pimpact_VectorField.hpp"
//#include "Pimpact_ModeField.hpp"
//#include "Pimpact_MultiField.hpp"
//#include "Pimpact_MultiHarmonicField.hpp"
#include "Pimpact_Fields.hpp"



/// \file
/// \todo make general templated create function



namespace Pimpact {


/// \relates MultiField
/// \relates ModeField
/// \relates VectorField
template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<VectorField<Scalar,Ordinal> > > > createMultiModeVectorField( int n=1 ) {

	auto space = Pimpact::createSpace();

  auto velc = createVectorField<Scalar,Ordinal>( space );
  auto vels = createVectorField<Scalar,Ordinal>( space );

  auto vel = createModeField( velc, vels );

  return( createMultiField<ModeField<VectorField<Scalar,Ordinal> > >( *vel, n ) );
}



/// \relates MultiField
/// \relates ModeField
/// \relates ScalarField
template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<ScalarField<Scalar,Ordinal> > > > createMultiModeScalarField( int n=1 ) {
	auto sVS = Pimpact::createSpace<Ordinal>();

  auto scalc = createScalarField<Scalar,Ordinal>(sVS);
  auto scals = createScalarField<Scalar,Ordinal>(sVS);

  auto scal = createModeField( scalc, scals );

  return( createMultiField<ModeField<ScalarField<Scalar,Ordinal> > >( *scal, n ) );
}



/// \relates MultiField
/// \relates ModeField
/// \relates VectorField
template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<VectorField<Scalar,Ordinal> > > > createInitMVF(
    EFlowType flowType,
    const Teuchos::RCP<const Space<Ordinal> >& space,
//    const Teuchos::RCP<const FieldSpace<Ordinal> >& fS,
//    const Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<Ordinal> > >& iIS,
//    const Teuchos::ArrayRCP< Teuchos::RCP< const IndexSpace<Ordinal> > >& fIS,
    Scalar re=1., Scalar omega=1., Scalar px=1.) {

  auto velc = createVectorField<double,int>(space);
  auto vels = createVectorField<double,int>(space);

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
  case Streaming2DFlow2 :
    u->getField(0).getCFieldPtr()->initField( Streaming2DC, re, omega, px );
    u->getField(0).getSFieldPtr()->initField( Streaming2DS, re, omega, px );
    break;
  case Streaming2DFlow3 :
    u->getField(0).getCFieldPtr()->initField( Streaming2DC, re, omega, px );
    u->getField(0).getSFieldPtr()->initField( Streaming2DS, re, omega, px );
    break;
  }
  return( u );
}



/// \relates MultiField
/// \relates ModeField
/// \relates ScalarField
template<class Scalar, class Ordinal>
Teuchos::RCP< MultiField<ModeField<ScalarField<Scalar,Ordinal> > > > createInitMSF(
    const Teuchos::RCP< const Space<Ordinal> > space ) {
//    const Teuchos::RCP< const FieldSpace<Ordinal> >& fS,
//    const Teuchos::RCP< const IndexSpace<Ordinal> >& iS ) {

  auto scac = Pimpact::createScalarField<Scalar,Ordinal>(space);
  auto scas = Pimpact::createScalarField<Scalar,Ordinal>(space);

  auto sca = Pimpact::createModeField( scac, scas );

  return( Pimpact::createMultiField<ModeField<ScalarField<Scalar,Ordinal> > >(*sca,1) );
}



/// \brief creates a multi-harmonic scalar field.
///
/// \relates MultiHarmonicField
/// \relates ScalarField
/// \param fS scalar Vector Space to which returned vector belongs
/// \param nf amount of modes
/// \return field vector
template<class S, class O>
Teuchos::RCP< MultiHarmonicField< ScalarField<S,O> > > createMultiHarmonicScalarField(
    const Teuchos::RCP<const Space<O> >& space,
    int nf) {
  auto field0 = createScalarField<S,O>( space );
  auto mfield = createModeField< ScalarField<S,O> >( field0, field0 );
  auto fields = createMultiField< ModeField< ScalarField<S,O> > >( *mfield, nf );
  return(
      Teuchos::rcp(
          new MultiHarmonicField< ScalarField<S,O> >( field0, fields ) ) );
}



/// \brief creates a multi-harmonic vector field.
///
/// \relates MultiHarmonicField
/// \relates VectorField
/// \param fieldS scalar Vector Space to which returned vector belongs
/// \param innerIS
/// \param fullIS
/// \param nf amount of modes
/// \return field vector
template<class S, class O>
Teuchos::RCP< MultiHarmonicField< VectorField<S,O> > >
createMultiHarmonicVectorField(
    const Teuchos::RCP< const Space<O> >& space,
    int nf ) {
  auto field0 = createVectorField<S,O>( space );
  auto mfield = createModeField< VectorField<S,O> >( field0, field0 );
  auto fields = createMultiField< ModeField< VectorField<S,O> > >( *mfield, nf );
  return Teuchos::rcp(
      new MultiHarmonicField< VectorField<S,O> >( field0, fields ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_FIELDFACTORY_HPP
