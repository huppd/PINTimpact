#pragma once
#ifndef PIMPACT_FIELDFACTORY_HPP
#define PIMPACT_FIELDFACTORY_HPP


#include "Pimpact_Space.hpp"

#include "Pimpact_Fields.hpp"



/// \file


namespace Pimpact {


/// \relates MultiField
/// \relates ModeField
/// \relates VectorField
template<class SpaceT>
Teuchos::RCP< MultiField<ModeField<VectorField<SpaceT> > > >
createMultiModeVectorField( const Teuchos::RCP<const SpaceT>& space, int n=1 ) {

  auto vel = create< ModeField< VectorField<SpaceT> > >( space );

  return( createMultiField<ModeField<VectorField<SpaceT> > >( *vel, n ) );

}



/// \relates MultiField
/// \relates ModeField
/// \relates ScalarField
template<class SpaceT>
Teuchos::RCP< MultiField<ModeField<ScalarField<SpaceT> > > >
createMultiModeScalarField( const Teuchos::RCP<const SpaceT>& space, int n=1 ) {

  auto scal = create< ModeField<ScalarField<SpaceT> > >( space );

  return( createMultiField<ModeField<ScalarField<SpaceT> > >( *scal, n ) );

}



/// \relates MultiField
/// \relates ModeField
/// \relates VectorField
template<class SpaceT>
Teuchos::RCP< MultiField<ModeField<VectorField<SpaceT> > > > createInitMVF(
    EFlowType flowType,
    const Teuchos::RCP<const SpaceT>& space,
    typename SpaceT::Scalar re=1.,
    typename SpaceT::Scalar omega=1.,
    typename SpaceT::Scalar px=1.) {

  auto vel  = create< ModeField<VectorField<SpaceT> > >( space );

  auto u   = Pimpact::createMultiField< ModeField<VectorField<SpaceT> > >(*vel,1);

  switch( EFlowType( flowType ) ) {
  case Pimpact::Zero2DFlow :
    u->getField(0).getCFieldPtr()->initField( ZeroFlow );
    u->getField(0).getSFieldPtr()->initField( ZeroFlow );
    break;
  case Pimpact::Poiseuille_inX :
    u->getField(0).getCFieldPtr()->initField( PoiseuilleFlow2D_inX );
    u->getField(0).getSFieldPtr()->initField( ZeroFlow );
    break;
  case Poiseuille_inY :
    u->getField(0).getCFieldPtr()->initField( PoiseuilleFlow2D_inY );
    u->getField(0).getSFieldPtr()->initField( ZeroFlow );
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
    u->getField(0).getSFieldPtr()->initField( ZeroFlow, re, omega, px );
    break;
  case Streaming2DFlow2 :
    u->getField(0).getCFieldPtr()->initField( Streaming2DC, re, omega, px );
    u->getField(0).getSFieldPtr()->initField( Streaming2DS, re, omega, px );
    break;
  case Streaming2DFlow3 :
    u->getField(0).getCFieldPtr()->initField( Streaming2DC, re, omega, px );
    u->getField(0).getSFieldPtr()->initField( Streaming2DS, re, omega, px );
    break;
  case OscilatingDisc2D:
    break;
  case OscilatingDisc2DVel:
    break;
  }
  return( u );
}



/// \relates MultiField
/// \relates ModeField
/// \relates ScalarField
template<class SpaceT>
Teuchos::RCP< MultiField<ModeField<ScalarField<SpaceT> > > > createInitMSF(
    const Teuchos::RCP< const SpaceT >& space ) {

  auto sca = Pimpact::create< ModeField<ScalarField<SpaceT> > >( space );

  return( Pimpact::createMultiField<ModeField<ScalarField<SpaceT> > >(*sca,1) );

}



/// \brief creates a multi-harmonic scalar field.
///
/// \relates MultiHarmonicField
/// \relates ScalarField
/// \param space scalar Vector Space to which returned vector belongs
/// \param nf amount of modes are created in multi-harmonic field
/// \return a multi-harmonic \c ScalarField
template<class SpaceT>
Teuchos::RCP< MultiHarmonicField< ScalarField<SpaceT> > > createMultiHarmonicScalarField(
    const Teuchos::RCP<const SpaceT >& space,
    int nf) {
  auto field0 = createScalarField<SpaceT>( space );
  auto mfield = create< ModeField< ScalarField<SpaceT> > >( space );
  auto fields = createMultiField< ModeField< ScalarField<SpaceT> > >( *mfield, nf );
  return(
      Teuchos::rcp(
          new MultiHarmonicField< ScalarField<SpaceT> >( field0, fields ) ) );
}



/// \brief creates a multi-harmonic vector field.
///
/// \relates MultiHarmonicField
/// \relates VectorField
/// \param nf amount of modes
/// \param space
/// \return field vector
template<class SpaceT>
Teuchos::RCP< MultiHarmonicField< VectorField<SpaceT> > >
createMultiHarmonicVectorField(
    const Teuchos::RCP< const SpaceT>& space,
    int nf ) {
  auto field0 = create<Pimpact::VectorField>( space );
  auto mfield = create< ModeField< VectorField<SpaceT> > >( space );
  auto fields = createMultiField< ModeField< VectorField<SpaceT> > >( *mfield, nf );
  return( Teuchos::rcp(
      new MultiHarmonicField< VectorField<SpaceT> >( field0, fields ) ) );
}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_FIELDFACTORY_HPP
