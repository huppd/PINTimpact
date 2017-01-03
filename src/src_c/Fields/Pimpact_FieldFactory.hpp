#pragma once
#ifndef PIMPACT_FIELDFACTORY_HPP
#define PIMPACT_FIELDFACTORY_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_Fields.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"
#include "Pimpact_Space.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {


/// \relates MultiField
/// \relates ModeField
/// \relates VectorField
template<class SpaceT>
Teuchos::RCP< MultiField<ModeField<VectorField<SpaceT> > > >
createMultiModeVectorField( const Teuchos::RCP<const SpaceT>& space, int n=1 ) {

	Teuchos::RCP< ModeField<VectorField<SpaceT> > > vel = create< ModeField< VectorField<SpaceT> > >( space );

	return( createMultiField<ModeField<VectorField<SpaceT> > >( *vel, n ) );
}



/// \relates MultiField
/// \relates ModeField
/// \relates ScalarField
template<class SpaceT>
Teuchos::RCP< MultiField<ModeField<ScalarField<SpaceT> > > >
createMultiModeScalarField( const Teuchos::RCP<const SpaceT>& space, int n=1 ) {

  
	Teuchos::RCP< MultiField<ModeField<ScalarField<SpaceT> > > >
		scal = create< ModeField<ScalarField<SpaceT> > >( space );

  return( createMultiField<ModeField<ScalarField<SpaceT> > >( *scal, n ) );
}




/// \relates MultiField
/// \relates ModeField
/// \relates ScalarField
template<class SpaceT>
Teuchos::RCP< MultiField<ModeField<ScalarField<SpaceT> > > >
createInitMSF(
    const Teuchos::RCP< const SpaceT >& space ) {

  
	Teuchos::RCP< ModeField<ScalarField<SpaceT> > >
		sca = Pimpact::create< ModeField<ScalarField<SpaceT> > >( space );

  return( Pimpact::createMultiField<ModeField<ScalarField<SpaceT> > >(*sca,1) );
}



/// \brief creates a multi-harmonic scalar field.
///
/// \relates MultiHarmonicField
/// \relates ScalarField
/// \param space scalar Vector Space to which returned vector belongs
/// \param nf amount of modes are created in multi-harmonic field
/// \return a multi-harmonic \c ScalarField
/// \deprecated nf
template<class SpaceT>
Teuchos::RCP< MultiHarmonicField< ScalarField<SpaceT> > > createMultiHarmonicScalarField(
    const Teuchos::RCP<const SpaceT >& space,
    int nf) {
  return(
      Teuchos::rcp(
          new MultiHarmonicField< ScalarField<SpaceT> >(space) ) );
}



/// \brief creates a multi-harmonic vector field.
///
/// \relates MultiHarmonicField
/// \relates VectorField
/// \param space
/// \return field vector
template<class SpaceT>
Teuchos::RCP< MultiHarmonicField< VectorField<SpaceT> > >
createMultiHarmonicVectorField(
    const Teuchos::RCP< const SpaceT>& space ) {
  return( Teuchos::rcp(
      new MultiHarmonicField< VectorField<SpaceT> >(space) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_FIELDFACTORY_HPP
