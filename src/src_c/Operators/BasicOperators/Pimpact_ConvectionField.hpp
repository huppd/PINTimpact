#pragma once
#ifndef PIMPACT_CONVECTIONFIELD_HPP
#define PIMPACT_CONVECTIONFIELD_HPP


#include "Pimpact_VectorField.hpp"

#include "Pimpact_InterpolateS2VOp.hpp"




namespace Pimpact {



/// \brief Stores the wind on differnt grid types.
/// should become template parameter for others such that interplating can be moved from assign to get(different storage needed)
/// \relates ConvectionVOp
template<class ST>
class ConvectionField {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  static const int dimension = SpaceT::dimension;

  static const int dimNC = SpaceT::dimNC;

  typedef VectorField<SpaceT>  DomainFieldT;

protected:

  /// \todo remove
  Teuchos::RCP<const SpaceT> space_;

  Teuchos::RCP<const InterpolateS2V<SpaceT> > interpolateS2V_;
  Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension,dimNC> > interpolateV2S_;

  Teuchos::RCP< ScalarField<SpaceT> > temp_;
  Teuchos::Tuple< Teuchos::Tuple<Teuchos::RCP<ScalarField<SpaceT> >, 3>, 3> u_;

  typedef Teuchos::Tuple< Teuchos::Tuple<Teuchos::RCP<ScalarField<SpaceT> >, 3>, 3> FieldTensor;

public:

  ConvectionField( const Teuchos::RCP<const SpaceT>& space  ):
    space_(space),
    interpolateS2V_( create<InterpolateS2V>(space) ),
    interpolateV2S_( createInterpolateV2S( space ) ),
    temp_( createScalarField<SpaceT>( space ) ),
    u_(
        Teuchos::tuple(
            Teuchos::tuple(
                createScalarField<SpaceT>( space, U ),
                createScalarField<SpaceT>( space, U ),
                createScalarField<SpaceT>( space, U )
            ),
            Teuchos::tuple(
                createScalarField<SpaceT>( space, V ),
                createScalarField<SpaceT>( space, V ),
                createScalarField<SpaceT>( space, V )
            ),
            Teuchos::tuple(
                createScalarField<SpaceT>( space, W ),
                createScalarField<SpaceT>( space, W ),
                createScalarField<SpaceT>( space, W )
            )
        )
    ) {};

  ConvectionField(
      const Teuchos::RCP<const SpaceT>& space,
      const Teuchos::RCP< InterpolateS2V<SpaceT> >& interpolateS2V,
      const Teuchos::RCP< InterpolateV2S<Scalar,Ordinal,dimension,dimNC> >& interpolateV2S,
      const Teuchos::RCP< ConvectionSOp<SpaceT> >& convectionSOp ):
    space_(space),
    interpolateS2V_(interpolateS2V),
    interpolateV2S_(interpolateV2S),
    temp_( createScalarField<SpaceT>( space ) ),
    u_(
        Teuchos::tuple(
            Teuchos::tuple(
                createScalarField<SpaceT>( space, U ),
                createScalarField<SpaceT>( space, U ),
                createScalarField<SpaceT>( space, U )
            ),
            Teuchos::tuple(
                createScalarField<SpaceT>( space, V ),
                createScalarField<SpaceT>( space, V ),
                createScalarField<SpaceT>( space, V )
            ),
            Teuchos::tuple(
                createScalarField<SpaceT>( space, W ),
                createScalarField<SpaceT>( space, W ),
                createScalarField<SpaceT>( space, W )
            )
        )
    ) {};


  void assignField( const DomainFieldT& mv ) const {

    for( int i=0; i<space_->dim(); ++i ) {
      interpolateV2S_->apply( mv.getConstField(i), *temp_ );
      for( int j=0; j<space_->dim(); ++j ) {
        interpolateS2V_->apply( *temp_, *u_[j][i] );
      }
    }
  };


  const FieldTensor& get() const { return( u_ ); }

  FieldTensor& get() { return( u_ ); }


}; // end of class ConvectionField


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONFIELD_HPP
