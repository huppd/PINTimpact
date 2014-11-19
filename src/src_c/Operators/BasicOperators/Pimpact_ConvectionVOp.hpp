#pragma once
#ifndef PIMPACT_CONVECTIONVOP_HPP
#define PIMPACT_CONVECTIONVOP_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_Space.hpp"

#include "Pimpact_VectorField.hpp"

#include "Pimpact_InterpolateS2VOp.hpp"
#include "Pimpact_ConvectionSOp.hpp"




namespace Pimpact {



/// \brief Convection Operator for Velocity fields
/// \ingroup BaseOperator
/// \relates ConvectionSOp
template<class ST>
class ConvectionVOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  static const int dimension = SpaceT::dimension;

  static const int dimNC = SpaceT::dimNC;

  typedef VectorField<SpaceT>  DomainFieldT;
  typedef VectorField<SpaceT>  RangeFieldT;

private:

  Teuchos::RCP<const SpaceT> space_;

  Teuchos::RCP<const InterpolateS2V<SpaceT> > interpolateS2V_;
  Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension,dimNC> > interpolateV2S_;

  Teuchos::RCP<const ConvectionSOp<SpaceT> > convectionSOp_;

  Teuchos::RCP< ScalarField<SpaceT> > temp_;
  Teuchos::Tuple< Teuchos::Tuple<Teuchos::RCP<ScalarField<SpaceT> >, 3>, 3> u_;

public:

  ConvectionVOp( const Teuchos::RCP<const SpaceT>& space  ):
    space_(space),
    interpolateS2V_( create<InterpolateS2V>(space) ),
    interpolateV2S_( createInterpolateV2S( space ) ),
    convectionSOp_(  create<ConvectionSOp>(space) ),
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

  ConvectionVOp(
      const Teuchos::RCP<const SpaceT>& space,
      const Teuchos::RCP< InterpolateS2V<SpaceT> >& interpolateS2V,
      const Teuchos::RCP< InterpolateV2S<Scalar,Ordinal,dimension,dimNC> >& interpolateV2S,
      const Teuchos::RCP< ConvectionSOp<SpaceT> >& convectionSOp ):
    space_(space),
    interpolateS2V_(interpolateS2V),
    interpolateV2S_(interpolateV2S),
    convectionSOp_(convectionSOp),
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

  void assignField( const DomainFieldT& mv ) {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    apply( x, x, y);

  }

  void apply( const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z, Scalar mul=0. ) const {

    for( int i=0; i<space_->dim(); ++i ) {
      interpolateV2S_->apply( x.getConstField(i), *temp_ );
//      temp_->write( i );
      for( int j=0; j<space_->dim(); ++j ) {
        interpolateS2V_->apply( *temp_, *u_[j][i] );
      }
    }

    if( mul<1.e-12 ) {
      z.init(0);
      mul=1.;
    }

    for( int i=0; i<space_->dim(); ++i ) {
      convectionSOp_->apply( u_[i], y.getConstField(i), z.getField(i), mul );
    }

  }


  bool hasApplyTranspose() const { return( false ); }


}; // end of class ConvectionVOp



/// \relates ConvectionVOp
template<class SpaceT>
Teuchos::RCP<ConvectionVOp<SpaceT> > createConvectionVOp(
    const Teuchos::RCP<const SpaceT>& space ) {
  return( Teuchos::rcp( new ConvectionVOp<SpaceT>( space ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONVOP_HPP
