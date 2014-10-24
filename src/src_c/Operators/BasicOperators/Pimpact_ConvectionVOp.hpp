#pragma once
#ifndef PIMPACT_CONVECTIONVOP_HPP
#define PIMPACT_CONVECTIONVOP_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_Space.hpp"

#include "Pimpact_VectorField.hpp"

#include "Pimpact_InterpolateS2VOp.hpp"
#include "Pimpact_ConvectionSOp.hpp"



namespace Pimpact {


extern "C" {

void OP_nonlinear(
    double* const phi1U, double* const phi1V, double* const phi1W,
    double* const phi2U, double* const phi2V, double* const phi2W,
    double* const nl1,   double* const nl2,   double* const nl3,
    const double& mul );

}


/// \ingroup BaseOperator
/// this Operator should be consisten to the old one... but will be probably redundant soon
template<class Scalar,class Ordinal, int dimension=3>
class ConvectionVOp {

public:

  typedef VectorField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal,dimension>  RangeFieldT;

private:

  Teuchos::RCP<const Space<Scalar,Ordinal,dimension> > space_;

  Teuchos::RCP<const InterpolateS2V<Scalar,Ordinal,dimension> > interpolateS2V_;
  Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension> > interpolateV2S_;

  Teuchos::RCP<const ConvectionSOp<Scalar,Ordinal,dimension> > convectionSOp_;

  Teuchos::RCP< ScalarField<Scalar,Ordinal,dimension> > temp_;
  Teuchos::Tuple< Teuchos::Tuple<Teuchos::RCP<ScalarField<Scalar,Ordinal,dimension> >, 3>, 3> u_;

public:

  ConvectionVOp( const Teuchos::RCP<const Space<Scalar,Ordinal,dimension> >& space  ):
    space_(space),
    interpolateS2V_( createInterpolateS2V<Scalar,Ordinal,dimension>(space) ),
    interpolateV2S_( space->getInterpolateV2S() ),
    convectionSOp_( createConvectionSOp<Scalar,Ordinal,dimension>( space ) ),
    temp_( createScalarField<Scalar,Ordinal,dimension>( space ) ),
    u_(
        Teuchos::tuple(
            Teuchos::tuple(
                createScalarField<Scalar,Ordinal,dimension>( space, U ),
                createScalarField<Scalar,Ordinal,dimension>( space, U ),
                createScalarField<Scalar,Ordinal,dimension>( space, U )
            ),
            Teuchos::tuple(
                createScalarField<Scalar,Ordinal,dimension>( space, V ),
                createScalarField<Scalar,Ordinal,dimension>( space, V ),
                createScalarField<Scalar,Ordinal,dimension>( space, V )
            ),
            Teuchos::tuple(
                createScalarField<Scalar,Ordinal,dimension>( space, W ),
                createScalarField<Scalar,Ordinal,dimension>( space, W ),
                createScalarField<Scalar,Ordinal,dimension>( space, W )
            )
        )
    ) {};

  ConvectionVOp(
      const Teuchos::RCP<const Space<Scalar,Ordinal,dimension> >& space,
      const Teuchos::RCP< InterpolateS2V<Scalar,Ordinal,dimension> >& interpolateS2V,
      const Teuchos::RCP< InterpolateV2S<Scalar,Ordinal,dimension> >& interpolateV2S,
      const Teuchos::RCP< ConvectionSOp<Scalar,Ordinal,dimension> >& convectionSOp ):
    space_(space),
    interpolateS2V_(interpolateS2V),
    interpolateV2S_(interpolateV2S),
    convectionSOp_(convectionSOp),
    temp_( createScalarField<Scalar,Ordinal,dimension>( space ) ),
    u_(
        Teuchos::tuple(
            Teuchos::tuple(
                createScalarField<Scalar,Ordinal,dimension>( space, U ),
                createScalarField<Scalar,Ordinal,dimension>( space, U ),
                createScalarField<Scalar,Ordinal,dimension>( space, U )
            ),
            Teuchos::tuple(
                createScalarField<Scalar,Ordinal,dimension>( space, V ),
                createScalarField<Scalar,Ordinal,dimension>( space, V ),
                createScalarField<Scalar,Ordinal,dimension>( space, V )
            ),
            Teuchos::tuple(
                createScalarField<Scalar,Ordinal,dimension>( space, W ),
                createScalarField<Scalar,Ordinal,dimension>( space, W ),
                createScalarField<Scalar,Ordinal,dimension>( space, W )
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
      for( int j=0; j<space_->dim(); ++j ) {
        interpolateS2V_->apply( *temp_, *u_[j][i] );
        u_[j][i]->write((1+i)*(1+j));
      }
    }

    for( int i=0; i<space_->dim(); ++i ) {
      convectionSOp_->apply( u_[i], y.getConstField(i), z.getField(i) );

    }

  }


  bool hasApplyTranspose() const { return( false ); }


}; // end of class ConvectionVOp



/// \relates ConvectionVOp
template< class S=double, class O=int, int d=3 >
Teuchos::RCP<ConvectionVOp<S,O,d> > createConvectionVOp(
    const Teuchos::RCP<const Space<S,O,d> >& space ) {
  return( Teuchos::rcp( new ConvectionVOp<S,O,d>( space ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONVOP_HPP
