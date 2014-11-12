#pragma once
#ifndef PIMPACT_CONVECTIONOP_HPP
#define PIMPACT_CONVECTIONOP_HPP


#include "Pimpact_Types.hpp"

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
/// \deprecated depends on an initialized IMPACT. But can be used to test consistency with new code use \c ConvectionVOp instead.
template<class ST>
class ConvectionOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef VectorField<SpaceT>  DomainFieldT;
  typedef VectorField<SpaceT>  RangeFieldT;

private:

  Teuchos::RCP<const SpaceT> space_;

public:

  ConvectionOp( const Teuchos::RCP<const SpaceT>& space  ):
    space_(space) {};

  void assignField( const DomainFieldT& mv ) {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    apply( x, x, y);

  }

  void apply( const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z, Scalar mul=0. ) const {

    for( int vel_dir=0; vel_dir<space_->dim(); ++vel_dir )
      x.exchange( vel_dir, vel_dir );
    y.exchange();

    if( std::abs(mul) < 1.e-12 ) {
      z.init( 0. );
      mul = 1.;
    }
    OP_nonlinear(
        x.sFields_[0]->getRawPtr(),x.sFields_[1]->getRawPtr(),x.sFields_[2]->getRawPtr(),
        y.sFields_[0]->getRawPtr(),y.sFields_[1]->getRawPtr(),y.sFields_[2]->getRawPtr(),
        z.sFields_[0]->getRawPtr(),z.sFields_[1]->getRawPtr(),z.sFields_[2]->getRawPtr(),
        mul );

    z.changed();

  }


  bool hasApplyTranspose() const { return( false ); }


}; // end of class ConvectionOp



/// \relates ConvectionOp
template<class SpaceT>
Teuchos::RCP<ConvectionOp<SpaceT> > createConvectionOp(
    const Teuchos::RCP<const SpaceT>& space ) {
  return( Teuchos::rcp( new ConvectionOp<SpaceT>( space ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONOP_HPP
