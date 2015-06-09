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
    const double* const phi1U, const double* const phi1V, const double* const phi1W,
    const double* const phi2U, const double* const phi2V, const double* const phi2W,
    double* const nl1,   double* const nl2,   double* const nl3,
    const double& mul );

}


/// \ingroup BaseOperator
/// \deprecated depends on an initialized IMPACT. Stays only to test consistency with new code use \c ConvectionVOp instead.
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
        x.getConstRawPtr(0),x.getConstRawPtr(1),x.getConstRawPtr(2),
        y.getConstRawPtr(0),y.getConstRawPtr(1),y.getConstRawPtr(2),
        z.getRawPtr(0),z.getRawPtr(1),z.getRawPtr(2),
        mul );

    z.changed();

  }


  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(space_); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}


}; // end of class ConvectionOp



/// \relates ConvectionOp
template<class SpaceT>
Teuchos::RCP<ConvectionOp<SpaceT> > createConvectionOp(
    const Teuchos::RCP<const SpaceT>& space ) {
  return( Teuchos::rcp( new ConvectionOp<SpaceT>( space ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONOP_HPP
