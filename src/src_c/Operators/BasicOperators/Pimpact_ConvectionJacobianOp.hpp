#pragma once
#ifndef PIMPACT_NONLINEARJACOBIANOP_HPP
#define PIMPACT_NONLINEARJACOBIANOP_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"

#include "Pimpact_ConvectionVOp.hpp"




namespace Pimpact {



/// \ingroup BaseOperator
/// \note u_ has to contain appropriate BC, temp_ and y doesnot matter, x should
///   have zero boundary conditions if used for linear solver.
/// \note if heavily used one should do the interpolation in the assignField method.
/// \relates ConvectionVOp
/// \todo remodel
/// \deprecated implement Newton in stencils
template<class ST>
class ConvectionJacobianOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef VectorField<SpaceT>  DomainFieldT;
  typedef VectorField<SpaceT>  RangeFieldT;

protected:

  Teuchos::RCP<DomainFieldT> u_;

  Teuchos::RCP<const ConvectionVOp<SpaceT> > convectionVOp_;
  Teuchos::RCP<const ConvectionVOp<SpaceT> > convectionVOp2_;

  const bool isNewton_;

public:

  ConvectionJacobianOp(
      const Teuchos::RCP<const SpaceT>& space,
      const bool& isNewton=true ):
        u_( Teuchos::null ),
        convectionVOp_( createConvectionVOp( space ) ),
        isNewton_(isNewton) {
    if( isNewton_ )
      convectionVOp2_ = createConvectionVOp( space );
  };


  /// \todo danger different behavior BC are not necessarily treated, assign
  ///doesnot care about BC( fixed, but always whole copy is not efficient)
  void assignField( const DomainFieldT& mv ) {

    convectionVOp_->assignField( mv );

    if( isNewton_ ) {
      if( u_.is_null() )
        u_ = mv.clone(DeepCopy);
      else
        u_->assign( mv );

      u_->exchange();
    }
  };

  void apply( const DomainFieldT& x, RangeFieldT& y, Scalar mul=0. ) const {

    if( mul<1.e-12 ) {
      y.init(0.);
      mul=1.;
    }

    convectionVOp_->apply( x, y, mul );

    if( isNewton_ )
      convectionVOp2_->apply(  x, *u_, y, mul );

  }


	Teuchos::RCP<const SpaceT> space() const { return(convectionVOp_->space()); };

  bool hasApplyTranspose() const { return( false ); }


}; // end of class ConvectionJacobianOp



/// \relates ConvectionJacobianOp
template< class SpaceT>
Teuchos::RCP<ConvectionJacobianOp<SpaceT> > createConvectionJacobianOp(
    const Teuchos::RCP<const SpaceT>& space,
    const bool& isNewton=true ) {

  return( Teuchos::rcp( new ConvectionJacobianOp<SpaceT>( space, isNewton ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_NONLINEARJACOBIANOP_HPP
