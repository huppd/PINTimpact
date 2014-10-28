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
template<class Scalar,class Ordinal, int dimension=3>
class ConvectionJacobianOp {

public:

  typedef VectorField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal,dimension>  RangeFieldT;

protected:

  Teuchos::RCP<DomainFieldT> u_;

  Teuchos::RCP<const ConvectionVOp<Scalar,Ordinal,dimension> > convectionVOp_;

  const bool isNewton_;

public:

  ConvectionJacobianOp(
      const Teuchos::RCP<const Space<Scalar,Ordinal,dimension> >& space,
      const bool& isNewton=true ):
    u_( Teuchos::null ),
    convectionVOp_( createConvectionVOp( space ) ),
    isNewton_(isNewton) {};

  ConvectionJacobianOp(
      const Teuchos::RCP<const ConvectionVOp<Scalar,Ordinal,dimension> >& convectionVOp,
      const Teuchos::RCP<DomainFieldT>& u,
      const bool& isNewton=true ):
    u_( Teuchos::null ),
    convectionVOp_(convectionVOp),
    isNewton_(isNewton) {};


  /// \todo danger different behavior BC are not necessarily treated, assign
  ///doesnot care about BC( fixed, but always whole copy is not efficient)
  void assignField( const DomainFieldT& mv ) {

    if( u_.is_null() )
      u_ = mv.clone(DeepCopy);
    else
      u_->assign( mv );

    u_->exchange();

  };

  void apply( const DomainFieldT& x, RangeFieldT& y, Scalar mul=0. ) const {

    if( mul<1.e-12 ) {
      y.init(0.);
      mul=1.;
    }

    convectionVOp_->apply( *u_, x, y, mul );

    if( isNewton_ )
      convectionVOp_->apply(  x, *u_, y, mul );

  }

  bool hasApplyTranspose() const { return( false ); }


}; // end of class ConvectionJacobianOp



/// \relates ConvectionJacobianOp
template< class S=double, class O=int, int d=3>
Teuchos::RCP<ConvectionJacobianOp<S,O,d> > createConvectionJacobianOp(
    const Teuchos::RCP<const Space<S,O,d> >& space,
    const bool& isNewton=true ) {

    return( Teuchos::rcp( new ConvectionJacobianOp<S,O,d>( space, isNewton ) ) );

}

/// \relates ConvectionJacobianOp
template< class S=double, class O=int, int d=3>
Teuchos::RCP<ConvectionJacobianOp<S,O,d> > createConvectionJacobianOp(
    const Teuchos::RCP<const ConvectionVOp<S,O,d> >& convectionVOp,
    const bool& isNewton=true ) {

    return( Teuchos::rcp( new ConvectionJacobianOp<S,O,d>( convectionVOp, isNewton ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_NONLINEARJACOBIANOP_HPP
