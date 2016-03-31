#pragma once
#ifndef PIMPACT_MULTIDIAGCONVECTIONJACOBIANOP_HPP
#define PIMPACT_MULTIDIAGCONVECTIONJACOBIANOP_HPP


#include "Pimpact_FieldFactory.hpp"
#include "Pimpact_MultiHarmonicField.hpp"
#include "Pimpact_NonlinearOp.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {



/// \ingroup MultiHarmonicOperator
template<class ST>
class MultiDiagConvectionJacobianOp {

  using Scalar = typename ST::Scalar;
  using Ordinal = typename ST::Ordinal;

public:

  using SpaceT = ST;

  using DomainFieldT = MultiHarmonicField< VectorField<SpaceT> >;
  using RangeFieldT = MultiHarmonicField< VectorField<SpaceT> >;

protected:

  Teuchos::RCP<DomainFieldT> u_;
  Teuchos::RCP<NonlinearOp<SpaceT> > op_;

  const bool isNewton_;

public:

  MultiDiagConvectionJacobianOp(
      const Teuchos::RCP<const SpaceT>& space,
      const bool& isNewton=true ):
        u_(Teuchos::null),
        op_( createNonlinearOp<SpaceT>( space ) ),
        isNewton_(isNewton) {};

  void assignField( const DomainFieldT& mv ) {
    if( Teuchos::is_null( u_ ) )
      u_ = mv.clone();
    else
      u_->assign( mv );
  };

protected:

  void apply(const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {
    int Nf = space()->nGlo(3);
    if( init_yes )
      z.init( 0. );

    // computing zero mode of y
    op_->apply( x.getConst0Field(), y.getConst0Field(), z.get0Field(), 1. );
    //    z.get0Field().add( 1., z.getConst0Field(), 1., *temp_ );


    for( int i=1; i<=Nf; ++i ) {
      // computing cos mode of y
      op_->apply( x.getConst0Field(), y.getConstCField(i-1), z.getCField(i-1), 1. );
      //      z.getCField(i-1).add( 1., z.getConstCField(i-1), 1., *temp_ );

      // computing sin mode of y
      op_->apply( x.getConst0Field(), y.getConstSField(i-1), z.getSField(i-1), 1. );
      //      z.getSField(i-1).add( 1., z.getConstSField(i-1), 1., *temp_ );

      if( 2*i<=Nf ) {

        // computing cos mode of y
        op_->apply( x.getConstCField(i+i-1), y.getConstCField(i-1), z.getCField(i-1), 0.5 );
        //        z.getCField(i-1).add( 1., z.getConstCField(i-1), 0.5, *temp_ );

        op_->apply( x.getConstSField(i+i-1), y.getConstSField(i-1), z.getCField(i-1), 0.5 );
        //        z.getCField(i-1).add( 1., z.getConstCField(i-1), 0.5, *temp_ );

        // computing sin mode of y
        op_->apply( x.getConstCField(i+i-1), y.getConstSField(i-1), z.getSField(i-1), -0.5 );
        //        z.getSField(i-1).add( 1., z.getConstSField(i-1), -0.5, *temp_ );

        op_->apply( x.getConstSField(i+i-1), y.getConstCField(i-1), z.getSField(i-1), 0.5 );
        //        z.getSField(i-1).add( 1., z.getConstSField(i-1), 0.5, *temp_ );

      }

    }
  }

public:

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    apply( *u_,  x,  y, true );

    if( isNewton_ )
      apply(  x,  *u_, y, false );

  }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "MultiDiagConvectionJacobianOp "); };

}; // end of class MultiDiagConvectionJacobianOp



/// \relates MultiDiagConvectionJacobianOp
template<class SpaceT>
Teuchos::RCP<MultiDiagConvectionJacobianOp<SpaceT> >
createMultiDiagConvectionJacobianOp(
    const Teuchos::RCP<const SpaceT>& space,
    const bool& isNewton=true ) {

  return( Teuchos::rcp( new MultiDiagConvectionJacobianOp<SpaceT>( space, isNewton ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULTIDIAGCONVECTIONJACOBIANOP_HPP
