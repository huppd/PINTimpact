#pragma once
#ifndef PIMPACT_MULTICONVECTIONOP_HPP
#define PIMPACT_MULTICONVECTIONOP_HPP


#include "Pimpact_Types.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"

#include "Pimpact_ConvectionVOp.hpp"



namespace Pimpact {



/// \ingroup MultiHarmonicOperator
template<class ST>
class MultiHarmonicConvectionOp {

//  typedef typename SpaceT::Scalar Scalar;
//  typedef typename SpaceT::Ordinal Ordinal;

public:

  typedef ST SpaceT;

  typedef MultiHarmonicField< VectorField<SpaceT> >  DomainFieldT;
  typedef MultiHarmonicField< VectorField<SpaceT> >  RangeFieldT;

protected:

  Teuchos::RCP<const ConvectionVOp<SpaceT> > op_;

public:

  MultiHarmonicConvectionOp( const Teuchos::RCP<const SpaceT>& space ):
    op_( createConvectionVOp(space) ) {};

  void assignField( const DomainFieldT& mv ) {
  };

protected:

  void apply( const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {

    int Nf = x.getNumberModes();
    if( init_yes )
      z.init( 0. );

    // computing zero mode of y
    op_->apply( x.getConst0Field(), y.getConst0Field(), z.get0Field(), 1.);

    for( int i=1; i<=Nf; ++i ) {
      op_->apply( x.getConstCField(i-1), y.getConstCField(i-1), z.get0Field(), 0.5 );
      op_->apply( x.getConstSField(i-1), y.getConstSField(i-1), z.get0Field(), 0.5 );
    }


    // computing cos mode of y
    for( int i=1; i<=Nf; ++i ) {
      op_->apply( x.getConst0Field(), y.getConstCField(i-1), z.getCField(i-1), 1. );

      op_->apply( x.getConstCField(i-1), y.getConst0Field(), z.getCField(i-1), 1. );

      for( int k=1; k+i<=Nf; ++k ) {
        op_->apply( x.getConstCField(k+i-1), y.getConstCField(k-1), z.getCField(i-1), 0.5 );

        op_->apply( x.getConstCField(k-1), y.getConstCField(k+i-1), z.getCField(i-1), 0.5 );

        op_->apply( x.getConstSField(k+i-1), y.getConstSField(k-1), z.getCField(i-1), 0.5 );

        op_->apply( x.getConstSField(k-1), y.getConstSField(k+i-1), z.getCField(i-1), 0.5 );
      }
    }

    // computing sin mode of y
    for( int i=1; i<=Nf; ++i ) {
      op_->apply( x.getConst0Field(), y.getConstSField(i-1), z.getSField(i-1), 1. );

      op_->apply( x.getConstSField(i-1), y.getConst0Field(), z.getSField(i-1), 1. );

      for( int k=1; k+i<=Nf; ++k ) {
        op_->apply( x.getConstCField(k+i-1), y.getConstSField(k-1), z.getSField(i-1), -0.5 );

        op_->apply( x.getConstCField(k-1), y.getConstSField(k+i-1), z.getSField(i-1), 0.5 );

        op_->apply( x.getConstSField(k+i-1), y.getConstCField(k-1), z.getSField(i-1), 0.5 );

        op_->apply( x.getConstSField(k-1), y.getConstCField(k+i-1), z.getSField(i-1), -0.5 );
      }
    }

    // strange terms
    int i;
    for( int k=1; k<=Nf; ++k ) {
      for( int l=1; l<=Nf; ++l ) {
        i = k+l;
        if( i<=Nf ) {
          op_->apply( x.getConstCField(k-1), y.getConstCField(l-1), z.getCField(i-1), 0.5 );
          op_->apply( x.getConstSField(k-1), y.getConstSField(l-1), z.getCField(i-1), -0.5 );

          op_->apply( x.getConstCField(k-1), y.getConstSField(l-1), z.getSField(i-1), 0.5 );
          op_->apply( x.getConstSField(k-1), y.getConstCField(l-1), z.getSField(i-1), 0.5 );
        }
      }
    }
  }

public:

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    apply(x,x,y);
  }


  bool hasApplyTranspose() const { return( false ); }

}; // end of class MultiHarmonicNonlinearOp



/// \relates MultiHarmonicConvectionOp
template<class SpaceT>
Teuchos::RCP<MultiHarmonicConvectionOp<SpaceT> >
createMultiHarmonicConvectionOp( const Teuchos::RCP<const SpaceT>& space ) {

  return( Teuchos::rcp( new MultiHarmonicConvectionOp<SpaceT>( space ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULTICONVECTIONOP_HPP
