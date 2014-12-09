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
template<class ConvVWrapT>
class MultiHarmonicConvectionOp {

public:

  typedef typename ConvVWrapT::SpaceT SpaceT;

  typedef MultiHarmonicField< VectorField<SpaceT> >  DomainFieldT;
  typedef MultiHarmonicField< VectorField<SpaceT> >  RangeFieldT;

  typedef typename ConvVWrapT::FieldTensor FieldTensor;

protected:

  Teuchos::RCP<const ConvVWrapT> op_;

  Teuchos::RCP< ConvectionField<SpaceT> > wind0_;
  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > windc_;
  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > winds_;

public:

  /// \todo get nf from grid
  MultiHarmonicConvectionOp( const Teuchos::RCP<const ConvVWrapT>& op, int nf ):
    op_( op ),
    wind0_( create<ConvectionField>(op->space()) ),
    windc_( nf ),
    winds_( nf ) {

    for( int i=0; i<nf; ++i ) {
      windc_[i] = create<ConvectionField>( op->space() );
      winds_[i] = create<ConvectionField>( op->space() );
    }

  };

  void assignField( const DomainFieldT& mv ) {

    wind0_->assignField( mv.getConst0Field() );
    int Nf = mv.getNumberModes();

    for( int i=0; i<Nf; ++i ) {
      windc_[i]->assignField( mv.getConstCField(i) );
      winds_[i]->assignField( mv.getConstSField(i) );
    }

  };

//protected:
//public:

//  void apply( const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {
  void apply( const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {

    int Nf = z.getNumberModes();
    if( init_yes )
      z.init( 0. );

    // computing zero mode of y
//    op_->apply( x.getConst0Field(), y.getConst0Field(), z.get0Field(), 1.);
    op_->apply( wind0_->get(), y.getConst0Field(), z.get0Field(), 1.);

    for( int i=1; i<=Nf; ++i ) {
//      op_->apply( x.getConstCField(i-1), y.getConstCField(i-1), z.get0Field(), 0.5 );
//      op_->apply( x.getConstSField(i-1), y.getConstSField(i-1), z.get0Field(), 0.5 );
      op_->apply( windc_[i-1]->get(), y.getConstCField(i-1), z.get0Field(), 0.5 );
      op_->apply( winds_[i-1]->get(), y.getConstSField(i-1), z.get0Field(), 0.5 );
    }


    // computing cos mode of y
    for( int i=1; i<=Nf; ++i ) {
//      op_->apply( x.getConst0Field(), y.getConstCField(i-1), z.getCField(i-1), 1. );
      op_->apply( wind0_->get(), y.getConstCField(i-1), z.getCField(i-1), 1. );

//      op_->apply( x.getConstCField(i-1), y.getConst0Field(), z.getCField(i-1), 1. );
      op_->apply( windc_[i-1]->get(), y.getConst0Field(), z.getCField(i-1), 1. );

      for( int k=1; k+i<=Nf; ++k ) {
//        op_->apply( x.getConstCField(k+i-1), y.getConstCField(k-1), z.getCField(i-1), 0.5 );
        op_->apply( windc_[k+i-1]->get(), y.getConstCField(k-1), z.getCField(i-1), 0.5 );

//        op_->apply( x.getConstCField(k-1), y.getConstCField(k+i-1), z.getCField(i-1), 0.5 );
        op_->apply( windc_[k-1]->get(), y.getConstCField(k+i-1), z.getCField(i-1), 0.5 );

//        op_->apply( x.getConstSField(k+i-1), y.getConstSField(k-1), z.getCField(i-1), 0.5 );
        op_->apply( winds_[k+i-1]->get(), y.getConstSField(k-1), z.getCField(i-1), 0.5 );

//        op_->apply( x.getConstSField(k-1), y.getConstSField(k+i-1), z.getCField(i-1), 0.5 );
        op_->apply( winds_[k-1]->get(), y.getConstSField(k+i-1), z.getCField(i-1), 0.5 );
      }
    }

    // computing sin mode of y
    for( int i=1; i<=Nf; ++i ) {
//      op_->apply( x.getConst0Field(), y.getConstSField(i-1), z.getSField(i-1), 1. );
      op_->apply( wind0_->get(), y.getConstSField(i-1), z.getSField(i-1), 1. );

//      op_->apply( x.getConstSField(i-1), y.getConst0Field(), z.getSField(i-1), 1. );
      op_->apply( winds_[i-1]->get(), y.getConst0Field(), z.getSField(i-1), 1. );

      for( int k=1; k+i<=Nf; ++k ) {
//        op_->apply( x.getConstCField(k+i-1), y.getConstSField(k-1), z.getSField(i-1), -0.5 );
        op_->apply( windc_[k+i-1]->get(), y.getConstSField(k-1), z.getSField(i-1), -0.5 );

//        op_->apply( x.getConstCField(k-1), y.getConstSField(k+i-1), z.getSField(i-1), 0.5 );
        op_->apply( windc_[k-1]->get(), y.getConstSField(k+i-1), z.getSField(i-1), 0.5 );

//        op_->apply( x.getConstSField(k+i-1), y.getConstCField(k-1), z.getSField(i-1), 0.5 );
        op_->apply( winds_[k+i-1]->get(), y.getConstCField(k-1), z.getSField(i-1), 0.5 );

//        op_->apply( x.getConstSField(k-1), y.getConstCField(k+i-1), z.getSField(i-1), -0.5 );
        op_->apply( winds_[k-1]->get(), y.getConstCField(k+i-1), z.getSField(i-1), -0.5 );
      }
    }

    // strange terms
    int i;
    for( int k=1; k<=Nf; ++k ) {
      for( int l=1; l<=Nf; ++l ) {
        i = k+l;
        if( i<=Nf ) {
//          op_->apply( x.getConstCField(k-1), y.getConstCField(l-1), z.getCField(i-1), 0.5 );
//          op_->apply( x.getConstSField(k-1), y.getConstSField(l-1), z.getCField(i-1), -0.5 );
          op_->apply( windc_[k-1]->get(), y.getConstCField(l-1), z.getCField(i-1), 0.5 );
          op_->apply( winds_[k-1]->get(), y.getConstSField(l-1), z.getCField(i-1), -0.5 );

//          op_->apply( x.getConstCField(k-1), y.getConstSField(l-1), z.getSField(i-1), 0.5 );
//          op_->apply( x.getConstSField(k-1), y.getConstCField(l-1), z.getSField(i-1), 0.5 );
          op_->apply( windc_[k-1]->get(), y.getConstSField(l-1), z.getSField(i-1), 0.5 );
          op_->apply( winds_[k-1]->get(), y.getConstCField(l-1), z.getSField(i-1), 0.5 );
        }
      }
    }
  }


//  void apply(const DomainFieldT& x, RangeFieldT& y) const {
//    apply(x,x,y);
//  }


  bool hasApplyTranspose() const { return( false ); }

}; // end of class MultiHarmonicNonlinearOp



/// \relates MultiHarmonicConvectionOp
template<class SpaceT>
Teuchos::RCP<MultiHarmonicConvectionOp<ConvectionVWrap<ConvectionSOp<SpaceT> > > >
createMultiHarmonicConvectionOp( const Teuchos::RCP<const SpaceT>& space, int nf ) {

  auto sop = Pimpact::create<Pimpact::ConvectionSOp>( space ) ;
  auto wrap = Pimpact::create<Pimpact::ConvectionVWrap>( sop );

  return( Teuchos::rcp( new MultiHarmonicConvectionOp<ConvectionVWrap<ConvectionSOp<SpaceT> > >( wrap, nf ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULTICONVECTIONOP_HPP
