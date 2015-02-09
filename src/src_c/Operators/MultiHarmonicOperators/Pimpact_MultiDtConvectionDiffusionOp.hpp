#pragma once
#ifndef PIMPACT_MULTIDTCONVECTIONDIFFUSIONOP_HPP
#define PIMPACT_MULTIDTCONVECTIONDIFFUSIONOP_HPP


#include "Pimpact_Types.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"

//#include "Pimpact_ConvectionVOp.hpp"
#include "Pimpact_ConvectionDiffusionSOp.hpp"
#include "Pimpact_ConvectionVWrap.hpp"



namespace Pimpact {



/// \ingroup MultiHarmonicOperator
template<class ST>
class MultiDtConvectionDiffusionOp {

public:

	typedef ST SpaceT;

	typedef typename SpaceT::Scalar Scalar;
	typedef typename SpaceT::Ordinal Ordinal;

  typedef MultiHarmonicField< VectorField<SpaceT> >  DomainFieldT;
  typedef MultiHarmonicField< VectorField<SpaceT> >  RangeFieldT;

  //typedef typename ConvVWrapT::FieldTensor FieldTensor;

protected:

  Teuchos::RCP<const ConvectionVWrap< ConvectionDiffusionSOp<SpaceT> > > op_;

  Teuchos::RCP< ConvectionField<SpaceT> > wind0_;
  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > windc_;
  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > winds_;

public:

  /// \todo get nf from grid
  MultiDtConvectionDiffusionOp( const Teuchos::RCP<const SpaceT>& space, int nf ):
    op_( create< ConvectionVWrap< ConvectionDiffusionSOp<SpaceT> > >(space) ),
    wind0_( create<ConvectionField>(space) ),
    windc_( nf ),
    winds_( nf ) {

    for( int i=0; i<nf; ++i ) {
      windc_[i] = create<ConvectionField>( space );
      winds_[i] = create<ConvectionField>( space );
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


//  void apply( const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {
  void apply( const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {

    int Nf = z.getNumberModes();
		Scalar iRe = 1./op_->space()->getDomain()->getDomainSize()->getRe();
		Scalar a2 = 1./op_->space()->getDomain()->getDomainSize()->getAlpha2();

    // computing zero mode of z
//    op_->apply( x.getConst0Field(), y.getConst0Field(), z.get0Field(), 1.);
    op_->apply( wind0_->get(), y.getConst0Field(), z.get0Field(), 0., 0., 1., iRe );

    for( int i=0; i<Nf; ++i ) {
//      op_->apply( x.getConstCField(i-1), y.getConstCField(i-1), z.get0Field(), 0.5 );
//      op_->apply( x.getConstSField(i-1), y.getConstSField(i-1), z.get0Field(), 0.5 );
      op_->apply( windc_[i]->get(), y.getConstCField(i), z.get0Field(), 1., 0., 0.5, 0. );
      op_->apply( winds_[i]->get(), y.getConstSField(i), z.get0Field(), 1., 0., 0.5, 0. );
    }


    // computing cos mode of z
    for( int i=1; i<=Nf; ++i ) {
      //op_->apply( wind0_->get(), y.getConstCField(i-1), z.getCField(i-1), 1. );
      op_->apply( wind0_->get(), y.getConstCField(i-1), z.getCField(i-1),  0., 0., 1., iRe );

//      op_->apply( x.getConstCField(i-1), y.getConst0Field(), z.getCField(i-1), 1. );
      op_->apply( windc_[i-1]->get(), y.getConst0Field(), z.getCField(i-1), 1., 0., 1., 0. );

      for( int k=1; k+i<=Nf; ++k ) {
//        op_->apply( x.getConstCField(k+i-1), y.getConstCField(k-1), z.getCField(i-1), 0.5 );
        op_->apply( windc_[k+i-1]->get(), y.getConstCField(k-1), z.getCField(i-1), 1., 0., 0.5, 0. );

//        op_->apply( x.getConstCField(k-1), y.getConstCField(k+i-1), z.getCField(i-1), 0.5 );
        op_->apply( windc_[k-1]->get(), y.getConstCField(k+i-1), z.getCField(i-1), 1., 0., 0.5, 0. );

//        op_->apply( x.getConstSField(k+i-1), y.getConstSField(k-1), z.getCField(i-1), 0.5 );
        op_->apply( winds_[k+i-1]->get(), y.getConstSField(k-1), z.getCField(i-1), 1., 0., 0.5, 0. );

//        op_->apply( x.getConstSField(k-1), y.getConstSField(k+i-1), z.getCField(i-1), 0.5 );
        op_->apply( winds_[k-1]->get(), y.getConstSField(k+i-1), z.getCField(i-1), 1., 0., 0.5, 0. );
      }
    }

    // computing sin mode of y
    for( int i=1; i<=Nf; ++i ) {
//      op_->apply( x.getConst0Field(), y.getConstSField(i-1), z.getSField(i-1), 1. );
      op_->apply( wind0_->get(), y.getConstSField(i-1), z.getSField(i-1), 0., 0., 1., iRe );

//      op_->apply( x.getConstSField(i-1), y.getConst0Field(), z.getSField(i-1), 1. );
      op_->apply( winds_[i-1]->get(), y.getConst0Field(), z.getSField(i-1), 1., 0., 1., 0. );

      for( int k=1; k+i<=Nf; ++k ) {
//        op_->apply( x.getConstCField(k+i-1), y.getConstSField(k-1), z.getSField(i-1), -0.5 );
        op_->apply( windc_[k+i-1]->get(), y.getConstSField(k-1), z.getSField(i-1), 1., 0., -0.5, 0. );

//        op_->apply( x.getConstCField(k-1), y.getConstSField(k+i-1), z.getSField(i-1), 0.5 );
        op_->apply( windc_[k-1]->get(), y.getConstSField(k+i-1), z.getSField(i-1), 1., 0.,  0.5, 0. );

//        op_->apply( x.getConstSField(k+i-1), y.getConstCField(k-1), z.getSField(i-1), 0.5 );
        op_->apply( winds_[k+i-1]->get(), y.getConstCField(k-1), z.getSField(i-1), 1., 0.,  0.5, 0. );

//        op_->apply( x.getConstSField(k-1), y.getConstCField(k+i-1), z.getSField(i-1), -0.5 );
        op_->apply( winds_[k-1]->get(), y.getConstCField(k+i-1), z.getSField(i-1), 1., 0., -0.5, 0. );
      }
    }

    // strange terms
    int i;
		Scalar mulI;
    for( int k=1; k<=Nf; ++k ) {
      for( int l=1; l<=Nf; ++l ) {
        i = k+l;
				mulI = (l-1==i-1)?(a2*l):0;
        if( i<=Nf ) {
//          op_->apply( x.getConstCField(k-1), y.getConstCField(l-1), z.getCField(i-1), 0.5 );
//          op_->apply( x.getConstSField(k-1), y.getConstSField(l-1), z.getCField(i-1), -0.5 );
          op_->apply( windc_[k-1]->get(), y.getConstCField(l-1), z.getCField(i-1), 1., 0.,  0.5, 0. );
          op_->apply( winds_[k-1]->get(), y.getConstSField(l-1), z.getCField(i-1), 1., mulI, -0.5, 0. );

//          op_->apply( x.getConstCField(k-1), y.getConstSField(l-1), z.getSField(i-1), 0.5 );
//          op_->apply( x.getConstSField(k-1), y.getConstCField(l-1), z.getSField(i-1), 0.5 );
          op_->apply( windc_[k-1]->get(), y.getConstSField(l-1), z.getSField(i-1), 1., 0., 0.5, 0. );
          op_->apply( winds_[k-1]->get(), y.getConstCField(l-1), z.getSField(i-1), 1., -mulI, 0.5, 0. );
        }
				else {
					if( l-1==i-1 ) {
						z.getCField(i-1)->add( 1., z.getCField(i-1),  mulI, y.getConstSField(i-1) );
						z.getSField(i-1)->add( 1., z.getSField(i-1), -mulI, y.getConstCField(i-1) );
					}

				}
      }
    }
  }



//  void apply(const DomainFieldT& x, RangeFieldT& y) const {
//    apply(x,x,y);
//  }


  bool hasApplyTranspose() const { return( false ); }

}; // end of class MultiDtConvectionDiffusionOp



/// \relates MultiDtConvectionDiffusionOp
template<class SpaceT>
Teuchos::RCP<MultiDtConvectionDiffusionOp<SpaceT> >
createMultiDtConvectionDiffusionOp( const Teuchos::RCP<const SpaceT>& space, int nf ) {

  return( Teuchos::rcp( new MultiDtConvectionDiffusionOp<SpaceT>( space, nf ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULTIDTCONVECTIONDIFFUSIONOP_HPP