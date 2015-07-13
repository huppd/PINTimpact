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

  Teuchos::RCP<ConvectionVWrap< ConvectionDiffusionSOp<SpaceT> > > op_;

  Teuchos::RCP< ConvectionField<SpaceT> > wind0_;
  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > windc_;
  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > winds_;

public:

  /// \todo get nf from grid
  MultiDtConvectionDiffusionOp( const Teuchos::RCP<const SpaceT>& space ):
    op_( create<ConvectionVWrap>( create<ConvectionDiffusionSOp<SpaceT> >(space) ) ),
    wind0_( create<ConvectionField>(space) ),
    windc_( space->nGlo(3) ),
    winds_( space->nGlo(3) ) {

    for( int i=0; i<space->nGlo(3); ++i ) {
      windc_[i] = create<ConvectionField>( space );
      winds_[i] = create<ConvectionField>( space );
    }

  };

  void assignField( const DomainFieldT& mv ) {

//		mv.write(99);

    wind0_->assignField( mv.getConst0Field() );
    int Nf = mv.getNumberModes();

    for( int i=0; i<Nf; ++i ) {
      windc_[i]->assignField( mv.getConstCField(i) );
      winds_[i]->assignField( mv.getConstSField(i) );
    }

  };


  void apply( const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {
		
		int Nf = z.getNumberModes();
		Scalar iRe = 1./op_->space()->getDomain()->getDomainSize()->getRe();
		Scalar a2 = op_->space()->getDomain()->getDomainSize()->getAlpha2()*iRe;

		Scalar mulI;


    // computing zero mode of z
    op_->apply( wind0_->get(), y.getConst0Field(), z.get0Field(), 0., 0., 1., iRe );

    for( int i=0; i<Nf; ++i ) {
      op_->apply( windc_[i]->get(), y.getConstCField(i), z.get0Field(), 1., 0., 0.5, 0. );
      op_->apply( winds_[i]->get(), y.getConstSField(i), z.get0Field(), 1., 0., 0.5, 0. );
    }


    // computing cos mode of z
    for( int i=1; i<=Nf; ++i ) {
      op_->apply( wind0_->get(), y.getConstCField(i-1), z.getCField(i-1),  0., 0., 1., iRe );

      op_->apply( windc_[i-1]->get(), y.getConst0Field(), z.getCField(i-1), 1., 0., 1., 0. );

      for( int k=1; k+i<=Nf; ++k ) {
				mulI = (k-1==i-1)?(a2*i):0;
        op_->apply( windc_[k+i-1]->get(), y.getConstCField(k-1), z.getCField(i-1), 1., 0., 0.5, 0. );

        op_->apply( windc_[k-1]->get(), y.getConstCField(k+i-1), z.getCField(i-1), 1., 0., 0.5, 0. );

        op_->apply( winds_[k+i-1]->get(), y.getConstSField(k-1), z.getCField(i-1), 1., mulI, 0.5, 0. );

        op_->apply( winds_[k-1]->get(), y.getConstSField(k+i-1), z.getCField(i-1), 1., 0., 0.5, 0. );
      }
    }

    // computing sin mode of y
    for( int i=1; i<=Nf; ++i ) {
      op_->apply( wind0_->get(), y.getConstSField(i-1), z.getSField(i-1), 0., 0., 1., iRe );

      op_->apply( winds_[i-1]->get(), y.getConst0Field(), z.getSField(i-1), 1., 0., 1., 0. );

      for( int k=1; k+i<=Nf; ++k ) {
				mulI = (k-1==i-1)?(a2*i):0;
        op_->apply( windc_[k+i-1]->get(), y.getConstSField(k-1), z.getSField(i-1), 1., 0., -0.5, 0. );

        op_->apply( windc_[k-1]->get(), y.getConstSField(k+i-1), z.getSField(i-1), 1., 0.,  0.5, 0. );

        op_->apply( winds_[k+i-1]->get(), y.getConstCField(k-1), z.getSField(i-1), 1., -mulI,  0.5, 0. );

        op_->apply( winds_[k-1]->get(), y.getConstCField(k+i-1), z.getSField(i-1), 1., 0., -0.5, 0. );
      }
    }

		// rest of time
		for( int i=Nf/2+1; i<=Nf; ++i ) {
				mulI = a2*i;
				z.getCFieldPtr(i-1)->add( 1., z.getCField(i-1),  mulI, y.getConstSField(i-1) );
				z.getSFieldPtr(i-1)->add( 1., z.getSField(i-1), -mulI, y.getConstCField(i-1) );
		}


    // strange terms
    int i;
    for( int k=1; k<=Nf; ++k ) {
      for( int l=1; l<=Nf; ++l ) {
        i = k+l;
        if( i<=Nf ) {
          op_->apply( windc_[k-1]->get(), y.getConstCField(l-1), z.getCField(i-1), 1., 0.,  0.5, 0. );
          op_->apply( winds_[k-1]->get(), y.getConstSField(l-1), z.getCField(i-1), 1., 0., -0.5, 0. );

          op_->apply( windc_[k-1]->get(), y.getConstSField(l-1), z.getSField(i-1), 1., 0., 0.5, 0. );
          op_->apply( winds_[k-1]->get(), y.getConstCField(l-1), z.getSField(i-1), 1., 0., 0.5, 0. );
        }
      }
    }
  }



//  void apply(const DomainFieldT& x, RangeFieldT& y) const {
//    apply(x,x,y);
//  }


	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "MHDtConvectionDiffusion" ); };

  void print( std::ostream& out=std::cout ) const {
		out <<  getLabel() << ":\n";
		op_->print( out );
  }

}; // end of class MultiDtConvectionDiffusionOp



/// \relates MultiDtConvectionDiffusionOp
template<class SpaceT>
Teuchos::RCP<MultiDtConvectionDiffusionOp<SpaceT> >
createMultiDtConvectionDiffusionOp( const Teuchos::RCP<const SpaceT>& space ) {

  return( Teuchos::rcp( new MultiDtConvectionDiffusionOp<SpaceT>( space ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULTIDTCONVECTIONDIFFUSIONOP_HPP
