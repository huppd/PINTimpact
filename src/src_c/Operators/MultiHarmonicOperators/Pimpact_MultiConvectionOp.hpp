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


  void apply( const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {

    int Nf = z.getNumberModes();

    // computing zero mode of y
    op_->apply( wind0_->get(), y.getConst0Field(), z.get0Field(), 0., 0., 1., 0.);

    for( int i=1; i<=Nf; ++i ) {
      op_->apply( windc_[i-1]->get(), y.getConstCField(i-1), z.get0Field(), 1., 0., 0.5, 0. );
      op_->apply( winds_[i-1]->get(), y.getConstSField(i-1), z.get0Field(), 1., 0., 0.5, 0. );
    }


    // computing cos mode of y
    for( int i=1; i<=Nf; ++i ) {
      op_->apply( wind0_->get(), y.getConstCField(i-1), z.getCField(i-1), 0., 0., 1., 0. );
      op_->apply( windc_[i-1]->get(), y.getConst0Field(), z.getCField(i-1), 1., 0., 1., 0. );

      for( int k=1; k+i<=Nf; ++k ) {
        op_->apply( windc_[k+i-1]->get(), y.getConstCField(k-1), z.getCField(i-1), 1., 0., 0.5, 0. );

        op_->apply( windc_[k-1]->get(), y.getConstCField(k+i-1), z.getCField(i-1), 1., 0., 0.5, 0. );

        op_->apply( winds_[k+i-1]->get(), y.getConstSField(k-1), z.getCField(i-1), 1., 0., 0.5, 0. );

        op_->apply( winds_[k-1]->get(), y.getConstSField(k+i-1), z.getCField(i-1), 1., 0., 0.5, 0. );
      }
    }

    // computing sin mode of y
    for( int i=1; i<=Nf; ++i ) {
      op_->apply( wind0_->get(),   y.getConstSField(i-1), z.getSField(i-1), 0., 0., 1., 0. );
      op_->apply( winds_[i-1]->get(), y.getConst0Field(), z.getSField(i-1), 1., 0., 1., 0. );

      for( int k=1; k+i<=Nf; ++k ) {
        op_->apply( windc_[k+i-1]->get(), y.getConstSField(k-1), z.getSField(i-1), 1., 0., -0.5, 0. );
        op_->apply( windc_[k-1]->get(), y.getConstSField(k+i-1), z.getSField(i-1), 1., 0.,  0.5, 0. );
        op_->apply( winds_[k+i-1]->get(), y.getConstCField(k-1), z.getSField(i-1), 1., 0., 0.5, 0. );
        op_->apply( winds_[k-1]->get(), y.getConstCField(k+i-1), z.getSField(i-1), 1., 0., -0.5, 0. );
      }
    }

    // strange terms
    int i;
    for( int k=1; k<=Nf; ++k ) {
      for( int l=1; l<=Nf; ++l ) {
        i = k+l;
        if( i<=Nf ) {
          op_->apply( windc_[k-1]->get(), y.getConstCField(l-1), z.getCField(i-1), 1., 0., 0.5, 0. );
          op_->apply( winds_[k-1]->get(), y.getConstSField(l-1), z.getCField(i-1), 1., 0.,-0.5, 0. );

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

	const std::string getLabel() const { return( "MultiHarmonicConvectionOp " ); };

  void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
		op_->print( out );
  }

}; // end of class MultiHarmonicConvectionOp



/// \relates MultiHarmonicConvectionOp
template<class SpaceT>
Teuchos::RCP<MultiHarmonicConvectionOp<ConvectionVWrap<ConvectionSOp<SpaceT> > > >
createMultiHarmonicConvectionOp( const Teuchos::RCP<const SpaceT>& space, int nf=-1 ) {

  auto sop = Pimpact::create<Pimpact::ConvectionSOp>( space ) ;
  auto wrap = Pimpact::create<Pimpact::ConvectionVWrap>( sop );

	if( nf<0 )
		return(
				Teuchos::rcp(
					new MultiHarmonicConvectionOp<ConvectionVWrap<ConvectionSOp<SpaceT> > >( wrap, space->nGlo(3) ) ) );

  return(
			Teuchos::rcp(
				new MultiHarmonicConvectionOp<ConvectionVWrap<ConvectionSOp<SpaceT> > >( wrap, nf ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULTICONVECTIONOP_HPP
