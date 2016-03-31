#pragma once
#ifndef PIMPACT_MULTIDTCONVECTIONDIFFUSIONOP_HPP
#define PIMPACT_MULTIDTCONVECTIONDIFFUSIONOP_HPP


#include "Pimpact_ConvectionDiffusionSOp.hpp"
#include "Pimpact_FieldFactory.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"
#include "Pimpact_NonlinearVWrap.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact {



/// \ingroup MultiHarmonicOperator
template<class ST>
class MultiDtConvectionDiffusionOp {

public:

  using SpaceT = ST;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using DomainFieldT = MultiHarmonicField< VectorField<SpaceT> >;
  using RangeFieldT = MultiHarmonicField< VectorField<SpaceT> >;

protected:

  Teuchos::RCP<NonlinearWrap< ConvectionDiffusionSOp<SpaceT> > > op_;

  Teuchos::RCP< ConvectionField<SpaceT> > wind0_;
  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > windc_;
  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > winds_;

public:

  /// \todo get nf from grid
  MultiDtConvectionDiffusionOp( const Teuchos::RCP<const SpaceT>& space ):
    op_( create<NonlinearWrap>( create<ConvectionDiffusionSOp<SpaceT> >(space) ) ),
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

		mv.exchange();

    wind0_->assignField( mv.getConst0Field() );

    int Nf = space()->nGlo(3);
    for( int i=0; i<Nf; ++i ) {
      windc_[i]->assignField( mv.getConstCField(i) );
      winds_[i]->assignField( mv.getConstSField(i) );
    }

  };


	/// \todo change loops such that only local stuff is computed
  void apply( const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {
		
		int Nf = space()->nGlo(3);
		Scalar iRe = 1./op_->space()->getDomainSize()->getRe();
		Scalar a2 = op_->space()->getDomainSize()->getAlpha2()*iRe;

		Scalar mulI;


		y.exchange();

    // computing zero mode of z
		if( space()->sInd(U,3)<0 ) {

			op_->apply( wind0_->get(), y.getConst0Field(), z.get0Field(), 0., 0., 1., iRe );

			for( int i=0; i<Nf; ++i ) {
				op_->apply( windc_[i]->get(), y.getConstCField(i), z.get0Field(), 1., 0., 0.5, 0. );
				op_->apply( winds_[i]->get(), y.getConstSField(i), z.get0Field(), 1., 0., 0.5, 0. );
			}

		}


    // computing cos mode of z
		for( Ordinal i=std::max(space()->sInd(U,3),0)+1; i<=space()->eInd(U,3); ++i ) {
//    for( int i=1; i<=Nf; ++i ) { // change to seInd
//
      op_->apply( wind0_->get(),      y.getConstCField(i-1), z.getCField(i-1), 0., 0., 1., iRe );

      op_->apply( windc_[i-1]->get(), y.getConst0Field(),    z.getCField(i-1), 1., 0., 1., 0.  );

      for( int k=1; k+i<=Nf; ++k ) { // thats fine
				mulI = (k-1==i-1)?(a2*i):0;
        op_->apply( windc_[k+i-1]->get(), y.getConstCField(k-1),   z.getCField(i-1), 1.,   0., 0.5, 0. );

        op_->apply( windc_[k-1]->get(),   y.getConstCField(k+i-1), z.getCField(i-1), 1.,   0., 0.5, 0. );

        op_->apply( winds_[k+i-1]->get(), y.getConstSField(k-1),   z.getCField(i-1), 1., mulI, 0.5, 0. );

        op_->apply( winds_[k-1]->get(),   y.getConstSField(k+i-1), z.getCField(i-1), 1.,   0., 0.5, 0. );
      }
    }

    // computing sin mode of y
//    for( int i=1; i<=Nf; ++i ) { // change to seInd
		for( Ordinal i=std::max(space()->sInd(U,3),0)+1; i<=space()->eInd(U,3); ++i ) {

      op_->apply( wind0_->get(),      y.getConstSField(i-1), z.getSField(i-1), 0., 0., 1., iRe );

      op_->apply( winds_[i-1]->get(), y.getConst0Field(),    z.getSField(i-1), 1., 0., 1., 0.  );

      for( int k=1; k+i<=Nf; ++k ) { // that is fine
				mulI = (k-1==i-1)?(a2*i):0;
        op_->apply( windc_[k+i-1]->get(), y.getConstSField(k-1),   z.getSField(i-1), 1.,    0., -0.5, 0. );

        op_->apply( windc_[k-1]->get(),   y.getConstSField(k+i-1), z.getSField(i-1), 1.,    0.,  0.5, 0. );

        op_->apply( winds_[k+i-1]->get(), y.getConstCField(k-1),   z.getSField(i-1), 1., -mulI,  0.5, 0. );

        op_->apply( winds_[k-1]->get(),   y.getConstCField(k+i-1), z.getSField(i-1), 1.,    0., -0.5, 0. );
      }
    }

		// rest of time
//		for( int i=Nf/2+1; i<=Nf; ++i ) { // change to local
		for( Ordinal i=std::max(space()->sInd(U,3),0)+1; i<=space()->eInd(U,3); ++i ) {
			if( Nf/2+1<=i && i<=Nf ) {
				mulI = a2*i;
				z.getCFieldPtr(i-1)->add( 1., z.getCField(i-1),  mulI, y.getConstSField(i-1) );
				z.getSFieldPtr(i-1)->add( 1., z.getSField(i-1), -mulI, y.getConstCField(i-1) );
			}
		}


		// strange terms
		int i;
		for( int k=1; k<=Nf; ++k ) { 
			for( int l=1; l<=Nf; ++l ) { // that is fine
				i = k+l; 
				if( i<=Nf ) { // do something here
					if( std::max(space()->sInd(U,3),0)+1<=i && i<=space()->eInd(U,3)) {
						op_->apply( windc_[k-1]->get(), y.getConstCField(l-1), z.getCField(i-1), 1., 0.,  0.5, 0. );
						op_->apply( winds_[k-1]->get(), y.getConstSField(l-1), z.getCField(i-1), 1., 0., -0.5, 0. );

						op_->apply( windc_[k-1]->get(), y.getConstSField(l-1), z.getSField(i-1), 1., 0., 0.5, 0. );
						op_->apply( winds_[k-1]->get(), y.getConstCField(l-1), z.getSField(i-1), 1., 0., 0.5, 0. );
					}
				}
			}
		}

		z.changed();

  }



//  void apply(const DomainFieldT& x, RangeFieldT& y) const {
//    apply(x,x,y);
//  }


	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(op_->space()); };

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
