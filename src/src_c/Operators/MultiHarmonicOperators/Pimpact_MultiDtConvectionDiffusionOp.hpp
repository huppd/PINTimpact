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

  using DomainFieldT = MultiHarmonicField< VectorField<SpaceT> >;
  using RangeFieldT = MultiHarmonicField< VectorField<SpaceT> >;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  Teuchos::RCP<NonlinearWrap< ConvectionDiffusionSOp<SpaceT> > > op_;

  Teuchos::RCP< ConvectionField<SpaceT> > wind0_;
  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > windc_;
  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > winds_;

  using FieldTensorT = typename ConvectionField<SpaceT>::FieldTensor;

public:

  MultiDtConvectionDiffusionOp( const Teuchos::RCP<const SpaceT>& space ):
    op_( create<NonlinearWrap>( create<ConvectionDiffusionSOp<SpaceT> >(space) ) ),
    wind0_( create<ConvectionField>(space) ),
    windc_( space->nGlo(3) ),
    winds_( space->nGlo(3) ) {

    for( Ordinal i=0; i<space->nGlo(3); ++i ) {
      windc_[i] = create<ConvectionField>( space );
      winds_[i] = create<ConvectionField>( space );
    }

  };


	/// \todo check for mv global if not make global temp + exchange
  void assignField( const DomainFieldT& mv ) {

		mv.exchange();

    wind0_->assignField( mv.getConst0Field() );

    for( Ordinal i=1; i<=space()->nGlo(3); ++i ) {
      windc_[i-1]->assignField( mv.getConstCField(i) );
      winds_[i-1]->assignField( mv.getConstSField(i) );
    }

  };


	/// \todo check for y global if not make global temp + exchange
  void apply( const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {
		
		Ordinal Nf = space()->nGlo(3);
		Scalar iRe = 1./op_->space()->getDomainSize()->getRe();
		Scalar a2 = op_->space()->getDomainSize()->getAlpha2()*iRe;

		Scalar mulI;


		y.exchange();

    // computing zero mode of z
		if( 0==space()->sInd(U,3) ) {

			op_->apply( get0Wind(), y.getConst0Field(), z.get0Field(), 0., 0., 1., iRe );

			for( Ordinal i=1; i<=Nf; ++i ) {
				op_->apply( getCWind(i), y.getConstCField(i), z.get0Field(), 1., 0., 0.5, 0. );
				op_->apply( getSWind(i), y.getConstSField(i), z.get0Field(), 1., 0., 0.5, 0. );
			}
		}


    // computing cos mode of z
		for( Ordinal i=std::max(space()->sInd(U,3),1); i<=space()->eInd(U,3); ++i ) {

      op_->apply( get0Wind( ), y.getConstCField(i), z.getCField(i), 0., 0., 1., iRe );
      op_->apply( getCWind(i), y.getConst0Field( ), z.getCField(i), 1., 0., 1., 0.  );

      for( Ordinal k=1; k+i<=Nf; ++k ) { // thats fine

				mulI = (k==i)?(a2*i):0;

        op_->apply( getCWind(k+i), y.getConstCField( k ), z.getCField(i), 1.,   0., 0.5, 0. );
        op_->apply( getCWind( k ), y.getConstCField(k+i), z.getCField(i), 1.,   0., 0.5, 0. );
        op_->apply( getSWind(k+i), y.getConstSField( k ), z.getCField(i), 1., mulI, 0.5, 0. );
        op_->apply( getSWind( k ), y.getConstSField(k+i), z.getCField(i), 1.,   0., 0.5, 0. );
      }
    }

    // computing sin mode of y
		for( Ordinal i=std::max(space()->sInd(U,3),1); i<=space()->eInd(U,3); ++i ) {

      op_->apply( get0Wind(),  y.getConstSField(i), z.getSField(i), 0., 0., 1., iRe );
      op_->apply( getSWind(i), y.getConst0Field(),  z.getSField(i), 1., 0., 1., 0.  );

      for( Ordinal k=1; k+i<=Nf; ++k ) { // that is fine

				mulI = (k==i)?(a2*i):0;

        op_->apply( getCWind(k+i), y.getConstSField( k ), z.getSField(i), 1.,    0., -0.5, 0. );
        op_->apply( getCWind( k ), y.getConstSField(k+i), z.getSField(i), 1.,    0.,  0.5, 0. );
        op_->apply( getSWind(k+i), y.getConstCField( k ), z.getSField(i), 1., -mulI,  0.5, 0. );
        op_->apply( getSWind( k ), y.getConstCField(k+i), z.getSField(i), 1.,    0., -0.5, 0. );
      }
    }

		// rest of time
		for( Ordinal i=std::max(space()->sInd(U,3),1); i<=space()->eInd(U,3); ++i ) {
			if( Nf/2+1<=i && i<=Nf ) {
				mulI = a2*i;
				z.getCFieldPtr(i)->add( 1., z.getCField(i),  mulI, y.getConstSField(i) );
				z.getSFieldPtr(i)->add( 1., z.getSField(i), -mulI, y.getConstCField(i) );
			}
		}


		// strange terms
		Ordinal i;
		for( Ordinal k=1; k<=Nf; ++k ) { 
			for( Ordinal l=1; l<=Nf; ++l ) { // that is fine
				i = k+l; 
				if( i<=Nf ) { // do something here
					if( std::max(space()->sInd(U,3),1)<=i && i<=space()->eInd(U,3) ) {
						op_->apply( getCWind(k), y.getConstCField(l), z.getCField(i), 1., 0.,  0.5, 0. );
						op_->apply( getSWind(k), y.getConstSField(l), z.getCField(i), 1., 0., -0.5, 0. );

						op_->apply( getCWind(k), y.getConstSField(l), z.getSField(i), 1., 0.,  0.5, 0. );
						op_->apply( getSWind(k), y.getConstCField(l), z.getSField(i), 1., 0.,  0.5, 0. );
					}
				}
			}
		}

		z.changed();
  }


	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "MHDtConvectionDiffusion" ); };

  void print( std::ostream& out=std::cout ) const {
		out <<  getLabel() << ":\n";
		op_->print( out );
  }

protected:

	constexpr const FieldTensorT& get0Wind() const { return( wind0_->get() ); }

	constexpr const FieldTensorT& getCWind( const Ordinal& i) const { return( windc_[i-1]->get() ); }
	constexpr const FieldTensorT& getSWind( const Ordinal& i) const { return( winds_[i-1]->get() ); }

}; // end of class MultiDtConvectionDiffusionOp



/// \relates MultiDtConvectionDiffusionOp
template<class SpaceT>
Teuchos::RCP<MultiDtConvectionDiffusionOp<SpaceT> >
createMultiDtConvectionDiffusionOp( const Teuchos::RCP<const SpaceT>& space ) {

  return( Teuchos::rcp( new MultiDtConvectionDiffusionOp<SpaceT>( space ) ) );

}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULTIDTCONVECTIONDIFFUSIONOP_HPP
