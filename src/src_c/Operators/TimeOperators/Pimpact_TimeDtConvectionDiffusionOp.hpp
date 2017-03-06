#pragma once
#ifndef PIMPACT_TIMEDTCONVECTIONDIFFUSIONOP_HPP
#define PIMPACT_TIMEDTCONVECTIONDIFFUSIONOP_HPP


#include "Pimpact_ConvectionDiffusionSOp.hpp"
#include "Pimpact_NonlinearVWrap.hpp"
#include "Pimpact_TimeField.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {



/// \ingroup TimeHarmonicOperator
/// \ingroup NonliearOperator
/// \deprecated
template<class ST, bool CNyes=false>
class TimeDtConvectionDiffusionOp {

public:

  using SpaceT = ST;

  using DomainFieldT = TimeField< VectorField<SpaceT> >;
  using RangeFieldT = TimeField< VectorField<SpaceT> >;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  Teuchos::RCP<NonlinearWrap< ConvectionDiffusionSOp<SpaceT> > > op_;

  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > wind_;

public:


  TimeDtConvectionDiffusionOp( const Teuchos::RCP<const SpaceT>& space ):
    op_( create<NonlinearWrap>( create<ConvectionDiffusionSOp<SpaceT> >(space) ) ),
    wind_( space()->nLoc(3) + space()->bu(3) - space()->bl(3) ) {

    Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

    for( Ordinal i=0; i<nt; ++i ) 
      wind_[i] = create<ConvectionField>( space );
  };

  void assignField( const DomainFieldT& mv ) {

		mv.exchange();

    Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

    for( Ordinal i=0; i<nt; ++i ) {
      wind_[i]->assignField( mv.getField(i) );
    }
  };


  void apply( const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {
		
		Ordinal sInd = space()->begin(F::S,3);
		Ordinal eInd = space()->end(F::S,3);

    Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

		Scalar iRe = 1./space()->getDomainSize()->getRe();
		Scalar a2 = space()->getDomainSize()->getAlpha2()*iRe;

		Scalar pi = 4.*std::atan(1.);

		Scalar mulI = a2*iRe*(static_cast<Scalar>(space()->nGlo(3)))/2./pi;

		
		y.exchange();

		if( CNyes ) {

			VectorField<SpaceT> temp( space() );


			for( Ordinal i=sInd; i<=eInd; ++i ) {

				z.getField(i).add( mulI, y.getField(i), -mulI, y.getField(i-1), B::N );

				if( i>sInd )
					z.getField(i).add( 1., z.getField(i), 0.5, temp );

				op_->apply( wind_[i]->get(), y.getField(i), temp, 0., 1., iRe, Add::N );

				z.getField(i).add( 1., z.getField(i), 0.5, temp );
			}

			op_->apply( wind_[sInd-1]->get(), y.getField(sInd-1), z.getField(sInd), 0., 0.5, iRe*0.5, Add::Y );

		}
		else {
			for( Ordinal i=sInd; i<=eInd; ++i ) {
				op_->apply( wind_[i]->get(), y.getField(i), z.getField(i), mulI, 1., iRe, Add::N );
				z.getField(i).add( 1., z.getField(i), -mulI, y.getField(i-1), B::N );
			}
		}
		z.changed();
	}



	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "TimeDtConvectionDiffusionOp " ); };

  void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
		op_->print( out );
  }

}; // end of class TimeDtConvectionDiffusionOp



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_TIMEDTCONVECTIONDIFFUSIONOP_HPP
