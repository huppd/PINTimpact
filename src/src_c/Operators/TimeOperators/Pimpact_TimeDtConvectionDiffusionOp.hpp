#pragma once
#ifndef PIMPACT_TIMEDTCONVECTIONDIFFUSIONOP_HPP
#define PIMPACT_TIMEDTCONVECTIONDIFFUSIONOP_HPP


#include "Pimpact_Types.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_TimeField.hpp"

#include "Pimpact_ConvectionDiffusionSOp.hpp"
#include "Pimpact_ConvectionVWrap.hpp"



namespace Pimpact {



/// \ingroup TimeHarmonicOperator
/// \ingroup NonliearOperator
template<class ST, bool CNyes=false>
class TimeDtConvectionDiffusionOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef TimeField< VectorField<SpaceT> >  DomainFieldT;
  typedef TimeField< VectorField<SpaceT> >  RangeFieldT;

protected:

  Teuchos::RCP<ConvectionVWrap< ConvectionDiffusionSOp<SpaceT> > > op_;

  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > wind_;

public:

  /// \todo get nf from grid
  TimeDtConvectionDiffusionOp( const Teuchos::RCP<const SpaceT>& space ):
    op_( create<ConvectionVWrap>( create<ConvectionDiffusionSOp<SpaceT> >(space) ) ),
    wind_( space()->nLoc(3) + space()->bu(3) - space()->bl(3) ) {

    Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

    for( Ordinal i=0; i<nt; ++i ) 
      wind_[i] = create<ConvectionField>( space );

  };

  void assignField( const DomainFieldT& mv ) {

		mv.exchange();

    Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

    for( Ordinal i=0; i<nt; ++i ) {
      wind_[i]->assignField( mv.getConstField(i) );
    }

  };


  void apply( const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {
		
		Ordinal sInd = space()->sInd(EField::S,3);
		Ordinal eInd = space()->eInd(EField::S,3);

    Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

		Scalar iRe = 1./space()->getDomain()->getDomainSize()->getRe();
		Scalar a2 = space()->getDomain()->getDomainSize()->getAlpha2()*iRe;

		Scalar pi = 4.*std::atan(1.);

		Scalar mulI = a2*iRe*((Scalar)space()->nGlo(3))/2./pi;

		
		y.exchange();

		if( CNyes ) {

			auto temp = create<VectorField>( space() );


			for( Ordinal i=sInd; i<eInd; ++i ) {

				z.getField(i).add( mulI, y.getConstField(i), -mulI, y.getConstField(i-1) );

				if( i>sInd )
					z.getField(i).add( 1., z.getConstField(i), 0.5, *temp );

				op_->apply( wind_[i]->get(), y.getConstField(i), *temp, 0., 0., 1., iRe );

				z.getField(i).add( 1., z.getConstField(i), 0.5, *temp );

			}

			op_->apply( wind_[sInd-1]->get(), y.getConstField(sInd-1), z.getField(sInd), 1., 0., 0.5, iRe*0.5 );


		}
		else {
			for( Ordinal i=sInd; i<eInd; ++i ) {
				op_->apply( wind_[i]->get(), y.getConstField(i), z.getField(i), 0., mulI, 1., iRe );
				z.getField(i).add( 1., z.getConstField(i), -mulI, y.getConstField(i-1) );
			}
		}
		z.changed();

	}



	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

}; // end of class TimeDtConvectionDiffusionOp



///// \relates TimeDtConvectionDiffusionOp
//template<class SpaceT>
//Teuchos::RCP<TimeDtConvectionDiffusionOp<SpaceT> >
//createTimeDtConvectionDiffusionOp( const Teuchos::RCP<const SpaceT>& space ) {
//
//  return( Teuchos::rcp( new TimeDtConvectionDiffusionOp<SpaceT>( space ) ) );
//
//}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_TIMEDTCONVECTIONDIFFUSIONOP_HPP