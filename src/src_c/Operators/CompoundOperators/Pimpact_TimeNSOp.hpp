#pragma once
#ifndef PIMPACT_TIMENSOp_HPP
#define PIMPACT_TIMENSOp_HPP


#include "Teuchos_RCP.hpp"
#include "Pimpact_ScalarField.hpp"

#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_ConvectionSOp.hpp"
#include "Pimpact_CompoundField.hpp"




namespace Pimpact {

extern "C" {

	void OP_TimeNS( 
      const int& dimens,
      const int* const N,
			const int* const bl, const int* const bu,
			const int* const cL, const int* const cU,
			const int* const dl, const int* const du,
			const int* const gl, const int* const gu,
			const int* const ss, const int* const nn,
			const int* const su, const int* const nu,
			const int* const sv, const int* const nv,
			const int* const sw, const int* const nw,
			const double* const c1uD, const double* const c2vD, const double* const c3wD,
			const double* const c1uU, const double* const c2vU,	const double* const c3wU,
			const double* const c1pD,	const double* const c2pD,	const double* const c3pD,
			const double* const c1pU,	const double* const c2pU,	const double* const c3pU,
      const double* const c11p,	const double* const c22p,	const double* const c33p,       
      const double* const c11u,	const double* const c22v,	const double* const c33w,       
      const double* const cD1,                  
      const double* const cD2,                  
      const double* const cD3,                  
      const double* const cG1,                  
      const double* const cG2,                  
      const double* const cG3,                  
      const double& mulI,                 
      const double& mulL,                 
      const double* const windU,                 
      const double* const windV,                 
      const double* const windW,                 
      const double* const veln,                 
      const double* const pn,                   
            double* const r_vel,                
            double* const r_p );

}

/// \ingroup CompoundOperator
///
// \f[ \begin{bmatrix} opV2V & opS2V \\ opV2S & 0 \end{bmatrix} \mathbf{x} = \mathbf{y} \f]
template<class ST>
class TimeNSOp {

public:

  typedef ST SpaceT;

  typedef CompoundField< TimeField<VectorField<ST> >, TimeField<ScalarField<ST> > >  DomainFieldT;
  typedef CompoundField< TimeField<VectorField<ST> >, TimeField<ScalarField<ST> > >  RangeFieldT;


protected:

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  static const int dimension = SpaceT::dimension;

  static const int dimNC = SpaceT::dimNC;


	Teuchos::RCP< TimeField<VectorField<ST> > > windU_;
	Teuchos::RCP< TimeField<VectorField<ST> > > windV_;
	Teuchos::RCP< TimeField<VectorField<ST> > > windW_;

  Teuchos::RCP<const InterpolateS2V<ST> > interpolateS2V_;
  Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension,dimNC> > interpolateV2S_;
  Teuchos::RCP<const ConvectionSOp<ST> > conv_;
  Teuchos::RCP<const HelmholtzOp<ST>   > helm_;
  Teuchos::RCP<const GradOp<ST>        > grad_;
  Teuchos::RCP<const DivOp<ST>         > div_;

public:


	TimeNSOp( const Teuchos::RCP<const SpaceT>& space ):
		windU_( create< TimeField< VectorField<ST> > >( space) ),
		windV_( create< TimeField< VectorField<ST> > >( space) ),
		windW_( create< TimeField< VectorField<ST> > >( space) ),
    interpolateS2V_( create<InterpolateS2V>(space) ),
    interpolateV2S_( createInterpolateV2S( space ) ),
		conv_( create< ConvectionSOp<ST> >(space) ),
		helm_( create< HelmholtzOp<ST>   >(space) ),
		grad_( create< GradOp<ST>        >(space) ),
		div_ ( create< DivOp<ST>         >(space) ) {};


	void apply( const DomainFieldT& x, RangeFieldT& y ) const { 

		Scalar pi = 4.*std::atan(1.);
		Scalar idt = ((Scalar)space()->nGlo()[3])/2./pi;
		Scalar re = space()->getDomain()->getDomainSize()->getRe();
		Scalar mulI = space()->getDomain()->getDomainSize()->getAlpha2()*idt/re;

		auto xu = x.getConstVFieldPtr();
		auto xp = x.getConstSFieldPtr();
		auto yu = y.getVFieldPtr();
		auto yp = y.getSFieldPtr();

		xu->exchange();

		for( Ordinal i=space()->sInd(S,3)-1; i<space()->eInd(S,3); ++i ) {
			xu->getConstFieldPtr(i)->exchange();
			xp->getConstFieldPtr(i)->exchange();
		}

		OP_TimeNS( 
				space()->dim(),
				space()->nLoc(),
				space()->bl(),
				space()->bu(),
				space()->nl(),
				space()->nu(),
				space()->dl(),
				space()->du(),
				space()->gl(),
				space()->gu(),
				space()->sInd(S),
				space()->eInd(S),
				space()->sInd(U),
				space()->eInd(U),
				space()->sInd(V),
				space()->eInd(V),
				space()->sInd(W),
				space()->eInd(W),
        conv_->getCD(X,U),
        conv_->getCD(Y,V),
        conv_->getCD(Z,W),
        conv_->getCU(X,U),
        conv_->getCU(Y,V),
        conv_->getCU(Z,W),
        conv_->getCD(X,S),
        conv_->getCD(Y,S),
        conv_->getCD(Z,S),
        conv_->getCU(X,S),
        conv_->getCU(Y,S),
        conv_->getCU(Z,S),
				helm_->getC(X,S),
				helm_->getC(Y,S),
				helm_->getC(Z,S),
				helm_->getC(X,U),
				helm_->getC(Y,V),
				helm_->getC(Z,W),
				div_->getC(X),
				div_->getC(Y),
				div_->getC(Z),
				grad_->getC(X),
				grad_->getC(Y),
				grad_->getC(Z),
				mulI,                 
				1./re,                 
				windU_->getConstRawPtr(),
				windV_->getConstRawPtr(),
				windW_->getConstRawPtr(),
				xu->getConstRawPtr(),
				xp->getConstRawPtr(),
				yu->getRawPtr(),
				yp->getRawPtr() );

		for( Ordinal i=space()->sInd(S,3)-1; i<space()->eInd(S,3); ++i ) {
			yu->getFieldPtr(i)->changed();
			yp->getFieldPtr(i)->changed();
		}

		yu->changed();
		yp->changed();

	}


	void assignField( const DomainFieldT& cmv ) {

    Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

		auto mv = cmv.getConstVFieldPtr();

    auto temp = createScalarField<ST>( space() );

		mv->exchange();

		for( Ordinal it=0; it<nt; ++it ) {

			interpolateV2S_->apply( mv->getConstField(it).getConstField(U), *temp );
			for( int j=0; j<space()->dim(); ++j ) {
				interpolateS2V_->apply( *temp, windU_->getFieldPtr(it)->getField(j) );
			}
			interpolateV2S_->apply( mv->getConstField(it).getConstField(V), *temp );
			for( int j=0; j<space()->dim(); ++j ) {
				interpolateS2V_->apply( *temp, windV_->getFieldPtr(it)->getField(j) );
			}
			if( 3==space()->dim() ) {
				interpolateV2S_->apply( mv->getConstField(it).getConstField(W), *temp );
				for( int j=0; j<space()->dim(); ++j ) {
					interpolateS2V_->apply( *temp, windW_->getFieldPtr(it)->getField(j) );
				}
			}
		}
		windU_->changed();
		windV_->changed();
		if( 3==space()->dim() )
			windW_->changed();

  };

	Teuchos::RCP<const SpaceT> space() const { return( conv_->space() ); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "TimeNSOp" ); };

}; // end of class TimeNSOp



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_TIMENSOp_HPP
