#pragma once
#ifndef PIMPACT_TIMENSOp_HPP
#define PIMPACT_TIMENSOp_HPP


#include "Teuchos_RCP.hpp"
#include "Pimpact_ScalarField.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_ConvectionSOp.hpp"




namespace Pimpact {

extern "C" {

void OP_TimeNS( 
		const int dimens,
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

	using SpaceT = ST;

	using DomainFieldT = CompoundField< TimeField<VectorField<ST> >, TimeField<ScalarField<ST> > >;
	using RangeFieldT = CompoundField< TimeField<VectorField<ST> >, TimeField<ScalarField<ST> > >;

protected:

	using Scalar = typename SpaceT::Scalar;
	using Ordinal = typename SpaceT::Ordinal;

	static const int sdim = SpaceT::sdim;
	static const int dimension = SpaceT::dimension;

	static const int dimNC = SpaceT::dimNC;

	Teuchos::RCP< TimeField<VectorField<ST> > > windU_;
	Teuchos::RCP< TimeField<VectorField<ST> > > windV_;
	Teuchos::RCP< TimeField<VectorField<ST> > > windW_;

	Teuchos::RCP<const InterpolateS2V<ST> > interpolateS2V_;
	Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,sdim,dimension,dimNC> > interpolateV2S_;
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
		Scalar re = space()->getDomainSize()->getRe();
		Scalar mulI = space()->getDomainSize()->getAlpha2()*idt/re;

		auto& xu = x.getVField();
		auto& xp = x.getSField();
		auto& yu = y.getVField();
		auto& yp = y.getSField();

		xu.exchange();

		for( Ordinal i=space()->si(F::S,3)-1; i<space()->ei(F::S,3); ++i ) {
			xu(i).exchange();
			xp(i).exchange();
		}

		OP_TimeNS( 
				SpaceT::sdim,
				space()->nLoc(),
				space()->bl(),
				space()->bu(),
				space()->nl(),
				space()->nu(),
				space()->dl(),
				space()->du(),
				space()->gl(),
				space()->gu(),
				space()->sInd(F::S),
				space()->eInd(F::S),
				space()->sInd(F::U),
				space()->eInd(F::U),
				space()->sInd(F::V),
				space()->eInd(F::V),
				space()->sInd(F::W),
				space()->eInd(F::W),
				conv_->getCD(X,F::U),
				conv_->getCD(Y,F::V),
				conv_->getCD(Z,F::W),
				conv_->getCU(X,F::U),
				conv_->getCU(Y,F::V),
				conv_->getCU(Z,F::W),
				conv_->getCD(X,F::S),
				conv_->getCD(Y,F::S),
				conv_->getCD(Z,F::S),
				conv_->getCU(X,F::S),
				conv_->getCU(Y,F::S),
				conv_->getCU(Z,F::S),
				helm_->getC(X,F::S),
				helm_->getC(Y,F::S),
				helm_->getC(Z,F::S),
				helm_->getC(X,F::U),
				helm_->getC(Y,F::V),
				helm_->getC(Z,F::W),
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
				xu.getConstRawPtr(),
				xp.getConstRawPtr(),
				yu.getRawPtr(),
				yp.getRawPtr() );

		for( Ordinal i=space()->si(F::S,3)-1; i<space()->ei(F::S,3); ++i ) {
			yu(i).changed();
			yp(i).changed();
		}

		yu.changed();
		yp.changed();

	}


	void computeResidual( const RangeFieldT& b, const DomainFieldT& x, RangeFieldT& res ) const {
		apply( x, res );
		res.add( 1., b, -1., res );
	}

	void assignField( const DomainFieldT& cmv ) {

		Ordinal nt = space()->nLoc(3) + space()->bu(3) - space()->bl(3);

		auto& mv = cmv.getVField();

		ScalarField<ST> temp( space() );

		mv.exchange();

		for( Ordinal it=0; it<nt; ++it ) {

			interpolateV2S_->apply( mv(it)(F::U), temp );
			for( F j=F::U; j<SpaceT::sdim; ++j ) {
				interpolateS2V_->apply( temp, windU_->operator()(it)(j) );
			}
			interpolateV2S_->apply( mv(it)(F::V), temp );
			for( F j=F::U; j<SpaceT::sdim; ++j ) {
				interpolateS2V_->apply( temp, windV_->operator()(it)(j) );
			}
			if( 3==SpaceT::sdim ) {
				interpolateV2S_->apply( mv(it)(F::W), temp );
				for( F j=F::U; j<SpaceT::sdim; ++j ) {
					interpolateS2V_->apply( temp, windW_->operator()(it)(j) );
				}
			}
		}
		windU_->changed();
		windV_->changed();
		if( 3==SpaceT::sdim )
			windW_->changed();

	};

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return( conv_->space() ); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

	bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const HelmholtzOp<ST> > getHelmholtzOp() const { return( helm_ ); }
	Teuchos::RCP<const GradOp<ST> > getGradOp() const { return( grad_ ); }
	Teuchos::RCP<const DivOp<ST> > getDivOp() const { return( div_ ); }
	Teuchos::RCP<const ConvectionSOp<ST> > getConvOp() const { return( conv_ ); }

	Teuchos::RCP< TimeField<VectorField<ST> > > getWindU_() const { return( windU_ ); }
	Teuchos::RCP< TimeField<VectorField<ST> > > getWindV_() const { return( windV_ ); }
	Teuchos::RCP< TimeField<VectorField<ST> > > getWindW_() const { return( windW_ ); }



	void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
	}

	const std::string getLabel() const { return( "TimeNSOp" ); };

}; // end of class TimeNSOp



} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::TimeNSOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::TimeNSOp< Pimpact::Space<double,int,4,4> >;
#endif



#endif // end of #ifndef PIMPACT_TIMENSOp_HPP
