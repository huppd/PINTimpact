#pragma once
#ifndef PIMPACT_TIMESTOKESOP_HPP
#define PIMPACT_TIMESTOKESOP_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_CompoundField.hpp"




namespace Pimpact {

extern "C" {

	void OP_TimeStokes( 
      const int& dimens,
      const int* const N,
      const int* const bl,
      const int* const bu,
      const int* const dl,
      const int* const du,
      const int* const gl,
      const int* const gu,
      const int* const ss,
      const int* const nn,
      const int* const su,
      const int* const nu,
      const int* const sv,
      const int* const nv,
      const int* const sw,
      const int* const nw,
      const double* const c11p,
			const double* const c22p,
			const double* const c33p,       
      const double* const c11u,
			const double* const c22v,
			const double* const c33w,       
      const double* const cD1,                  
      const double* const cD2,                  
      const double* const cD3,                  
      const double* const cG1,                  
      const double* const cG2,                  
      const double* const cG3,                  
      const double& mulI,                 
      const double& mulL,                 
      const double* const velp,                 
      const double* const veln,                 
      const double* const pn,                   
            double* const r_vel,                
            double* const r_p );

}

/// \ingroup CompoundOperator
///
// \f[ \begin{bmatrix} opV2V & opS2V \\ opV2S & 0 \end{bmatrix} \mathbf{x} = \mathbf{y} \f]
template<class ST>
class TimeStokesOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef CompoundField< TimeField<VectorField<ST> >, TimeField<ScalarField<ST> > >  DomainFieldT;
  typedef CompoundField< TimeField<VectorField<ST> >, TimeField<ScalarField<ST> > >  RangeFieldT;


protected:

//  Teuchos::RCP<VF> temp_;

  Teuchos::RCP<const HelmholtzOp<ST> > helm_;
  Teuchos::RCP<const GradOp<ST> > grad_;
  Teuchos::RCP<const DivOp<ST> > div_;

public:

  /// \todo constructor from space
	TimeStokesOp(
			const Teuchos::RCP<const SpaceT>& space ):
		helm_( create< HelmholtzOp<ST> >(space) ),
		grad_( create< GradOp<ST> >(space) ),
		div_( create< DivOp<ST> >(space) ) {};

	void apply(const DomainFieldT& x, RangeFieldT& y,
			Belos::ETrans trans=Belos::NOTRANS  ) const {

		Scalar pi = 4.*std::atan(1.);
		Scalar idt = ((Scalar)space()->nGlo()[3])/2./pi;
		Scalar re = space()->getDomain()->getDomainSize()->getRe();
		Scalar mulI = space()->getDomain()->getDomainSize()->getAlpha2()*idt/re;

		auto xu = x.getConstVFieldPtr();
		auto xp = x.getConstSFieldPtr();
		auto yu = y.getVFieldPtr();
		auto yp = y.getSFieldPtr();

		xu->exchange();
//		xp->exchange();

		for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i ) {
			xu->getConstFieldPtr(i-1)->exchange();
			xu->getConstFieldPtr(i)->exchange();
			xp->getConstFieldPtr(i)->exchange();

			OP_TimeStokes( 
        space()->dim(),
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
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
				xu->getConstFieldPtr(i-1)->getConstRawPtr(),
				xu->getConstFieldPtr(i)->getConstRawPtr(),
				xp->getConstFieldPtr(i)->getConstRawPtr(),
				yu->getFieldPtr(i)->getRawPtr(),
				yp->getFieldPtr(i)->getRawPtr() );

			yu->getFieldPtr(i)->changed();
			yp->getFieldPtr(i)->changed();

		}

		yu->changed();
		yp->changed();
		
	}

	void assignField( const DomainFieldT& mv ) { };

	Teuchos::RCP<const SpaceT> space() const { return( helm_->space() ); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }


  Teuchos::RCP<const HelmholtzOp<ST> > getHelmholtzOp() const { return( helm_ ); }
  Teuchos::RCP<const GradOp<ST> > getGradOp() const { return( grad_ ); }
  Teuchos::RCP<const DivOp<ST> > getDivOp() const { return( div_ ); }
}; // end of class TimeStokesOp



} // end of namespace Pimpact



#ifdef COMPILE_ETI
extern template class Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,4> >;
#endif

#endif // end of #ifndef PIMPACT_TIMESTOKESOP_HPP