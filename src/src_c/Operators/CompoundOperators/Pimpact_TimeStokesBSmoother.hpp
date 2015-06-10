#pragma once
#ifndef PIMPACT_TIMESTOKESBSMOOTHER_HPP
#define PIMPACT_TIMESTOKESBSMOOTHER_HPP


#include "Teuchos_RCP.hpp"




namespace Pimpact {

extern "C" {

	void OP_TimeStokesBSmoother( 
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
template<class OperatorT>
class TimeStokesBSmoother {

public:

  typedef typename OperatorT::SpaceT SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef CompoundField< TimeField<VectorField<SpaceT> >, TimeField<ScalarField<SpaceT> > >  DomainFieldT;
  typedef CompoundField< TimeField<VectorField<SpaceT> >, TimeField<ScalarField<SpaceT> > >  RangeFieldT;


protected:

  const Teuchos::RCP<const OperatorT> op_;

public:

  /// \todo constructor from space
	TimeStokesBSmoother( const Teuchos::RCP<const OperatorT>& op ):
		op_( op ) {};

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
				op_->getHelmholtzOp()->getC(X,S),
				op_->getHelmholtzOp()->getC(Y,S),
				op_->getHelmholtzOp()->getC(Z,S),
				op_->getHelmholtzOp()->getC(X,U),
				op_->getHelmholtzOp()->getC(Y,V),
				op_->getHelmholtzOp()->getC(Z,W),
				op_->getDivOp()->getC(X),
				op_->getDivOp()->getC(Y),
				op_->getDivOp()->getC(Z),
				op_->getGradOp()->getC(X),
				op_->getGradOp()->getC(Y),
				op_->getGradOp()->getC(Z),
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

	Teuchos::RCP<const SpaceT> space() const { return( op_->space() ); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }


}; // end of class TimeStokesBSmoother



} // end of namespace Pimpact



#ifdef COMPILE_ETI
#include "Pimpact_TimeStokesOp.hpp"
extern template class Pimpact::TimeStokesBSmoother< Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::TimeStokesBSmoother< Pimpact::TimeStokesOp< Pimpact::Space<double,int,4,4> > >;
#endif

#endif // end of #ifndef PIMPACT_TIMESTOKESBSMOOTHER_HPP