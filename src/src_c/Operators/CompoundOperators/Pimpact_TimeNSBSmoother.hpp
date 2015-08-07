#pragma once
#ifndef PIMPACT_TIMENSBSMOOTHER_HPP
#define PIMPACT_TIMENSBSMOOTHER_HPP


#include "Teuchos_RCP.hpp"


namespace Pimpact {

extern "C" {

	void OP_TimeNSBSmoother( 
      const int& dimens,
      const int* const N,
      const int* const bl,
      const int* const bu,
      const int* const cl,
      const int* const cu,
      //const int* const BCL,
      //const int* const BCU,
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
    const double* const c1uD,
    const double* const c2vD,
    const double* const c3wD,
    const double* const c1uU,
    const double* const c2vU,
    const double* const c3wU,
    const double* const c1pD,
    const double* const c2pD,
    const double* const c3pD,
    const double* const c1pU,
    const double* const c2pU,
    const double* const c3pU,
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
      const double* const rhs_vel,                 
      const double* const rhs_p,                   
            double* const vel,                
            double* const p );

}

/// \ingroup CompoundOperator
///
// \f[ \begin{bmatrix} opV2V & opS2V \\ opV2S & 0 \end{bmatrix} \mathbf{x} = \mathbf{y} \f]
template<class OperatorT>
class TimeNSBSmoother {

public:

  typedef typename OperatorT::SpaceT SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef CompoundField< TimeField<VectorField<SpaceT> >, TimeField<ScalarField<SpaceT> > >  DomainFieldT;
  typedef CompoundField< TimeField<VectorField<SpaceT> >, TimeField<ScalarField<SpaceT> > >  RangeFieldT;


protected:

  const Teuchos::RCP<const OperatorT> op_;
	int numIters_;

public:

  /// \todo constructor from space
	TimeNSBSmoother( const Teuchos::RCP<const OperatorT>& op , Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList()):
		op_( op ),
		numIters_( pl->get<int>("numIters",4) )	{};

	void apply(const DomainFieldT& x, RangeFieldT& y ) const {

		Scalar pi = 4.*std::atan(1.);
		Scalar idt = ((Scalar)space()->nGlo()[3])/2./pi;
		Scalar re = space()->getDomain()->getDomainSize()->getRe();
		Scalar mulI = space()->getDomain()->getDomainSize()->getAlpha2()*idt/re;

		for( int iters=0; iters<numIters_; ++iters ) {

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
			}

			OP_TimeNSBSmoother(
					space()->dim(),
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					space()->nl(),
					space()->nu(),
                    //space()->getDomain()->getBCLocal()->getBCL(),
                    //space()->getDomain()->getBCLocal()->getBCU(),
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
                      			op_->getConvOp()->getCD(X,U),
                      			op_->getConvOp()->getCD(Y,V),
                      			op_->getConvOp()->getCD(Z,W),
                      			op_->getConvOp()->getCU(X,U),
                      			op_->getConvOp()->getCU(Y,V),
                      			op_->getConvOp()->getCU(Z,W),
                      			op_->getConvOp()->getCD(X,S),
                      			op_->getConvOp()->getCD(Y,S),
                      			op_->getConvOp()->getCD(Z,S),
                      			op_->getConvOp()->getCU(X,S),
                      			op_->getConvOp()->getCU(Y,S),
                      			op_->getConvOp()->getCU(Z,S),
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
					windU_->getConstRawPtr(),
					windV_->getConstRawPtr(),
					windW_->getConstRawPtr(),               
					xu->getConstRawPtr(),
					xp->getConstRawPtr(),
					yu->getRawPtr(),
					yp->getRawPtr() );

			for( Ordinal i=space()->sInd(S,3); i<space()->eInd(S,3); ++i ) {
				yu->getFieldPtr(i)->changed();
				yp->getFieldPtr(i)->changed();
			} 

			yu->changed();
			yp->changed();
		}
	}

	void assignField( const DomainFieldT& mv ) { };

	Teuchos::RCP<const SpaceT> space() const { return( op_->space() ); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

    bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "TimeNSBSmoother " ); };

}; // end of class TimeNSBSmoother



} // end of namespace Pimpact



#ifdef COMPILE_ETI
#include "Pimpact_TimeNSOp.hpp"
extern template class Pimpact::TimeNSBSmoother< Pimpact::TimeNSOp< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::TimeNSBSmoother< Pimpact::TimeNSOp< Pimpact::Space<double,int,4,4> > >;
#endif

#endif // end of #ifndef PIMPACT_TIMENSBSMOOTHER_HPP
