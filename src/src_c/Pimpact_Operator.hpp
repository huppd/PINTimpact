#pragma once
#ifndef PIMPACT_OPERATOR_HPP
#define PIMPACT_OPERATOR_HPP

#include"Teuchos_RCP.hpp"

#include"Pimpact_ScalarField.hpp"
#include"Pimpact_VectorField.hpp"
#include"Pimpact_CompoundField.hpp"
#include"Pimpact_OperatorMV.hpp"
#include"Pimpact_LinearProblem.hpp"

namespace Pimpact{


extern "C" {
	void OP_grad( const int& m, double* phi, double *grad );
  void OP_div( double* phiU, double* phiV, double* phiW, double* div);
  void OP_helmholtz( const int& m, const bool& exch_yes, const double& mulI, const double& multL, double* phi, double* Lap);
  void OP_div_grad( const bool& corner_yes, double* phi, double* lap );
  void OP_bc_extrapolation( const int& m, double* phi );
}


template<class Scalar,class Ordinal>
class Grad {
public:
	typedef ScalarField<Scalar,Ordinal>  DomainFieldT;
	typedef VectorField<Scalar,Ordinal>  RangeFieldT;
	typedef NonModeOp OpType;

	void apply(const DomainFieldT& x, RangeFieldT& y) const {
		int dim = x.getFieldSpace()->dim_;
		for( int i=0; i<dim; ++i) {
			OP_grad(i+1,x.s_,y.vec_[i]);
//				OP_bc_extrapolation( i+1, y.vec_[i] );
		}
	}
	bool hasApplyTranspose() const { return false; }
};


template<class Scalar,class Ordinal>
class Div {
public:
	typedef VectorField<Scalar,Ordinal>  DomainFieldT;
	typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
	typedef NonModeOp OpType;

	void apply(const DomainFieldT& x, RangeFieldT& y) const {
		OP_div(x.vec_[0],x.vec_[1],x.vec_[2],y.s_);
	} // return erro
	bool hasApplyTranspose() const { return false; }
};


/**
 * \brief Helmholtz operator
 * \todo change 2 to dim
 */
template<class Scalar,class Ordinal>
class Helmholtz {
protected:
	Scalar mulI_;
	Scalar mulL_;
public:
	Helmholtz():mulI_(1.),mulL_(1.) {};
	Helmholtz(Scalar mulI, Scalar mulL):mulI_(mulI),mulL_(mulL) {};

	void setMulI(Scalar mulI){ mulI_ = mulI;};
	void setMulL(Scalar mulL){ mulL_ = mulL;};

	typedef VectorField<Scalar,Ordinal>  DomainFieldT;
	typedef VectorField<Scalar,Ordinal>  RangeFieldT;
	typedef NonModeOp OpType;

	void apply(const DomainFieldT& x, RangeFieldT& y) const {
//  		y.assign(x);
//  		int d=0;
		for( int d=0; d<2; ++d )
			OP_helmholtz( d+1, true, mulI_, mulL_, x.vec_[d], y.vec_[d] ) ;
	}
	bool hasApplyTranspose() const { return false; }
};


template<class Scalar,class Ordinal>
Teuchos::RCP< OperatorMV<Helmholtz<Scalar,Ordinal> > > createHelmholtz( Scalar mulI, Scalar mulL ) {
	return Teuchos::rcp( new OperatorMV<Helmholtz<Scalar,Ordinal> >( Teuchos::rcp( new Helmholtz<Scalar,Ordinal>(mulI, mulL) ) ) );
}


/**
 * \brief "laplace" for pressure
 * \todo not workin properly?
 * \todo add temporary variable
 */
template<class Scalar,class Ordinal>
class Div_Grad {
public:
	typedef ScalarField<Scalar,Ordinal>  DomainFieldT;
	typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
	typedef NonModeOp OpType;

	void apply(const DomainFieldT& x, RangeFieldT& y) const {
		y.init();
		OP_div_grad( true, x.s_, y.s_ );
	}

	bool hasApplyTranspose() const { return false; }
};


template<class Scalar,class Ordinal>
class nonlinear_Udx {
public:
	typedef VectorField<Scalar,Ordinal>  DomainFieldT;
	typedef VectorField<Scalar,Ordinal>  RangeFieldT;
	typedef NonModeOp OpType;
private:
	Teuchos::RCP<DomainFieldT> u_;
public:

	void set_U(Teuchos::RCP<DomainFieldT> u) { u_=u; }

	void apply(const DomainFieldT& x, RangeFieldT& y) const { return; }
	bool hasApplyTranspose() const { return true; }
};


template<class Scalar,class Ordinal>
class nonlinear_dxU {
public:
	typedef VectorField<Scalar,Ordinal>  DomainFieldT;
	typedef VectorField<Scalar,Ordinal>  RangeFieldT;
	typedef NonModeOp OpType;
private:
	Teuchos::RCP<DomainFieldT> u_;
public:

	void set_U(Teuchos::RCP<DomainFieldT> u) { u_=u; }

	void apply(const DomainFieldT& x, RangeFieldT& y) const { return; }
	bool hasApplyTranspose() const { return true; }
};


template< class Scalar, class Ordinal >
class Div_Hinv_Grad {
public:
	typedef ScalarField<Scalar,Ordinal>  DomainFieldT;
	typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
	typedef MultiField<VectorField<Scalar,Ordinal> > MVF;
	typedef OperatorMV<Helmholtz<Scalar,Ordinal> > HType;
	typedef NonModeOp OpType;
private:
	Teuchos::RCP< MVF > temp0_;
	Teuchos::RCP< MVF > temp1_;
	Teuchos::RCP< Div<Scalar,Ordinal> > div_;
	Teuchos::RCP< Grad<Scalar,Ordinal> > grad_;
	Teuchos::RCP< LinearProblem<Scalar, MVF, HType > > H_;

public:
	Div_Hinv_Grad( Teuchos::RCP<MVF> temp,
//  			Teuchos::RCP<Div<Scalar,Ordinal> > div,
//  			Teuchos::RCP<Grad<Scalar,Ordinal> > grad,
			Teuchos::RCP< LinearProblem<Scalar, MVF, HType > > H ):
				temp0_(temp->Clone(1)), temp1_(temp->Clone(1)),
				div_(Teuchos::rcp( new Div<Scalar,Ordinal> ) ),
				grad_(Teuchos::rcp( new Grad<Scalar,Ordinal> ) ),
				H_(H) {};

	void apply(const DomainFieldT& x, RangeFieldT& y) const {
		grad_->apply(x,temp0_->GetVec(0) );
		H_->solve( temp1_, temp0_);
		div_->apply(temp1_->GetVec(0),y);
		return;
	}
	bool hasApplyTranspose() const { return false; }

};


template< class Scalar, class Ordinal>
Teuchos::RCP<OperatorMV< Div_Hinv_Grad<Scalar,Ordinal> > > createDivHinvGrad(
		Teuchos::RCP<MultiField<VectorField<Scalar,Ordinal> > > temp,
		Teuchos::RCP< LinearProblem<Scalar, MultiField<VectorField<Scalar,Ordinal> >, OperatorMV<Helmholtz<Scalar,Ordinal> > > > H ) {

	return Teuchos::rcp( new OperatorMV<Div_Hinv_Grad<Scalar,Ordinal> >( Teuchos::rcp( new Div_Hinv_Grad<Scalar,Ordinal>( temp, H ) ) ) );

}


template<class Scalar,class Ordinal>
class Dt {
	Scalar omega_;
public:
	Dt():omega_(1.) {};
	Dt(Scalar omega):omega_(omega) {};

	typedef VectorField<Scalar,Ordinal>  DomainFieldT;
	typedef VectorField<Scalar,Ordinal>  RangeFieldT;
	typedef ModeOp OpType;


	void apply(const ModeField<DomainFieldT>& x, ModeField<RangeFieldT>& y ) const {
		y.getFieldC()->add( 0., *x.getConstFieldC(), omega_, *x.getConstFieldS() );
		y.getFieldS()->add( -omega_, *x.getConstFieldC(), 0., *x.getConstFieldS() );
//		y.getFieldC()->assign( x.getConstFieldS() );
//		y.getFieldC()->scale( -omega_ );
//		y.getFieldS()->assign( x.getConstFieldC() );
//		y.getFieldS()->scale( omega_ );
	}
	bool hasApplyTranspose() const { return false; }
};


template< class Scalar, class Ordinal>
Teuchos::RCP<OperatorMV< Dt<Scalar,Ordinal> > > createDt( Scalar omega = 1. ) {
	return Teuchos::rcp( new OperatorMV<Dt<Scalar,Ordinal> >( Teuchos::rcp( new Dt<Scalar,Ordinal>( omega ) ) ) );
}


template<class Scalar,class Ordinal>
class DtL {
	Scalar omega_;
	Teuchos::RCP<Helmholtz<Scalar,Ordinal> > L_;
public:

	DtL():omega_(1.),L_(Teuchos::rcp(new Helmholtz<Scalar,Ordinal>( 0., 1. )) ) {};

	DtL( Scalar omega, Scalar mulI, Scalar mulL):omega_(omega),L_(Teuchos::rcp( new Helmholtz<Scalar,Ordinal>(mulI,mulL) ) ) {};

	DtL( Scalar omega, Teuchos::RCP<Helmholtz<Scalar,Ordinal> > L ):omega_(omega),L_(L) {};

	typedef VectorField<Scalar,Ordinal>  DomainFieldT;
	typedef VectorField<Scalar,Ordinal>  RangeFieldT;
	typedef ModeOp OpType;


	void apply(const ModeField<DomainFieldT>& x, ModeField<RangeFieldT>& y ) const {
		L_->apply( *x.getConstFieldC(), *y.getFieldC() );
		y.getFieldC()->add( 1., *y.getConstFieldC(), omega_, *x.getConstFieldS() );
		L_->apply( *x.getConstFieldS(), *y.getFieldS() );
		y.getFieldS()->add( -omega_, *x.getConstFieldC(), 1., *y.getConstFieldS() );
	}

	bool hasApplyTranspose() const { return false; }
};

template< class Scalar, class Ordinal>
Teuchos::RCP<OperatorMV< DtL<Scalar,Ordinal> > > createDtL( Scalar omega=1., Scalar mulI=0., Scalar mulL=1. ) {
	return Teuchos::rcp( new OperatorMV<DtL<Scalar,Ordinal> >( Teuchos::rcp( new DtL<Scalar,Ordinal>( omega, mulI, mulL ) ) ) );
}


template< class Scalar, class Ordinal >
class Div_DtLinv_Grad {
public:
	typedef ScalarField<Scalar,Ordinal>  DomainFieldT;
	typedef ScalarField<Scalar,Ordinal>  RangeFieldT;
	typedef MultiField<ModeField<VectorField<Scalar,Ordinal> > > MVF;
	typedef OperatorMV<DtL<Scalar,Ordinal> > HType;
	typedef ModeOp OpType;
private:
	Teuchos::RCP< MVF > temp0_;
	Teuchos::RCP< MVF > temp1_;
	Teuchos::RCP< Div<Scalar,Ordinal> > div_;
	Teuchos::RCP< Grad<Scalar,Ordinal> > grad_;
	Teuchos::RCP< LinearProblem<Scalar, MVF, HType > > H_;

public:
	Div_DtLinv_Grad( Teuchos::RCP<MVF> temp,
//  			Teuchos::RCP<Div<Scalar,Ordinal> > div,
//  			Teuchos::RCP<Grad<Scalar,Ordinal> > grad,
			Teuchos::RCP< LinearProblem<Scalar, MVF, HType > > H ):
				temp0_(temp->Clone(1)), temp1_(temp->Clone(1)),
				div_(Teuchos::rcp( new Div<Scalar,Ordinal> ) ),
				grad_(Teuchos::rcp( new Grad<Scalar,Ordinal> ) ),
				H_(H) {};

	void apply(const ModeField<DomainFieldT>& x, ModeField<RangeFieldT>& y) const {
		grad_->apply( *x.getConstFieldC(), *temp0_->GetVec(0).getFieldC() );
		grad_->apply( *x.getConstFieldS(), *temp0_->GetVec(0).getFieldS() );
		H_->solve( temp1_, temp0_);
		div_->apply( *temp1_->GetVec(0).getConstFieldC(), *y.getFieldC() );
		div_->apply( *temp1_->GetVec(0).getConstFieldS(), *y.getFieldS() );
	}

	bool hasApplyTranspose() const { return false; }

};


template< class Scalar, class Ordinal>
Teuchos::RCP<OperatorMV< Div_DtLinv_Grad<Scalar,Ordinal> > > createDivDtLinvGrad(
		Teuchos::RCP<MultiField<ModeField<VectorField<Scalar,Ordinal> > > > temp,
		Teuchos::RCP< LinearProblem<Scalar, MultiField<ModeField<VectorField<Scalar,Ordinal> > >, OperatorMV<DtL<Scalar,Ordinal> > > > H ) {

	return Teuchos::rcp( new OperatorMV<Div_DtLinv_Grad<Scalar,Ordinal> >( Teuchos::rcp( new Div_DtLinv_Grad<Scalar,Ordinal>( temp, H) ) ) );
}

template<class Scalar,class Ordinal>
class CompoundStokes {
	Scalar omega_;
	Teuchos::RCP<Helmholtz<Scalar,Ordinal> > L_;
	Teuchos::RCP<Div<Scalar,Ordinal> > div_;
	Teuchos::RCP<Grad<Scalar,Ordinal> > grad_;

	typedef ScalarField<Scalar,Ordinal>  SF;
	typedef VectorField<Scalar,Ordinal>  VF;
	typedef ModeField<SF>  	             MSF;
	typedef ModeField<VF>                MVF;
public:

//	CompoundStokes():omega_(1.),L_(Teuchos::rcp(new Helmholtz<Scalar,Ordinal>( 0., 1. )) ) {};

	CompoundStokes( Scalar omega, Scalar mulI, Scalar mulL, Teuchos::RCP<VF> temp ):
		omega_(omega),L_(Teuchos::rcp( new Helmholtz<Scalar,Ordinal>(mulI,mulL) ) ),
		div_(Teuchos::rcp(new Div<Scalar,Ordinal>())), grad_(Teuchos::rcp(new Grad<Scalar,Ordinal>() )), temp_(temp) {};

//	CompoundStokes( Scalar omega, Teuchos::RCP<Helmholtz<Scalar,Ordinal> > L ):omega_(omega),L_(L) {};

	typedef CompoundField<MVF,MSF>  DomainFieldT;
	typedef CompoundField<MVF,MSF>  RangeFieldT;
	typedef ModeOp OpType;


	void apply(const DomainFieldT& x_, RangeFieldT& y_ ) const {
		auto x = x_.getConstVField();
		auto y = y_.getVField();
		auto xp = x_.getConstSField();
		auto yp = y_.getSField();
		// H-blockz
		L_->apply( *x->getConstFieldC(), *y->getFieldC() );
		y->getFieldC()->add( 1., *y->getConstFieldC(), omega_, *x->getConstFieldS() );
		L_->apply( *x->getConstFieldS(), *y->getFieldS() );
		y->getFieldS()->add( -omega_, *x->getConstFieldC(), 1., *y->getConstFieldS() );
		// div
		div_->apply( *x->getConstFieldC(), *yp->getFieldC() );
		div_->apply( *x->getConstFieldS(), *yp->getFieldS() );
		// grad pressure
		grad_->apply( *xp->getConstFieldC(), *temp_ );
		y->getFieldC()->add( 1., *y->getFieldC(), 1., *temp_ );
		grad_->apply( *xp->getConstFieldS(), *temp_ );
		y->getFieldS()->add( 1., *y->getFieldS(), 1., *temp_ );
	}

	bool hasApplyTranspose() const { return false; }
protected:
	Teuchos::RCP<VF> temp_;
};

template< class Scalar, class Ordinal>
Teuchos::RCP<OperatorMV< CompoundStokes<Scalar,Ordinal> > > createCompoundStokes(
		Scalar omega, Scalar mulI, Scalar mulL, Teuchos::RCP<VectorField<Scalar,Ordinal> > temp ) {

	return Teuchos::rcp( new OperatorMV<CompoundStokes<Scalar,Ordinal> >( Teuchos::rcp( new CompoundStokes<Scalar,Ordinal>( omega,mulI,mulL,temp) ) ) );
}
} // end of namespace Pimpact

#endif
