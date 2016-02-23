#pragma once
#ifndef PIMPACT_DIVGRADO2JSMOOTHER_HPP
#define PIMPACT_DIVGRADO2JSMOOTHER_HPP


#include "Pimpact_DivGradO2Op.hpp"




namespace Pimpact{


extern "C" {

void OP_DivGradO2JSmoother(
    const int& dimens,
    const int* const N,
    const int* const BL,
    const int* const BU,
    const int* const SR,
    const int* const ER,
    const double* const cdg1,
    const double* const cdg2,
    const double* const cdg3,
    const double& omega,
    const double* const b,
    const double* const x,
          double* const temp );

void OP_DivGradO2JBCSmoother(
    const int* const N,
    const int* const BL,
    const int* const BU,
    const int* const BCL,
    const int* const BCU,
    const int* const SR,
    const int* const ER,
    const double* const cdg1,
    const double* const cdg2,
    const double* const cdg3,
    const double& omBC,
    const double* const b,
    const double* const x,
          double* const temp );

}


/// \brief \f$\omega\f$-Jacobian smoother for second Order DivGradOp.
///
///
/// \relates DivGradO2Op
/// \ingroup BaseOperator
/// \todo instead of hardcode 2nd Order it would be pretty to use new space with StencilWidth<3,2>
/// \todo handle corner
template<class OperatorT>
class DivGradO2JSmoother {

public:

  using SpaceT = typename OperatorT::SpaceT;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

  Scalar omega_;
  int nIter_;
	bool levelYes_;

  Teuchos::RCP<DomainFieldT> temp_;

  const Teuchos::RCP<const OperatorT> op_;

public:

  DivGradO2JSmoother(
      const Teuchos::RCP<const SpaceT>& space):
    omega_( (2==space->dim())?0.8:6./7. ),
    nIter_( 2 ),
    levelYes_( false ),
    temp_( create<DomainFieldT>(space) ),
    op_( Teuchos::rcp( new OperatorT(space) ) ) {}

	/// \brief constructor
	///
	/// \param[in] op pointer to operator that is smoothed
	/// \param[in] pl  Parameter list of options for the multi grid solver.
	///   These are the options accepted by the solver manager:
	///   - "omega" - a \c Scalar damping factor. Default: for 2D 0.8 for 3D 6./7.  /
	///   - "numIters" - a \c int number of smoothing steps . Default: 4  /
	///   - "level" - a \c bool number of smoothing steps . Default: false  /
  DivGradO2JSmoother(
      const Teuchos::RCP<const OperatorT>& op,
      const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList() ):
    omega_( pl->get<Scalar>("omega", (2==op->space()->dim())?0.8:6./7. ) ),
    nIter_( pl->get<int>( "numIters", 2 ) ),
    levelYes_( pl->get<bool>( "level", false ) ),
    temp_( create<DomainFieldT>( op->space() ) ),
    op_(op) {}


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}( x - A y_k ) \f]
	void apply(const DomainFieldT& x, RangeFieldT& y,
			Belos::ETrans trans=Belos::NOTRANS ) const {

		for( int i=0; i<nIter_; ++i) {

      y.exchange();

			OP_DivGradO2JSmoother(
					space()->dim(),
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					op_->getSR(),
					op_->getER(),
					op_->getC(X),
					op_->getC(Y),
					op_->getC(Z),
					omega_,
					x.getConstRawPtr(),
					y.getConstRawPtr(),
					temp_->getRawPtr() );

			applyJBCSmoothing( x, y, *temp_ );

			temp_->changed();
			temp_->exchange();

			OP_DivGradO2JSmoother(
					space()->dim(),
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					op_->getSR(),
					op_->getER(),
					op_->getC(X),
					op_->getC(Y),
					op_->getC(Z),
					omega_,
					x.getConstRawPtr(),
					temp_->getConstRawPtr(),
					y.getRawPtr() );

			applyJBCSmoothing( x, *temp_, y );

			y.changed();
		}
		if( levelYes_ )
			y.level();

	}

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    out << "\t omega: " << omega_ << "\n";
    out << "\t numIter: " << nIter_ << "\n";
    op_->print( out );
  }

	const std::string getLabel() const { return( "DivGradO2JSmoother" ); };

protected:

	void applyJBCSmoothing( const DomainFieldT& b, const DomainFieldT& x, RangeFieldT& y ) const {
		// boundary conditions in X
		if( space()->getBCLocal()->getBCL(X)>0 ) {
			Ordinal i = 1;
			for( Ordinal k=op_->getSR(Z); k<=op_->getER(Z); ++k )
				for( Ordinal j=op_->getSR(Y); j<=op_->getER(Y); ++j )
					y.at(i,j,k) =
						(1-omega_)*x.at(i,j,k) + omega_/op_->getC(X,i,0)*( b.at(i,j,k) - op_->getC(X,i,+1)*x.at(i+1,j,k) );

		}
		if( space()->getBCLocal()->getBCU(X)>0 ) {
			Ordinal i = space()->nLoc(X);
			for( Ordinal k=op_->getSR(Z); k<=op_->getER(Z); ++k )
				for( Ordinal j=op_->getSR(Y); j<=op_->getER(Y); ++j )
					y.at(i,j,k) =
						(1-omega_)*x.at(i,j,k) + omega_/op_->getC(X,i,0) *( b.at(i,j,k) - op_->getC(X,i,-1)*x.at(i-1,j,k) );

		}

		// boundary conditions in Y
		if( space()->getBCLocal()->getBCL(Y)>0 ) {
			Ordinal j = 1;
			for( Ordinal k=op_->getSR(Z); k<=op_->getER(Z); ++k )
				for( Ordinal i=op_->getSR(X); i<=op_->getER(X); ++i )
					y.at(i,j,k) =
						(1-omega_)*x.at(i,j,k) + omega_/op_->getC(Y,j,0) *( b.at(i,j,k) - op_->getC(Y,j,+1)*x.at(i,j+1,k) );

		}
		if( space()->getBCLocal()->getBCU(Y)>0 ) {
			Ordinal j = space()->nLoc(Y);
			for( Ordinal k=op_->getSR(Z); k<=op_->getER(Z); ++k )
				for( Ordinal i=op_->getSR(X); i<=op_->getER(X); ++i )
					y.at(i,j,k) =
						(1-omega_)*x.at(i,j,k) + omega_/op_->getC(Y,j,0) *( b.at(i,j,k) - op_->getC(Y,j,-1)*x.at(i,j-1,k) );
		}

		// boundary conditions in Z
		if( space()->getBCLocal()->getBCL(Z)>0 ) {
			Ordinal k = 1;
			for( Ordinal j=op_->getSR(Y); j<=op_->getER(Y); ++j )
				for( Ordinal i=op_->getSR(X); i<=op_->getER(X); ++i )
					y.at(i,j,k) =
						(1-omega_)*x.at(i,j,k) + omega_/op_->getC(Z,k,0) *( b.at(i,j,k) - op_->getC(Z,k,+1)*x.at(i,j,k+1) );
		}
		if( space()->getBCLocal()->getBCU(Z)>0 ) {
			Ordinal k = space()->nLoc(Z);
			for( Ordinal j=op_->getSR(Y); j<=op_->getER(Y); ++j )
				for( Ordinal i=op_->getSR(X); i<=op_->getER(X); ++i )
					y.at(i,j,k) =
						(1-omega_)*x.at(i,j,k) + omega_/op_->getC(Z,k,0) *( b.at(i,j,k) - op_->getC(Z,k,-1)*x.at(i,j,k-1) );
		}
	}

}; // end of class DivGradO2JSmoother



/// \todo move somewhere better
template<template<class> class SmootherT, class OperatorT>
Teuchos::RCP< SmootherT<OperatorT> >
create(
    const Teuchos::RCP<OperatorT>& op,
    Teuchos::RCP<Teuchos::ParameterList> pl ) {

  return(
      Teuchos::rcp( new SmootherT<OperatorT>( op, pl ) ) );

}


/// \todo move somewhere better
template<class SmootherT, class OperatorT>
Teuchos::RCP< SmootherT >
create(
    const Teuchos::RCP< OperatorT>& op,
    Teuchos::RCP<Teuchos::ParameterList> pl ) {

  return(
      Teuchos::rcp( new SmootherT( op, pl ) ) );

}



} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,4> > >;
#endif


#endif // end of #ifndef PIMPACT_DIVGRADO2JSMOOTHER_HPP
