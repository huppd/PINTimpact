#pragma once
#ifndef PIMPACT_DIVGRADO2JSMOOTHER_HPP
#define PIMPACT_DIVGRADO2JSMOOTHER_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_DivGradO2Op.hpp"
#include "Pimpact_TeuchosTransfer.hpp"




namespace Pimpact{



/// \brief \f$\omega\f$-Jacobian smoother for second Order DivGradOp.
///
///
/// \relates DivGradO2Op
/// \ingroup BaseOperator
/// \note todo instead of hardcode 2nd Order it would be pretty to use new space with StencilWidth<3,2>
template<class OperatorT>
class DivGradO2JSmoother {

public:

  using SpaceT = typename OperatorT::SpaceT;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

	using SolverT = TeuchosSolver<OperatorT>;

  Scalar omega_;
  int nIter_;

	bool levelYes_;

  const Teuchos::RCP<const OperatorT> op_;

public:

	DivGradO2JSmoother( const Teuchos::RCP<const SpaceT>& space ):
    omega_( (2==SpaceT::sdim)?0.8:6./7. ),
    nIter_( 4 ),
    levelYes_( false ),
    op_( Teuchos::rcp( new OperatorT(space) ) ) {}

	/// \brief constructor
	///
	/// \param[in] op pointer to operator that is smoothed
	/// \param[in] pl  Parameter list of options for the multi grid solver.
	///   These are the options accepted by the solver manager:
	///   - "omega" - a \c Scalar damping factor. Default: for 2D 0.8 for 3D 6./7.  /
	///   - "numIters" - a \c int number of smoothing steps. Default: 4  /
	///   - "BC smoothing" - a \c int type of BC smoothing 0, 0: Jacbobian, else: direct. Default: 0 /
	///   - "debth" - for direct BC smoothing only meaning depth in wand normal direction of 2D BC problems. Default: 2 /
	///   - "level" - a \c bool number of smoothing steps. Default: false  /
  DivGradO2JSmoother(
      const Teuchos::RCP<const OperatorT>& op,
      const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList() ):
    omega_( pl->get<Scalar>("omega", (2==SpaceT::sdim)?0.8:6./7. ) ),
    nIter_( pl->get<int>( "numIters", 2 ) ),
    levelYes_( pl->get<bool>( "level", false ) ),
		op_(op) {}


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}( x - A y_k ) \f]
	void apply(const DomainFieldT& b, RangeFieldT& y, const Add& add=Add::N ) const {

		DomainFieldT temp( space() );

		for( int i=0; i<nIter_; ++i) {

      y.exchange();

			if( 3==SpaceT::sdim )
				for( Ordinal k=space()->begin(F::S,Z); k<=space()->end(F::S,Z); ++k )
					for( Ordinal j=space()->begin(F::S,Y); j<=space()->end(F::S,Y); ++j )
						for( Ordinal i=space()->begin(F::S,X); i<=space()->end(F::S,X); ++i ) {
							temp(i,j,k) = innerStenc3D( b, y, i,j,k);
						}
			else
				for( Ordinal k=space()->begin(F::S,Z); k<=space()->end(F::S,Z); ++k )
					for( Ordinal j=space()->begin(F::S,Y); j<=space()->end(F::S,Y); ++j )
						for( Ordinal i=space()->begin(F::S,X); i<=space()->end(F::S,X); ++i ) {
							temp(i,j,k) = innerStenc2D( b, y, i,j,k);
						}
			
			temp.changed();
			temp.exchange();

			if( 3==SpaceT::sdim )
				for( Ordinal k=space()->begin(F::S,Z); k<=space()->end(F::S,Z); ++k )
					for( Ordinal j=space()->begin(F::S,Y); j<=space()->end(F::S,Y); ++j )
						for( Ordinal i=space()->begin(F::S,X); i<=space()->end(F::S,X); ++i ) {
							y(i,j,k) = innerStenc3D( b, temp, i,j,k);
						}
			else
				for( Ordinal k=space()->begin(F::S,Z); k<=space()->end(F::S,Z); ++k )
					for( Ordinal j=space()->begin(F::S,Y); j<=space()->end(F::S,Y); ++j )
						for( Ordinal i=space()->begin(F::S,X); i<=space()->end(F::S,X); ++i ) {
							y(i,j,k) = innerStenc2D( b, temp, i,j,k);
						}

			y.changed();
		}
		if( levelYes_ )
			y.level();
	}

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    out << "\t omega: " << omega_ << "\n";
    out << "\t numIter: " << nIter_ << "\n";
    op_->print( out );
  }

	const std::string getLabel() const { return( "DivGradO2JSmoother" ); };

protected:

	constexpr Scalar innerStenc3D( const DomainFieldT& b, const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const { 

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->begin(F::S,X) ) ||
		           (space()->getBCLocal()->getBCU(X) > 0 && i==space()->end(F::S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->begin(F::S,Y) ) ||
		           (space()->getBCLocal()->getBCU(Y) > 0 && j==space()->end(F::S,Y) ) ;
		const bool bcZ = (space()->getBCLocal()->getBCL(Z) > 0 && k==space()->begin(F::S,Z) ) ||
		           (space()->getBCLocal()->getBCU(Z) > 0 && k==space()->end(F::S,Z) ) ;

		const Scalar& eps = 0.1;

		const Scalar epsX = (bcY||bcZ)?eps:1.;
		const Scalar epsY = (bcX||bcZ)?eps:1.;
		const Scalar epsZ = (bcX||bcY)?eps:1.;

		return(
				(1.-omega_)*x(i,j,k) +
				omega_/( epsX*getC(X,i,0) + epsY*getC(Y,j,0) + epsZ*getC(Z,k,0) )*(
					b(i,j,k) - 
				epsX*getC(X,i,-1)*x(i-1,j  ,k  ) - epsX*getC(X,i,1)*x(i+1,j  ,k  ) - 
				epsY*getC(Y,j,-1)*x(i  ,j-1,k  ) - epsY*getC(Y,j,1)*x(i  ,j+1,k  ) - 
				epsZ*getC(Z,k,-1)*x(i  ,j  ,k-1) - epsZ*getC(Z,k,1)*x(i  ,j  ,k+1) ) );
	} 

	constexpr Scalar innerStenc2D( const DomainFieldT& b, const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const { 

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->begin(F::S,X) ) ||
		           (space()->getBCLocal()->getBCU(X) > 0 && i==space()->end(F::S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->begin(F::S,Y) ) ||
		           (space()->getBCLocal()->getBCU(Y) > 0 && j==space()->end(F::S,Y) ) ;

		const Scalar& eps = 0.1;

		const Scalar epsX = bcY?eps:1.;
		const Scalar epsY = bcX?eps:1.;

		Scalar diag = (epsX*getC(X,i,0) + epsY*getC(Y,j,0));
		diag = (diag!=0.)?diag:1.;
		return(
				(1.-omega_)*x(i,j,k) +
				omega_/diag*( b(i,j,k) - 
				epsX*getC(X,i,-1)*x(i-1,j  ,k  ) - epsX*getC(X,i,1)*x(i+1,j  ,k  ) - 
				epsY*getC(Y,j,-1)*x(i  ,j-1,k  ) - epsY*getC(Y,j,1)*x(i  ,j+1,k  ) ) );
	} 

	constexpr const Scalar* getC( const ECoord& dir) const  { 
		return( op_->getC( dir ) ); 
	} 

	constexpr const Scalar* getC( const int& dir) const  { 
		return( op_->getC( dir ) ); 
	} 

	constexpr const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const  { 
		return( op_->getC( dir, i, off ) ); 
	} 

	constexpr const Scalar& getC( const int& dir, Ordinal i, Ordinal off ) const  { 
		return( op_->getC( dir, i, off ) ); 
	} 
	
}; // end of class DivGradO2JSmoother



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
