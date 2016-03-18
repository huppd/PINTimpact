#pragma once
#ifndef PIMPACT_DIVGRADO2JSMOOTHER_HPP
#define PIMPACT_DIVGRADO2JSMOOTHER_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_DivGradO2Op.hpp"
#include "Pimpact_TeuchosTransfer.hpp"




namespace Pimpact{


//extern "C"
//void OP_DivGradO2JSmoother(
    //const int& dimens,
    //const int* const N,
    //const int* const BL,
    //const int* const BU,
    //const int* const SR,
    //const int* const ER,
    //const double* const cdg1,
    //const double* const cdg2,
    //const double* const cdg3,
    //const double& omega,
    //const double* const b,
    //const double* const x,
          //double* const temp );




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

	using SolverT = TeuchosSolver<OperatorT>;

  Scalar omega_;
  int nIter_;
	int bcSmoothing_;
	Ordinal depth_;

	bool levelYes_;


  const Teuchos::RCP<const OperatorT> op_;

	Teuchos::Tuple< Teuchos::RCP<SolverT>, 3 > AsovL_;
	Teuchos::Tuple< Teuchos::RCP<SolverT>, 3 > AsovU_;

public:

	DivGradO2JSmoother( const Teuchos::RCP<const SpaceT>& space ):
    omega_( (2==space->dim())?0.8:6./7. ),
    nIter_( 2 ),
    levelYes_( false ),
		bcSmoothing_( 0 ),
		depth_( 2 ),
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
    omega_( pl->get<Scalar>("omega", (2==op->space()->dim())?0.8:6./7. ) ),
    nIter_( pl->get<int>( "numIters", 2 ) ),
		bcSmoothing_( pl->get<int>( "BC smoothing", 0 ) ),
		depth_( pl->get<Ordinal>( "depth", 2 ) ),
    levelYes_( pl->get<bool>( "level", false ) ),
    op_(op) {
		
			if( 0<bcSmoothing_ ) {
				Ordinal SS[3];
				Ordinal NN[3];

				// boundary conditions in X
				if( space()->getBCLocal()->getBCL(X)>0 ) {
					SS[X] = space()->sInd(EField::S,X);
					SS[Y] = space()->sInd(EField::S,Y);
					SS[Z] = space()->sInd(EField::S,Z);

					NN[X] = space()->sInd(EField::S,X)+depth_;
					NN[Y] = space()->eInd(EField::S,Y);
					NN[Z] = space()->eInd(EField::S,Z);

					// set solver and solve
					AsovL_[0] = Teuchos::rcp( new SolverT( op_, SS, NN) );
				}
				if( space()->getBCLocal()->getBCU(X)>0 ) {
					SS[X] = space()->eInd(EField::S,X)-depth_;
					SS[Y] = space()->sInd(EField::S,Y);
					SS[Z] = space()->sInd(EField::S,Z);

					NN[X] = space()->eInd(EField::S,X);
					NN[Y] = space()->eInd(EField::S,Y);
					NN[Z] = space()->eInd(EField::S,Z);

					// set solver and solve
					AsovU_[0] = Teuchos::rcp( new SolverT( op_, SS, NN ) );

				}

				// boundary conditions in Y
				if( space()->getBCLocal()->getBCL(Y)>0 ) {
					SS[X] = space()->sInd(EField::S,X);
					SS[Y] = space()->sInd(EField::S,Y);
					SS[Z] = space()->sInd(EField::S,Z);

					NN[X] = space()->eInd(EField::S,X);
					NN[Y] = space()->sInd(EField::S,Y)+depth_;
					NN[Z] = space()->eInd(EField::S,Z);

					// set solver and solve
					AsovL_[1] = Teuchos::rcp( new SolverT( op_, SS, NN ) );

				}
				if( space()->getBCLocal()->getBCU(Y)>0 ) {
					SS[X] = space()->sInd(EField::S,X);
					SS[Y] = space()->eInd(EField::S,Y)-depth_;
					SS[Z] = space()->sInd(EField::S,Z);

					NN[X] = space()->eInd(EField::S,X);
					NN[Y] = space()->eInd(EField::S,Y);
					NN[Z] = space()->eInd(EField::S,Z);

					// set solver and solve
					AsovU_[1] = Teuchos::rcp( new SolverT( op_, SS, NN ) );

				}

				// boundary conditions in Z
				if( space()->getBCLocal()->getBCL(Z)>0 ) {
					SS[X] = space()->sInd(EField::S,X);
					SS[Y] = space()->sInd(EField::S,Y);
					SS[Z] = space()->sInd(EField::S,Z);

					NN[X] = space()->eInd(EField::S,X);
					NN[Y] = space()->eInd(EField::S,Y);
					NN[Z] = space()->sInd(EField::S,Z)+depth_;

					// set solver and solve
					AsovL_[2] = Teuchos::rcp( new SolverT( op_, SS, NN ) );

				}
				if( space()->getBCLocal()->getBCU(Z)>0 ) {
					SS[X] = space()->sInd(EField::S,X);
					SS[Y] = space()->sInd(EField::S,Y);
					SS[Z] = space()->eInd(EField::S,Z)-depth_;

					NN[X] = space()->eInd(EField::S,X);
					NN[Y] = space()->eInd(EField::S,Y);
					NN[Z] = space()->eInd(EField::S,Z);

					// set solver and solve
					AsovU_[2] = Teuchos::rcp( new SolverT( op_, SS,NN ) );

				}

			}

		}


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}( x - A y_k ) \f]
	void apply(const DomainFieldT& b, RangeFieldT& y,
			Belos::ETrans trans=Belos::NOTRANS ) const {

		Teuchos::RCP<DomainFieldT> temp =
			Teuchos::rcp( new DomainFieldT(space()) );

		for( int i=0; i<nIter_; ++i) {

      y.exchange();

			if( 3==space()->dim() )
				for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
					for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
						for( Ordinal i=getSR(X); i<=getER(X); ++i ) {
							temp->at(i,j,k) = innerStenc3D( b, y, i,j,k);
						}
			else
				for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
					for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
						for( Ordinal i=getSR(X); i<=getER(X); ++i ) {
							temp->at(i,j,k) = innerStenc2D( b, y, i,j,k);
           }
			
			if( 0==bcSmoothing_ )
				applyJBCSmoothing( b, y, *temp, omega_ );
			else
				applyDBCSmoothing( b, *temp, omega_ );

			temp->changed();
			temp->exchange();

			if( 3==space()->dim() )
				for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
					for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
						for( Ordinal i=getSR(X); i<=getER(X); ++i ) {
							y.at(i,j,k) = innerStenc3D( b, *temp, i,j,k);
						}
			else
				for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
					for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
						for( Ordinal i=getSR(X); i<=getER(X); ++i ) {
							y.at(i,j,k) = innerStenc2D( b, *temp, i,j,k);
           }

			if( 0==bcSmoothing_ )
				applyJBCSmoothing( b, *temp, y, omega_ );
			else
				applyDBCSmoothing( b, y, omega_ );

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

	void applyDBCSmoothing( const DomainFieldT& b, DomainFieldT& x, const Scalar& omega ) const {

		// boundary conditions in 
		for( int i=0; i<3; ++i ) {
			if( space()->getBCLocal()->getBCL(i)>0 )
				AsovL_[i]->apply( b, x, omega );

			if( space()->getBCLocal()->getBCU(i)>0 )
				AsovU_[i]->apply( b, x, omega );
		}

	}


	void applyJBCSmoothing( const DomainFieldT& b, const DomainFieldT& x, RangeFieldT& y ) const {
		applyJBCSmoothing( b, x, y, omega_ );
	}


	void applyJBCSmoothing( const DomainFieldT& b, const DomainFieldT& x,
			RangeFieldT& y, const Scalar& omega ) const {

		for( int d=0; d<3; ++d ) {

			int d1 = getDir1( d );
			int d2 = getDir2( d );

			Teuchos::Tuple< Ordinal, 3> i;
			// boundary conditions in X
			if( space()->getBCLocal()->getBCL( d )>0 ) {
				i[d] = 1;
				for( i[d2]=getSR(d2); i[d2]<=getER(d2); ++i[d2] )
					for( i[d1]=getSR(d1); i[d1]<=getER(d1); ++i[d1] ) {

						Teuchos::Tuple< Ordinal, 3> ip = i;

						++ip[d];

						y.at(i[0],i[1],i[2]) =
							(1-omega)*x.at(i[0],i[1],i[2]) + omega/getC(d,i[d],0) *( b.at(i[0],i[1],i[2]) - getC(d,i[d],+1)*x.at(ip[0],ip[1],ip[2]) );
					}

			}
			if( space()->getBCLocal()->getBCU(d)>0 ) {
				i[d] = space()->nLoc(d);
				for( i[d2]=getSR(d2); i[d2]<=getER(d2); ++i[d2] )
					for( i[d1]=getSR(d1); i[d1]<=getER(d1); ++i[d1] ) {

						Teuchos::Tuple< Ordinal, 3> ip = i;

						--ip[d];

						y.at(i[0],i[1],i[2]) =
							(1-omega)*x.at(i[0],i[1],i[2]) + omega/getC(d,i[d],0) *( b.at(i[0],i[1],i[2]) - getC(d,i[d],-1)*x.at(ip[0],ip[1],ip[2]) );
					}

			}
		}

	}


	inline Scalar innerStenc3D( const DomainFieldT& b, const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const { 

		return(
				(1.-omega_)*x.at(i,j,k) +
				omega_/( getC(X,i,0) + getC(Y,j,0) + getC(Z,k,0) )*(
					b.at(i,j,k) - 
				getC(X,i,-1)*x.at(i-1,j  ,k  ) - getC(X,i,1)*x.at(i+1,j  ,k  ) - 
				getC(Y,j,-1)*x.at(i  ,j-1,k  ) - getC(Y,j,1)*x.at(i  ,j+1,k  ) - 
				getC(Z,k,-1)*x.at(i  ,j  ,k-1) - getC(Z,k,1)*x.at(i  ,j  ,k+1) ) );

	} 

	inline Scalar innerStenc2D( const DomainFieldT& b, const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const { 

		return(
				(1.-omega_)*x.at(i,j,k) +
				omega_/( getC(X,i,0) + getC(Y,j,0) )*(
					b.at(i,j,k) - 
				getC(X,i,-1)*x.at(i-1,j  ,k  ) - getC(X,i,1)*x.at(i+1,j  ,k  ) - 
				getC(Y,j,-1)*x.at(i  ,j-1,k  ) - getC(Y,j,1)*x.at(i  ,j+1,k  ) ) );

	} 

	inline const Scalar* getC( const ECoord& dir) const  { 
		return( op_->getC( dir ) ); 
	} 

	inline const Scalar* getC( const int& dir) const  { 
		return( op_->getC( dir ) ); 
	} 

	inline const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const  { 
		return( op_->getC( dir, i, off ) ); 
	} 

	inline const Scalar& getC( const int& dir, Ordinal i, Ordinal off ) const  { 
		return( op_->getC( dir, i, off ) ); 
	} 

	inline const Ordinal* getSR() const { return( op_->getSR() ); } 
	inline const Ordinal* getER() const { return( op_->getER() ); } 

	inline const Ordinal& getSR( const Ordinal& coord ) const { return( op_->getSR( coord ) ); } 
	inline const Ordinal& getER( const Ordinal& coord ) const { return( op_->getER( coord ) ); } 

	inline const Ordinal& getSR( const ECoord& coord ) const { return( op_->getSR( coord ) ); } 
	inline const Ordinal& getER( const ECoord& coord ) const { return( op_->getER( coord ) ); } 

	
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
