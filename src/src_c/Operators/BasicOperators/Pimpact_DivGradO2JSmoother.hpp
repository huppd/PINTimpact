#pragma once
#ifndef PIMPACT_DIVGRADO2JSMOOTHER_HPP
#define PIMPACT_DIVGRADO2JSMOOTHER_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialQRDenseSolver.hpp"

#include "Pimpact_DivGradO2Op.hpp"
#include "Pimpact_TeuchosTransfer.hpp"




namespace Pimpact{


extern "C"
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

	using VectorT = Teuchos::SerialDenseVector<Ordinal,Scalar>;
	using MatrixT = Teuchos::SerialDenseMatrix<Ordinal,Scalar>;

	//using SolverT = Teuchos::SerialQRDenseSolver<Ordinal,Scalar>;
	using SolverT = Teuchos::SerialDenseSolver<Ordinal,Scalar>;

  Scalar omega_;
  int nIter_;
	int bcSmoothing_;
	Ordinal depth_;

	bool levelYes_;

  Teuchos::RCP<DomainFieldT> temp_;

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
    temp_( Teuchos::rcp( new DomainFieldT(space) ) ),
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
    temp_( create<DomainFieldT>( op->space() ) ),
    op_(op) {
		
			if( 0<bcSmoothing_ ) {
				Ordinal SS[3];
				Ordinal NN[3];

				Teuchos::RCP<TeuchosTransfer<SpaceT> > trans;

				// boundary conditions in X
				if( space()->getBCLocal()->getBCL(X)>0 ) {
					SS[X] = space()->sInd(EField::S,X);
					SS[Y] = space()->sInd(EField::S,Y);
					SS[Z] = space()->sInd(EField::S,Z);

					NN[X] = space()->sInd(EField::S,X)+depth_;
					NN[Y] = space()->eInd(EField::S,Y);
					NN[Z] = space()->eInd(EField::S,Z);

					trans = Teuchos::rcp( new TeuchosTransfer<SpaceT>( space(), SS, NN ) );

					Teuchos::RCP<MatrixT> A = Teuchos::rcp( new MatrixT( trans->getN(), trans->getN(), true ) );
					trans->apply( op_, A );

					// set solver and solve
					AsovL_[0] = Teuchos::rcp( new SolverT() );
					AsovL_[0]->factorWithEquilibration( true );
					AsovL_[0]->setMatrix( A );
					AsovL_[0]->factor();

				}
				if( space()->getBCLocal()->getBCU(X)>0 ) {
					SS[X] = space()->eInd(EField::S,X)-depth_;
					SS[Y] = space()->sInd(EField::S,Y);
					SS[Z] = space()->sInd(EField::S,Z);

					NN[X] = space()->eInd(EField::S,X);
					NN[Y] = space()->eInd(EField::S,Y);
					NN[Z] = space()->eInd(EField::S,Z);

					trans = Teuchos::rcp( new TeuchosTransfer<SpaceT>( space(), SS, NN ) );

					Teuchos::RCP<MatrixT> A = Teuchos::rcp( new MatrixT( trans->getN(), trans->getN(), true ) );
					trans->apply( op_, A );

					// set solver and solve
					AsovU_[0] = Teuchos::rcp( new SolverT() );
					AsovU_[0]->factorWithEquilibration( true );
					AsovU_[0]->setMatrix( A );
					AsovU_[0]->factor();

				}

				// boundary conditions in Y
				if( space()->getBCLocal()->getBCL(Y)>0 ) {
					SS[X] = space()->sInd(EField::S,X);
					SS[Y] = space()->sInd(EField::S,Y);
					SS[Z] = space()->sInd(EField::S,Z);

					NN[X] = space()->eInd(EField::S,X);
					NN[Y] = space()->sInd(EField::S,Y)+depth_;
					NN[Z] = space()->eInd(EField::S,Z);

					trans = Teuchos::rcp( new TeuchosTransfer<SpaceT>( space(), SS, NN ) );

					Teuchos::RCP<MatrixT> A = Teuchos::rcp( new MatrixT( trans->getN(), trans->getN(), true ) );
					trans->apply( op_, A );

					// set solver and solve
					AsovL_[1] = Teuchos::rcp( new SolverT() );
					AsovL_[1]->factorWithEquilibration( true );
					AsovL_[1]->setMatrix( A );
					AsovL_[1]->factor();

				}
				if( space()->getBCLocal()->getBCU(Y)>0 ) {
					SS[X] = space()->sInd(EField::S,X);
					SS[Y] = space()->eInd(EField::S,Y)-depth_;
					SS[Z] = space()->sInd(EField::S,Z);

					NN[X] = space()->eInd(EField::S,X);
					NN[Y] = space()->eInd(EField::S,Y);
					NN[Z] = space()->eInd(EField::S,Z);

					trans = Teuchos::rcp( new TeuchosTransfer<SpaceT>( space(), SS, NN ) );

					Teuchos::RCP<MatrixT> A = Teuchos::rcp( new MatrixT( trans->getN(), trans->getN(), true ) );
					trans->apply( op_, A );

					// set solver and solve
					AsovU_[1] = Teuchos::rcp( new SolverT() );
					AsovU_[1]->factorWithEquilibration( true );
					AsovU_[1]->setMatrix( A );
					AsovU_[1]->factor();

				}

				// boundary conditions in Z
				if( space()->getBCLocal()->getBCL(Z)>0 ) {
					SS[X] = space()->sInd(EField::S,X);
					SS[Y] = space()->sInd(EField::S,Y);
					SS[Z] = space()->sInd(EField::S,Z);

					NN[X] = space()->eInd(EField::S,X);
					NN[Y] = space()->eInd(EField::S,Y);
					NN[Z] = space()->sInd(EField::S,Z)+depth_;

					trans = Teuchos::rcp( new TeuchosTransfer<SpaceT>( space(), SS, NN ) );

					Teuchos::RCP<MatrixT> A = Teuchos::rcp( new MatrixT( trans->getN(), trans->getN(), true ) );
					trans->apply( op_, A );

					// set solver and solve
					AsovL_[2] = Teuchos::rcp( new SolverT() );
					AsovL_[2]->factorWithEquilibration( true );
					AsovL_[2]->setMatrix( A );
					AsovL_[2]->factor();

				}
				if( space()->getBCLocal()->getBCU(Z)>0 ) {
					SS[X] = space()->sInd(EField::S,X);
					SS[Y] = space()->sInd(EField::S,Y);
					SS[Z] = space()->eInd(EField::S,Z)-depth_;

					NN[X] = space()->eInd(EField::S,X);
					NN[Y] = space()->eInd(EField::S,Y);
					NN[Z] = space()->eInd(EField::S,Z);

					trans = Teuchos::rcp( new TeuchosTransfer<SpaceT>( space(), SS, NN ) );

					Teuchos::RCP<MatrixT> A = Teuchos::rcp( new MatrixT( trans->getN(), trans->getN(), true ) );
					trans->apply( op_, A );

					// set solver and solve
					AsovU_[2] = Teuchos::rcp( new SolverT() );
					AsovU_[2]->factorWithEquilibration( true );
					AsovU_[2]->setMatrix( A );
					AsovU_[2]->factor();

				}

			}

		}


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}( x - A y_k ) \f]
	void apply(const DomainFieldT& b, RangeFieldT& y,
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
					b.getConstRawPtr(),
					y.getConstRawPtr(),
					temp_->getRawPtr() );

			if( 0==bcSmoothing_ )
				applyJBCSmoothing( b, y, *temp_, omega_ );
			else
				applyDBCSmoothing( b, *temp_, omega_ );

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
					b.getConstRawPtr(),
					temp_->getConstRawPtr(),
					y.getRawPtr() );

			if( 0==bcSmoothing_ )
				applyJBCSmoothing( b, *temp_, y, omega_ );
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


		Ordinal SS[3];
		Ordinal NN[3];

		Teuchos::RCP<TeuchosTransfer<SpaceT> > trans;

		// boundary conditions in X
		if( space()->getBCLocal()->getBCL(X)>0 ) {
			SS[X] = space()->sInd(EField::S,X);
			SS[Y] = space()->sInd(EField::S,Y);
			SS[Z] = space()->sInd(EField::S,Z);

			NN[X] = space()->sInd(EField::S,X)+depth_;
			NN[Y] = space()->eInd(EField::S,Y);
			NN[Z] = space()->eInd(EField::S,Z);

			trans = Teuchos::rcp( new TeuchosTransfer<SpaceT>( space(), SS, NN ) );

			Teuchos::RCP<MatrixT> A = Teuchos::rcp( new MatrixT( trans->getN(), trans->getN(), true ) );

			Teuchos::RCP<VectorT> X = Teuchos::rcp( new VectorT(trans->getN(), false) );
			Teuchos::RCP<VectorT>	B = Teuchos::rcp( new VectorT(trans->getN(), false) );

			x.exchange(); // ???
			trans->apply( op_, A );
			trans->apply( b, B );
			trans->updateRHS( op_, x, B );

			// set solver and solve
			AsovL_[0]->setVectors( X, B );
			AsovL_[0]->solve();

			trans->apply( X, x, omega);
			//y.setCornersZero(); // ???
			x.changed();

		}
		if( space()->getBCLocal()->getBCU(X)>0 ) {
			SS[X] = space()->eInd(EField::S,X)-depth_;
			SS[Y] = space()->sInd(EField::S,Y);
			SS[Z] = space()->sInd(EField::S,Z);

			NN[X] = space()->eInd(EField::S,X);
			NN[Y] = space()->eInd(EField::S,Y);
			NN[Z] = space()->eInd(EField::S,Z);

			trans = Teuchos::rcp( new TeuchosTransfer<SpaceT>( space(), SS, NN ) );

			Teuchos::RCP<MatrixT> A = Teuchos::rcp( new MatrixT( trans->getN(), trans->getN(), true ) );

			Teuchos::RCP<VectorT> X = Teuchos::rcp( new VectorT(trans->getN(), false) );
			Teuchos::RCP<VectorT>	B = Teuchos::rcp( new VectorT(trans->getN(), false) );

			x.exchange(); // ???
			trans->apply( b, B );
			trans->updateRHS( op_, x, B );

			// set solver and solve
			AsovU_[0]->setVectors( X, B );
			AsovU_[0]->solve();

			trans->apply( X, x, omega);
			//y.setCornersZero(); // ???
			x.changed();
		}

		// boundary conditions in Y
		if( space()->getBCLocal()->getBCL(Y)>0 ) {
			SS[X] = space()->sInd(EField::S,X);
			SS[Y] = space()->sInd(EField::S,Y);
			SS[Z] = space()->sInd(EField::S,Z);

			NN[X] = space()->eInd(EField::S,X);
			NN[Y] = space()->sInd(EField::S,Y)+depth_;
			NN[Z] = space()->eInd(EField::S,Z);

			trans = Teuchos::rcp( new TeuchosTransfer<SpaceT>( space(), SS, NN ) );

			Teuchos::RCP<MatrixT> A = Teuchos::rcp( new MatrixT( trans->getN(), trans->getN(), true ) );

			Teuchos::RCP<VectorT> X = Teuchos::rcp( new VectorT(trans->getN(), false) );
			Teuchos::RCP<VectorT>	B = Teuchos::rcp( new VectorT(trans->getN(), false) );

			x.exchange(); // ???
			trans->apply( b, B );
			trans->updateRHS( op_, x, B );

			// set solver and solve
			AsovL_[1]->setVectors( X, B );
			AsovL_[1]->solve();

			trans->apply( X, x, omega);
			//y.setCornersZero(); // ???
			x.changed();
		}
		if( space()->getBCLocal()->getBCU(Y)>0 ) {
			SS[X] = space()->sInd(EField::S,X);
			SS[Y] = space()->eInd(EField::S,Y)-depth_;
			SS[Z] = space()->sInd(EField::S,Z);

			NN[X] = space()->eInd(EField::S,X);
			NN[Y] = space()->eInd(EField::S,Y);
			NN[Z] = space()->eInd(EField::S,Z);

			trans = Teuchos::rcp( new TeuchosTransfer<SpaceT>( space(), SS, NN ) );

			Teuchos::RCP<VectorT> X = Teuchos::rcp( new VectorT(trans->getN(), false) );
			Teuchos::RCP<VectorT>	B = Teuchos::rcp( new VectorT(trans->getN(), false) );

			x.exchange(); // ???
			trans->apply( b, B );
			trans->updateRHS( op_, x, B );

			// set solver and solve
			AsovU_[1]->setVectors( X, B );
			AsovU_[1]->solve();

			trans->apply( X, x, omega);
			//y.setCornersZero(); // ???
			x.changed();
		}

		// boundary conditions in Z
		if( space()->getBCLocal()->getBCL(Z)>0 ) {
			SS[X] = space()->sInd(EField::S,X);
			SS[Y] = space()->sInd(EField::S,Y);
			SS[Z] = space()->sInd(EField::S,Z);

			NN[X] = space()->eInd(EField::S,X);
			NN[Y] = space()->eInd(EField::S,Y);
			NN[Z] = space()->sInd(EField::S,Z)+depth_;

			trans = Teuchos::rcp( new TeuchosTransfer<SpaceT>( space(), SS, NN ) );

			Teuchos::RCP<VectorT> X = Teuchos::rcp( new VectorT(trans->getN(), false) );
			Teuchos::RCP<VectorT>	B = Teuchos::rcp( new VectorT(trans->getN(), false) );

			x.exchange(); // ???
			trans->apply( b, B );
			trans->updateRHS( op_, x, B );

			// set solver and solve
			AsovL_[2]->setVectors( X, B );
			AsovL_[2]->solve();

			trans->apply( X, x, omega);
			//y.setCornersZero(); // ???
			x.changed();
		}
		if( space()->getBCLocal()->getBCU(Z)>0 ) {
			SS[X] = space()->sInd(EField::S,X);
			SS[Y] = space()->sInd(EField::S,Y);
			SS[Z] = space()->eInd(EField::S,Z)-depth_;

			NN[X] = space()->eInd(EField::S,X);
			NN[Y] = space()->eInd(EField::S,Y);
			NN[Z] = space()->eInd(EField::S,Z);

			trans = Teuchos::rcp( new TeuchosTransfer<SpaceT>( space(), SS, NN ) );

			Teuchos::RCP<VectorT> X = Teuchos::rcp( new VectorT(trans->getN(), false) );
			Teuchos::RCP<VectorT>	B = Teuchos::rcp( new VectorT(trans->getN(), false) );

			x.exchange(); // ???
			trans->apply( b, B );
			trans->updateRHS( op_, x, B );

			// set solver and solve
			AsovU_[2]->setVectors( X, B );
			AsovU_[2]->solve();

			trans->apply( X, x, omega );
			//y.setCornersZero(); // ???
			x.changed();
		}

	}

	void applyJBCSmoothing( const DomainFieldT& b, const DomainFieldT& x, RangeFieldT& y ) const {
		applyJBCSmoothing( b, x, y, omega_ );
	}

	void applyJBCSmoothing( const DomainFieldT& b, const DomainFieldT& x, RangeFieldT& y, const Scalar& omega ) const {
		// boundary conditions in X
		if( space()->getBCLocal()->getBCL(X)>0 ) {
			Ordinal i = 1;
			for( Ordinal k=op_->getSR(Z); k<=op_->getER(Z); ++k )
				for( Ordinal j=op_->getSR(Y); j<=op_->getER(Y); ++j )
					y.at(i,j,k) =
						(1-omega)*x.at(i,j,k) + omega/op_->getC(X,i,0)*( b.at(i,j,k) - op_->getC(X,i,+1)*x.at(i+1,j,k) );

		}
		if( space()->getBCLocal()->getBCU(X)>0 ) {
			Ordinal i = space()->nLoc(X);
			for( Ordinal k=op_->getSR(Z); k<=op_->getER(Z); ++k )
				for( Ordinal j=op_->getSR(Y); j<=op_->getER(Y); ++j )
					y.at(i,j,k) =
						(1-omega)*x.at(i,j,k) + omega/op_->getC(X,i,0) *( b.at(i,j,k) - op_->getC(X,i,-1)*x.at(i-1,j,k) );

		}

		// boundary conditions in Y
		if( space()->getBCLocal()->getBCL(Y)>0 ) {
			Ordinal j = 1;
			for( Ordinal k=op_->getSR(Z); k<=op_->getER(Z); ++k )
				for( Ordinal i=op_->getSR(X); i<=op_->getER(X); ++i )
					y.at(i,j,k) =
						(1-omega)*x.at(i,j,k) + omega/op_->getC(Y,j,0) *( b.at(i,j,k) - op_->getC(Y,j,+1)*x.at(i,j+1,k) );

		}
		if( space()->getBCLocal()->getBCU(Y)>0 ) {
			Ordinal j = space()->nLoc(Y);
			for( Ordinal k=op_->getSR(Z); k<=op_->getER(Z); ++k )
				for( Ordinal i=op_->getSR(X); i<=op_->getER(X); ++i )
					y.at(i,j,k) =
						(1-omega)*x.at(i,j,k) + omega/op_->getC(Y,j,0) *( b.at(i,j,k) - op_->getC(Y,j,-1)*x.at(i,j-1,k) );
		}

		// boundary conditions in Z
		if( space()->getBCLocal()->getBCL(Z)>0 ) {
			Ordinal k = 1;
			for( Ordinal j=op_->getSR(Y); j<=op_->getER(Y); ++j )
				for( Ordinal i=op_->getSR(X); i<=op_->getER(X); ++i )
					y.at(i,j,k) =
						(1-omega)*x.at(i,j,k) + omega/op_->getC(Z,k,0) *( b.at(i,j,k) - op_->getC(Z,k,+1)*x.at(i,j,k+1) );
		}
		if( space()->getBCLocal()->getBCU(Z)>0 ) {
			Ordinal k = space()->nLoc(Z);
			for( Ordinal j=op_->getSR(Y); j<=op_->getER(Y); ++j )
				for( Ordinal i=op_->getSR(X); i<=op_->getER(X); ++i )
					y.at(i,j,k) =
						(1-omega)*x.at(i,j,k) + omega/op_->getC(Z,k,0) *( b.at(i,j,k) - op_->getC(Z,k,-1)*x.at(i,j,k-1) );
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
