#pragma once
#ifndef PIMPACT_TEUCHOSUTILS_HPP
#define PIMPACT_TEUCHOSUTILS_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialQRDenseSolver.hpp"

#include "Pimpact_DivGradO2Op.hpp"
#include "Pimpact_ScalarField.hpp"




namespace Pimpact{




/// \brief helper class to transfer ScalarField's to Teuchos::SerialDenseVector
/// and DivGradO2Op to Teuchos::SerialDenseMatrix 
///
/// \note it would be wishful to have a more general apply function for general operators
/// 
/// \tparam SpaceT Space 
template<class SpaceT>
class TeuchosTransfer {

public:

	using DomainFieldT= ScalarField<SpaceT>;
	using RangeFieldT = ScalarField<SpaceT>;

	using Scalar  = typename SpaceT::Scalar;
	using Ordinal = typename SpaceT::Ordinal;

	using VectorT = Teuchos::SerialDenseVector<Ordinal,Scalar>;
	using MatrixT = Teuchos::SerialDenseMatrix<Ordinal,Scalar>; /// \todo having sparse matrix

protected:

	using TO = const Teuchos::Tuple<Ordinal,SpaceT::sdim>;

	const Teuchos::RCP<const SpaceT> space_;

	Ordinal N_;

	TO SS_;
	TO NN_;

	Teuchos::Tuple<Ordinal,SpaceT::sdim-1> cw_;

	constexpr Ordinal getI( const Ordinal& i, const Ordinal& j, const Ordinal& k ) {
		return(
				(SpaceT::sdim==2)?
				(i-SS_[X]) + (j-SS_[Y])*cw_[0] :
				(i-SS_[X]) + (j-SS_[Y])*cw_[0] + (k-SS_[Z])*cw_[1]
				);
	}

public:

	constexpr const Ordinal& getN() { return( N_ ); }

	TeuchosTransfer(
			const Teuchos::RCP<const SpaceT>& space,
			const Ordinal* SS,
			const Ordinal* NN ):
		space_(space) {

			for( int i=0; i<SpaceT::sdim; ++i ) {
				SS_[i] = SS[i];
				NN_[i] = NN[i]; 
			}

			N_ = 1;

			for( int i=0; i<SpaceT::sdim; ++i ) 
				N_ *= ( NN_[i] - SS_[i] + 1 );

			cw_[0] = ( NN_[0] - SS_[0] + 1 );
			if( 3==SpaceT::sdim )
				cw_[1] = ( NN_[1] - SS_[1] + 1 )*cw_[0];
		}

	TeuchosTransfer( const Teuchos::RCP<const SpaceT>& space ):
		TeuchosTransfer( space, space->sInd(F::S), space->eInd(F::S) ) {}


	void apply( const DomainFieldT& x, VectorT& v ) const {

		assert( N_== v.numRows() * v.numCols() );

		if( 2==SpaceT::sdim )
			for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j )
				for( Ordinal i=SS_[X]; i<=NN_[X]; ++i )
						v( getI(i,j,1) ) = x(i,j,1);
		else 
			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k )
				for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j )
					for( Ordinal i=SS_[X]; i<=NN_[X]; ++i )
						v( getI(i,j,k) ) = x(i,j,k);
	}



	void apply( const VectorT& v, DomainFieldT& x ) const {

		assert( N_ == v.numRows()*v.numCols() );

		if( 2==SpaceT::sdim )
			for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j )
				for( Ordinal i=SS_[X]; i<=NN_[X]; ++i )
					x(i,j,1) = v( getI(i,j,1) );
		else
			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k )
				for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j )
					for( Ordinal i=SS_[X]; i<=NN_[X]; ++i )
						x(i,j,k) = v( getI(i,j,k) );
		//x.print();
	}


	void apply( const VectorT& v, DomainFieldT& x, const Scalar& om ) const {

		assert( N_==v.numRows()*v.numCols() );

		if( 2==SpaceT::sdim )
			for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j )
				for( Ordinal i=SS_[X]; i<=NN_[X]; ++i )
					x(i,j,1) = (1-om)*x(i,j,1) + om*v( getI(i,j,1) );
		else
			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k )
				for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j )
					for( Ordinal i=SS_[X]; i<=NN_[X]; ++i )
						x(i,j,k) = (1-om)*x(i,j,k) + om*v( getI(i,j,k) );
	}


	void apply( const Teuchos::RCP<const DivGradO2Op<SpaceT> >& op,
			Teuchos::RCP<MatrixT> A ) const {

		if( A.is_null() )
			A = Teuchos::rcp( new MatrixT( N_, N_, true ) );
		else
			*A = 0.;

		assert( N_*N_==A->numRows()*A->numCols() );

		// inner stencil
		for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k )
			for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j )
				for( Ordinal i=SS_[X]; i<=NN_[X]; ++i ) {

					const bool bcX = (op->space()->getBCLocal()->getBCL(X) > 0 && i==op->space()->begin(F::S,X) ) ||
						(op->space()->getBCLocal()->getBCU(X) > 0 && i==op->space()->end(F::S,X) ) ;
					const bool bcY = (op->space()->getBCLocal()->getBCL(Y) > 0 && j==op->space()->begin(F::S,Y) ) ||
						(op->space()->getBCLocal()->getBCU(Y) > 0 && j==op->space()->end(F::S,Y) ) ;
					const bool bcZ = (op->space()->getBCLocal()->getBCL(Z) > 0 && k==op->space()->begin(F::S,Z) ) ||
						(op->space()->getBCLocal()->getBCU(Z) > 0 && k==op->space()->end(F::S,Z) ) ;

					const Scalar& eps = 1.e-1;

					const Scalar epsX = (bcY||bcZ)?eps:1.;
					const Scalar epsY = (bcX||bcZ)?eps:1.;
					const Scalar epsZ = (bcX||bcY)?eps:1.;

					const Ordinal I = getI( i, j, k );

					for( int o=-1; o<=1; ++o ) {
						if( (i+o)>=SS_[X] && (i+o)<=NN_[X] ) 
							(*A)( I, getI( i+o, j, k ) ) += epsX*op->getC( X, i, o) ;

						if( (j+o)>=SS_[Y] && (j+o)<=NN_[Y] )
							(*A)( I, getI( i, j+o, k ) ) += epsY*op->getC( Y, j, o) ;

						if( 3==SpaceT::sdim && (k+o)>=SS_[Z] && (k+o)<=NN_[Z] )
							(*A)( I, getI( i, j, k+o ) ) += epsZ*op->getC( Z, k, o) ;
					}
				}
	}


	// \note not BC ready
	void updateRHS( const Teuchos::RCP<const DivGradO2Op<SpaceT> > op,
			const DomainFieldT& x, Teuchos::RCP<VectorT> b ) const {

		assert( N_==b->numRows()*b->numCols() );

		// boundary conditions in X
		if( space_->getBCLocal()->getBCL(X)<=0 || SS_[X]>space_->begin(F::S,X,true) ) {

			Ordinal i = SS_[X];

			Ordinal SZ = SS_[Z];
			if( space_->getBCLocal()->getBCL(Z)>0 && SS_[Z]==space_->begin(F::S,Z,true) )
				++SZ;
			Ordinal NZ = NN_[Z];
			if( space_->getBCLocal()->getBCU(Z)>0 && NN_[Z]==space_->end(F::S,Z,true) )
				--NZ;
			Ordinal SY = SS_[Y];
			if( space_->getBCLocal()->getBCL(Y)>0 && SS_[Y]==space_->begin(F::S,Y,true) )
				++SY;
			Ordinal NY = NN_[Y];
			if( space_->getBCLocal()->getBCU(Y)>0 && NN_[Y]==space_->end(F::S,Y,true) )
				--NY;

			for( Ordinal k=SZ; k<=NZ; ++k )
				for( Ordinal j=SY; j<=NY; ++j )
					(*b)( getI(i,j,k) ) -= x(i-1,j,k)*op->getC( X, i, -1 );
		}

		if( space_->getBCLocal()->getBCU(X)<=0 || NN_[X]<space_->end(F::S,X,true) ) {

			Ordinal i = NN_[X];

			Ordinal SZ = SS_[Z];
			if( space_->getBCLocal()->getBCL(Z)>0 && SS_[Z]==space_->begin(F::S,Z,true) )
				++SZ;
			Ordinal NZ = NN_[Z];
			if( space_->getBCLocal()->getBCU(Z)>0 && NN_[Z]==space_->end(F::S,Z,true) )
				--NZ;
			Ordinal SY = SS_[Y];
			if( space_->getBCLocal()->getBCL(Y)>0 && SS_[Y]==space_->begin(F::S,Y,true) )
				++SY;
			Ordinal NY = NN_[Y];
			if( space_->getBCLocal()->getBCU(Y)>0 && NN_[Y]==space_->end(F::S,Y,true) )
				--NY;

			for( Ordinal k=SZ; k<=NZ; ++k )
				for( Ordinal j=SY; j<=NY; ++j )
					(*b)( getI(i,j,k) ) -= x(i+1,j,k)*op->getC( X, i, +1 );
		}

		// boundary conditions in Y
		if( space_->getBCLocal()->getBCL(Y)<=0 || SS_[Y]>space_->begin(F::S,Y,true) ) {

			Ordinal j = SS_[Y];

			Ordinal SZ = SS_[Z];
			if( space_->getBCLocal()->getBCL(Z)>0 && SS_[Z]==space_->begin(F::S,Z,true) )
				++SZ;
			Ordinal NZ = NN_[Z];
			if( space_->getBCLocal()->getBCU(Z)>0 && NN_[Z]==space_->end(F::S,Z,true) )
				--NZ;
			Ordinal SX = SS_[X];
			if( space_->getBCLocal()->getBCL(X)>0 && SS_[X]==space_->begin(F::S,X,true) )
				++SX;
			Ordinal NX = NN_[X];
			if( space_->getBCLocal()->getBCU(X)>0 && NN_[X]==space_->end(F::S,X,true) )
				--NX;

			for( Ordinal k=SZ; k<=NZ; ++k )
				for( Ordinal i=SX; i<=NX; ++i )
					(*b)( getI(i,j,k) ) -= x(i,j-1,k)*op->getC( Y, j, -1 );
		}

		if( space_->getBCLocal()->getBCU(Y)<=0 || NN_[Y]<space_->end(F::S,Y,true) ) {

			Ordinal j = NN_[Y];

			Ordinal SZ = SS_[Z];
			if( space_->getBCLocal()->getBCL(Z)>0 && SS_[Z]==space_->begin(F::S,Z,true) )
				++SZ;
			Ordinal NZ = NN_[Z];
			if( space_->getBCLocal()->getBCU(Z)>0 && NN_[Z]==space_->end(F::S,Z,true) )
				--NZ;
			Ordinal SX = SS_[X];
			if( space_->getBCLocal()->getBCL(X)>0 && SS_[X]==space_->begin(F::S,X,true) )
				++SX;
			Ordinal NX = NN_[X];
			if( space_->getBCLocal()->getBCU(X)>0 && NN_[X]==space_->end(F::S,X,true) )
				--NX;

			for( Ordinal k=SZ; k<=NZ; ++k )
				for( Ordinal i=SX; i<=NX; ++i )
					(*b)( getI(i,j,k) ) -= x(i,j+1,k)*op->getC( Y, j, +1 );
		}

		// boundary conditions in Z
		if( space_->getBCLocal()->getBCL(Z)<=0 || SS_[Z]>space_->begin(F::S,Z,true) ) {

			Ordinal k = SS_[Z];

			Ordinal SX = SS_[X];
			if( space_->getBCLocal()->getBCL(X)>0 && SS_[X]==space_->begin(F::S,X,true) )
				++SX;
			Ordinal NX = NN_[X];
			if( space_->getBCLocal()->getBCU(X)>0 && NN_[X]==space_->end(F::S,X,true) )
				--NX;
			Ordinal SY = SS_[Y];
			if( space_->getBCLocal()->getBCL(Y)>0 && SS_[Y]==space_->begin(F::S,Y,true) )
				++SY;
			Ordinal NY = NN_[Y];
			if( space_->getBCLocal()->getBCU(Y)>0 && NN_[Y]==space_->end(F::S,Y,true) )
				--NY;

			for( Ordinal j=SY; j<=NY; ++j )
				for( Ordinal i=SX; i<=NX; ++i )
					(*b)( getI(i,j,k) ) -= x(i,j,k-1)*op->getC( Z, k, -1 );
		}

		if( space_->getBCLocal()->getBCU(Z)<=0 || NN_[Z]<space_->end(F::S,Z,true) ) {

			Ordinal k = NN_[Z];

			Ordinal SX = SS_[X];
			if( space_->getBCLocal()->getBCL(X)>0 && SS_[X]==space_->begin(F::S,X,true) )
				++SX;
			Ordinal NX = NN_[X];
			if( space_->getBCLocal()->getBCU(X)>0 && NN_[X]==space_->end(F::S,X,true) )
				--NX;
			Ordinal SY = SS_[Y];
			if( space_->getBCLocal()->getBCL(Y)>0 && SS_[Y]==space_->begin(F::S,Y,true) )
				++SY;
			Ordinal NY = NN_[Y];
			if( space_->getBCLocal()->getBCU(Y)>0 && NN_[Y]==space_->end(F::S,Y,true) )
				--NY;

			for( Ordinal j=SY; j<=NY; ++j )
				for( Ordinal i=SX; i<=NX; ++i )
					(*b)( getI(i,j,k) ) -= x(i,j,k+1)*op->getC( Z, k, +1 );
		}
	}

	constexpr const std::string getLabel() const { return( "TeuchosTransfer" ); };

	void print( std::ostream& out=std::cout ) const {
		out << "--- " << getLabel() << " ---\n";
		out << "N: " << N_ << "\n";
		out << "cw: " << cw_ << "\n";
		out << "SS: \t";
		 out << SS_ << "\n";
		out << "NN: \t";
		out << NN_ << "\n";
		out << "\n";
	}

}; // end of class TeuchosTransfer



template<class OperatorT>
class TeuchosSolver {

public:

	using SpaceT = typename OperatorT::SpaceT;

	using Scalar  = typename SpaceT::Scalar;
	using Ordinal = typename SpaceT::Ordinal;

	using DomainFieldT= typename OperatorT::DomainFieldT;
	using RangeFieldT = typename OperatorT::RangeFieldT;

	using VectorT = Teuchos::SerialDenseVector<Ordinal,Scalar>;
	using MatrixT = Teuchos::SerialDenseMatrix<Ordinal,Scalar>;

	//using SolverT = Teuchos::SerialDenseSolver<Ordinal,Scalar>;
	using SolverT = Teuchos::SerialQRDenseSolver<Ordinal,Scalar>;

protected:

	Teuchos::RCP< const OperatorT > op_;
	Teuchos::RCP< TeuchosTransfer<SpaceT> > trans_;
	Teuchos::RCP< SolverT > Asov_;


	void init() {

		Teuchos::RCP<MatrixT> A = Teuchos::rcp( new MatrixT( getN(), getN(), true ) );
		Asov_ = Teuchos::rcp( new SolverT() );

		trans_->apply( op_, A );

		// set solver
		Asov_->factorWithEquilibration( true );
		Asov_->setMatrix( A );
		Asov_->factor();
	}

public:

	TeuchosSolver( const Teuchos::RCP<const OperatorT>& op ):
		op_( op ),
		trans_( Teuchos::rcp( new TeuchosTransfer<SpaceT>( op->space() ) ) ) {

			init();
		}

	TeuchosSolver(
			const Teuchos::RCP<const OperatorT>& op,
			const Ordinal* SS,
			const Ordinal* NN ):
		op_(op),
		trans_( Teuchos::rcp( new TeuchosTransfer<SpaceT>( op->space(), SS, NN ) ) ) {

			init();
		}

	void apply( const DomainFieldT& x, RangeFieldT& y )  const {

		VectorT X_( getN(), false );
		VectorT	B_( getN(), false );

		trans_->apply( x, B_ );

		y.exchange();
		//trans_->updateRHS( op_, y, B_ );

		Asov_->setVectors( Teuchos::rcpFromRef(X_), Teuchos::rcpFromRef(B_) );
		Asov_->solve();

		trans_->apply( X_, y );

		y.changed();
	}

	void apply( const DomainFieldT& x, RangeFieldT& y, const Scalar& omega )  const {

		VectorT X_( getN(), false );
		VectorT	B_( getN(), false );

		trans_->apply( x, B_ );

		Asov_->setVectors( X_, B_ );
		Asov_->solve();

		trans_->apply( X_, y, omega );

		y.changed();
	}

	/// \name getter
	/// @{ 
	
protected:

	constexpr const Ordinal& getN() const { return( trans_->getN() ); }

public:

	constexpr const Teuchos::RCP< const OperatorT >& getOperator() const { return( op_ ); }
	constexpr const Teuchos::RCP< TeuchosTransfer<SpaceT> >& getTeuchosTransfer() const { return( trans_ ); }

	///  @} 
};



template<class OperatorT>
class TeuchosEigenvalues {

	using SpaceT = typename OperatorT::SpaceT;

	using Scalar  = typename SpaceT::Scalar;
	using Ordinal = typename SpaceT::Ordinal;

	using DomainFieldT= typename OperatorT::DomainFieldT;
	using RangeFieldT = typename OperatorT::RangeFieldT;

	using VectorT = Teuchos::SerialDenseVector<Ordinal,Scalar>;
	using MatrixT = Teuchos::SerialDenseMatrix<Ordinal,Scalar>;

protected:

	Teuchos::RCP< const OperatorT > op_;
	Teuchos::RCP< TeuchosTransfer<SpaceT> > trans_;

public:

	TeuchosEigenvalues( const Teuchos::RCP<const OperatorT>& op ):
		op_( op ),
		trans_( Teuchos::rcp( new TeuchosTransfer<SpaceT>( op->space() ) ) ) { }

	/// \name getter
	/// @{ 
	
protected:

	const Ordinal& getN() const { return( trans_->getN() ); }

public:

	Teuchos::RCP< const OperatorT > getOperator() const { return( op_ ); }
	Teuchos::RCP< TeuchosTransfer<SpaceT> > getTeuchosTransfer() const { return( trans_ ); }

	///  @} 
	//
	void computeEV( const ECoord& dir, Scalar& evMax, Scalar& evMin ) const {


		Ordinal n = op_->space()->nLoc(dir)-1+1;

		Teuchos::RCP< MatrixT > A;

		A = Teuchos::rcp( new MatrixT( n, n, true ) );
		// construct A
		for( Ordinal i=1; i<=op_->space()->nLoc(dir); ++i )
			for( int o=-1; o<=1; ++o ) {
				Ordinal I = i-1;
				if( (i+o)>=op_->space()->begin(F::S,dir) && (i+o)<=op_->space()->end(F::S,dir) ) 
					(*A)( I, I+o ) += op_->getC( dir, i, o) ;
			}


		// lower boundaries // deperacted?
		if( op_->space()->getBCLocal()->getBCL(dir)>0 ) {
			Ordinal i = 1;
			Ordinal I = i-1;
			// set row zero
			for( Ordinal l=0; l<A->numCols(); ++l ) 
				(*A)(I,l) = 0.;
			(*A)( I, I )   += op_->getC( dir, i, 0) ;
			(*A)( I, I+1 ) += op_->getC( dir, i,+1) ;
		}

		// upper boundaries
		if( op_->space()->getBCLocal()->getBCU(dir)>0 ) {
			Ordinal i = op_->space()->nLoc(dir);
			Ordinal I = i-1;
			// set row zero
			for( Ordinal l=0; l<A->numCols(); ++l ) 
				(*A)(I,l) = 0.;
			(*A)( I, I )   += op_->getC( dir, i, 0) ;
			(*A)( I, I-1 ) += op_->getC( dir, i,-1) ;
		}

		// compute EV, todo move to TeuchosBla
		Ordinal sdim = 0;
		Teuchos::RCP< VectorT > evr = Teuchos::rcp( new VectorT(n,false) );
		Teuchos::RCP< VectorT > evi = Teuchos::rcp( new VectorT(n,false) );
		Teuchos::RCP< VectorT > work = Teuchos::rcp( new VectorT(4*n,false) );
		Teuchos::RCP< VectorT > rwork = Teuchos::rcp( new VectorT(3*n,false) );

		Teuchos::ArrayRCP<Ordinal> bwork = Teuchos::arcp<Ordinal>( n );
		Ordinal info = 0;

		Teuchos::LAPACK<Ordinal,Scalar> lapack;
		lapack.GEES(
				'N', 						// Schur vectors are not computed
				n,              // The order of the matrix A. N >= 0
				A->values(),    // A is DOUBLE PRECISION array, dimension (LDA,N); On entry, the N-by-N matrix A.; On exit, A has been overwritten by its real Schur form T.;
				A->stride(),   // The leading dimension of the array A.
				&sdim,           // If SORT = 'N', SDIM = 0.; If SORT = 'S', SDIM = number of eigenvalues (after sorting); for which SELECT is true. (Complex conjugate pairs for which SELECT is true for either eigenvalue count as 2.)
				evr->values(),  // 
				evi->values(),
				0,              // If JOBVS = 'N', VS is not referenced.
				1,              // The leading dimension of the array VS.
				work->values(), // 
				4*n,           //
				//-1,           //
				rwork->values(), //
				bwork.getRawPtr(),
				&info );

		//std::cout << "info: " << info << "\n";
		//std::cout << "opti work: " << (*work)[0]/N_ << "\n";

		Scalar eps = 1.e-10;

		Teuchos::RCP<std::ostream> out = Pimpact::createOstream( "ev_"+toString(dir)+".txt" );
		for( Ordinal i=0; i<n; ++i ) {
			*out << (*evr)[i] << "\n";
			if( 0==i ) {
				for( Ordinal j=0; j<n; ++ j ) {
					if( std::fabs( (*evr)[j] )>eps ) {
						evMax = (*evr)[j] ;
						j=n;
					}
				}
				for( Ordinal j=0; j<n; ++ j ) {
					if( std::fabs( (*evr)[j] )>eps ) {
						evMin = (*evr)[j] ;
						j=n;
					}
				}
			}
			else {
				if( std::fabs( (*evr)[i] )>eps ) 
					evMax = std::max( evMax, (*evr)[i] );
				if( std::fabs( (*evr)[i] )>eps ) 
					evMin = std::min( evMin, (*evr)[i] );
			}
		}

			/*
			 * = 0: successful exit
       *    < 0: if INFO = -i, the i-th argument had an illegal value.
       *    > 0: if INFO = i, and i is
       *       <= N: the QR algorithm failed to compute all the
       *             eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI
       *             contain those eigenvalues which have converged; if
       *             JOBVS = 'V', VS contains the matrix which reduces A
       *             to its partially converged Schur form.
       *       = N+1: the eigenvalues could not be reordered because some
       *             eigenvalues were too close to separate (the problem
       *             is very ill-conditioned);
       *       = N+2: after reordering, roundoff changed values of some
       *             complex eigenvalues so that leading eigenvalues in
       *             the Schur form no longer satisfy SELECT=.TRUE.  This
       *             could also be caused by underflow due to scaling.
			 */
	}

	void computeEV( Scalar& evMax, Scalar& evMin ) const {

		evMax = 0.;
		evMin = 0.;
		
		for( int i=0; i<SpaceT::sdim; ++i ) {
			Scalar evMaxTemp, evMinTemp;
			computeEV( static_cast<Pimpact::ECoord>(i), evMaxTemp, evMinTemp );
			evMax += evMaxTemp;
			evMin += evMinTemp;
		}
	}

	// move to EV computer
	void computeFullEV( Scalar& evMax, Scalar& evMin ) const {

		Teuchos::RCP<const TeuchosTransfer<SpaceT> > trans =
			Teuchos::rcp( new TeuchosTransfer<SpaceT>(
						op_->space(),
						op_->space()->sIndB(F::S),
						op_->space()->eIndB(F::S) ) );

		Ordinal N = trans->getN();

		Teuchos::RCP<MatrixT>	A = Teuchos::rcp( new MatrixT( N, N, true ) );

		trans->apply( op_, A );
		{
			Teuchos::RCP<std::ostream> out = Pimpact::createOstream( "divGrad.txt" );
			*out << *A;
		}

		
		// compute EV, todo move to TeuchosBla
		Ordinal sdim = 0;
		Teuchos::RCP< VectorT > evr = Teuchos::rcp( new VectorT(N,false) );
		Teuchos::RCP< VectorT > evi = Teuchos::rcp( new VectorT(N,false) );
		Teuchos::RCP< VectorT > work = Teuchos::rcp( new VectorT(4*N,false) );
		Teuchos::RCP< VectorT > rwork = Teuchos::rcp( new VectorT(3*N,false) );

		Teuchos::ArrayRCP<Ordinal> bwork = Teuchos::arcp<Ordinal>( N );
		Ordinal info = 0;

		Teuchos::LAPACK<Ordinal,Scalar> lapack;
		lapack.GEES(
				'N', 						// Schur vectors are not computed
				trans->getN(), // The order of the matrix A. N >= 0
				A->values(),   // A is DOUBLE PRECISION array, dimension (LDA,N); On entry, the N-by-N matrix A.; On exit, A has been overwritten by its real Schur form T.;
				A->stride(),   // The leading dimension of the array A.
				&sdim,           // If SORT = 'N', SDIM = 0.; If SORT = 'S', SDIM = number of eigenvalues (after sorting); for which SELECT is true. (Complex conjugate pairs for which SELECT is true for either eigenvalue count as 2.)
				evr->values(),  // 
				evi->values(),
				0,              // If JOBVS = 'N', VS is not referenced.
				1,              // The leading dimension of the array VS.
				work->values(), // 
				4*N,           //
				//-1,           //
				rwork->values(), //
				bwork.getRawPtr(),
				&info );

		std::cout << "info: " << info << "\n";
		//std::cout << "opti work: " << (*work)[0]/N << "\n";

		Teuchos::RCP<std::ostream> out = Pimpact::createOstream( "evfull.txt" );
		for( Ordinal i=0; i<N; ++i ) {
			*out << (*evr)[i] << "\t" << (*evi)[i] << "\n";
			if( 0==i ) {
				evMax = (*evr)[i] ;
				evMin = (*evr)[i] ;
			}
			else {
				evMax = std::max( evMax, (*evr)[i] );
				evMin = std::min( evMin, (*evr)[i] );
			}
		}

			/*
			 * = 0: successful exit
       *    < 0: if INFO = -i, the i-th argument had an illegal value.
       *    > 0: if INFO = i, and i is
       *       <= N: the QR algorithm failed to compute all the
       *             eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI
       *             contain those eigenvalues which have converged; if
       *             JOBVS = 'V', VS contains the matrix which reduces A
       *             to its partially converged Schur form.
       *       = N+1: the eigenvalues could not be reordered because some
       *             eigenvalues were too close to separate (the problem
       *             is very ill-conditioned);
       *       = N+2: after reordering, roundoff changed values of some
       *             complex eigenvalues so that leading eigenvalues in
       *             the Schur form no longer satisfy SELECT=.TRUE.  This
       *             could also be caused by underflow due to scaling.
			 */
	}
};


} // end of namespace Pimpact



#ifdef COMPILE_ETI
extern template class Pimpact::TeuchosTransfer< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::TeuchosTransfer< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::TeuchosTransfer< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::TeuchosTransfer< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_TEUCHOSUTILS_HPP
