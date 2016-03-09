#pragma once
#ifndef PIMPACT_TEUCHOSTRANSFER_HPP
#define PIMPACT_TEUCHOSTRANSFER_HPP


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

	using Scalar  = typename SpaceT::Scalar;
	using Ordinal = typename SpaceT::Ordinal;

	using DomainFieldT= ScalarField<SpaceT>;
	using RangeFieldT = ScalarField<SpaceT>;

protected:

	using VectorT = Teuchos::SerialDenseVector<Ordinal,Scalar>;
	using MatrixT = Teuchos::SerialDenseMatrix<Ordinal,Scalar>;

	Teuchos::RCP<const SpaceT> space_;

	Ordinal N_;

	const Ordinal* SS_;
	const Ordinal* NN_;

	Teuchos::Tuple<Ordinal,3> cw_;

	Ordinal getI( const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {
		return( (i-SS_[X]) + (j-SS_[Y])*cw_[0] + (k-SS_[Z])*cw_[0]*cw_[1] );
	}

public:

	const Ordinal& getN() const { return( N_ ); }

	TeuchosTransfer(
			const Teuchos::RCP<const SpaceT>& space,
			const Ordinal* SS, const Ordinal* NN ):
		space_(space),SS_(SS),NN_(NN) {

			N_ = 1;

			for( int i=0; i<3; ++i ) {
				N_ *= ( NN_[i] - SS_[i] + 1 );
				cw_[i] = ( NN_[i] - SS_[i] + 1 );
			}

		}



	void apply( const DomainFieldT& x, Teuchos::RCP<VectorT> v ) const {

		if( v.is_null() )
			v = Teuchos::rcp( new VectorT( N_, false ) );

		TEUCHOS_TEST_FOR_EXCEPT( N_!= v->numRows() * v->numCols() )

			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k )
				for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j )
					for( Ordinal i=SS_[X]; i<=NN_[X]; ++i )
						(*v)( getI(i,j,k) ) = x.at(i,j,k);

	}



	void apply( const Teuchos::RCP<const VectorT>& v, DomainFieldT& x ) const {

		TEUCHOS_TEST_FOR_EXCEPT( N_!=v->numRows()*v->numCols() )

			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k )
				for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j )
					for( Ordinal i=SS_[X]; i<=NN_[X]; ++i )
						x.at(i,j,k) = (*v)( getI(i,j,k) );

	}


	void apply( const Teuchos::RCP<const VectorT>& v, DomainFieldT& x, const Scalar& om ) const {

		TEUCHOS_TEST_FOR_EXCEPT( N_!=v->numRows()*v->numCols() )

			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k )
				for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j )
					for( Ordinal i=SS_[X]; i<=NN_[X]; ++i )
						x.at(i,j,k) = (1-om)*x.at(i,j,k) + om*(*v)( getI(i,j,k) );

	}



	void apply( Teuchos::RCP<const DivGradO2Op<SpaceT> > op,
			Teuchos::RCP<MatrixT> A ) const {

		if( A.is_null() )
			A = Teuchos::rcp( new MatrixT( N_, N_, true ) );
		else
			*A = 0.;

		TEUCHOS_TEST_FOR_EXCEPT( N_*N_!=A->numRows()*A->numCols() );

		// inner stencil
		for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k )
			for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j )
				for( Ordinal i=SS_[X]; i<=NN_[X]; ++i ) {

					Ordinal I = getI( i, j, k );

					for( int o=-1; o<=1; ++o ) {
						if( (i+o)>=SS_[X] && (i+o)<=NN_[X] ) 
							(*A)( I, getI( i+o, j, k ) ) += op->getC( X, i, o) ;

						if( (j+o)>=SS_[Y] && (j+o)<=NN_[Y] )
							(*A)( I, getI( i, j+o, k ) ) += op->getC( Y, j, o) ;

						if( (k+o)>=SS_[Z] && (k+o)<=NN_[Z] )
							(*A)( I, getI( i, j, k+o ) ) += op->getC( Z, k, o) ;
					}

				}

		// boundary conditions
		if( space_->getBCLocal()->getBCL(X)>0 && SS_[X]==space_->sIndB(EField::S,X) ) {

			Ordinal i = space_->sIndB(EField::S,X);

			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k )
				for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j ) {

					Ordinal I = getI( i, j, k );

					// set row zero
					for( Ordinal l=0; l<A->numCols(); ++ l ) 
						(*A)(I,l) = 0.;

					for( int o=-1; o<=1; ++o )
						if( (i+o)>=SS_[X] && (i+o)<=NN_[X] )
							(*A)( I, getI( i+o, j, k ) ) = op->getC( X, i, o) ;
				}
		}

		if( space_->getBCLocal()->getBCU(X)>0 && NN_[X]==space_->eIndB(EField::S,X) ) {

			Ordinal i = space_->eIndB(EField::S,X);

			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k )
				for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j ) {

					Ordinal I = getI( i, j, k );

					// set row zero
					for( Ordinal l=0; l<A->numCols(); ++ l ) 
						(*A)(I,l) = 0.;

					// set stencil
					for( int o=-1; o<=1; ++o )
						if( (i+o)>=SS_[X] && (i+o)<=NN_[X] )
							(*A)( I, getI( i+o, j, k ) ) = op->getC( X, i, o) ;
				}
		}

		if( space_->getBCLocal()->getBCL(Y)>0 && SS_[Y]==space_->sIndB(EField::S,Y) ) {

			Ordinal j = space_->sIndB(EField::S,Y);

			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k )
				for( Ordinal i=SS_[X]; i<=NN_[X]; ++i ) {

					Ordinal I = getI( i, j, k );

					// set row zero
					for( Ordinal l=0; l<A->numCols(); ++ l ) 
						(*A)(I,l) = 0.;

					for( int o=-1; o<=1; ++o )
						if( (j+o)>=SS_[Y] && (j+o)<=NN_[Y] )
							(*A)( I, getI( i, j+o, k ) ) = op->getC( Y, j, o) ;
				}
		}

		if( space_->getBCLocal()->getBCU(Y)>0 && NN_[Y]==space_->eIndB(EField::S,Y) ) {

			Ordinal j = space_->eIndB(EField::S,Y);

			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k )
				for( Ordinal i=SS_[X]; i<=NN_[X]; ++i ) {

					Ordinal I = getI( i, j, k );

					// set row zero
					for( Ordinal l=0; l<A->numCols(); ++ l ) 
						(*A)(I,l) = 0.;

					for( int o=-1; o<=1; ++o )
						if( (j+o)>=SS_[Y] && (j+o)<=NN_[Y] )
							(*A)( I, getI( i, j+o, k ) ) = op->getC( Y, j, o);
				}
		}

		if( space_->getBCLocal()->getBCL(Z)>0 && SS_[Z]==space_->sIndB(EField::S,Z) ) {

			Ordinal k = space_->sIndB(EField::S,Z);

			for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j )
				for( Ordinal i=SS_[X]; i<=NN_[X]; ++i ) {

					Ordinal I = getI( i, j, k );

					// set row zero
					for( Ordinal l=0; l<A->numCols(); ++ l ) 
						(*A)(I,l) = 0.;

					for( int o=-1; o<=1; ++o ) {
						if( (k+o)>=SS_[Z] && (k+o)<=NN_[Z] )
							(*A)( I, getI( i, j, k+o ) ) += op->getC( Z, k, o) ;
					}

				}
		}

		if( space_->getBCLocal()->getBCU(Z)>0 && NN_[Z]==space_->eIndB(EField::S,Z) ) {

			Ordinal k = space_->eIndB(EField::S,Z);

			for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j )
				for( Ordinal i=SS_[X]; i<=NN_[X]; ++i ) {

					Ordinal I = getI( i, j, k );

					// set row zero
					for( Ordinal l=0; l<A->numCols(); ++ l ) 
						(*A)(I,l) = 0.;

					for( int o=-1; o<=1; ++o ) 
						if( (k+o)>=SS_[Z] && (k+o)<=NN_[Z] )
							(*A)( I, getI( i, j, k+o ) ) += op->getC( Z, k, o) ;
				}
		}

		// --- corners ---
		if( space_->getBCLocal()->getBCL(X)>0 && space_->getBCLocal()->getBCL(Y)>0 &&
				SS_[X]==space_->sIndB(EField::S,X) && SS_[Y]==space_->sIndB(EField::S,Y) ) {

			Ordinal i = space_->sIndB(EField::S,X);
			Ordinal j = space_->sIndB(EField::S,Y);

			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k ) {

				Ordinal I = getI( i, j, k );

				// set row zero
				for( Ordinal l=0; l<A->numCols(); ++ l ) (*A)(I,l) = 0.;

				(*A)(I,I) = 1.;
			}
		}
		if( space_->getBCLocal()->getBCL(X)>0 && space_->getBCLocal()->getBCU(Y)>0 &&
				SS_[X]==space_->sIndB(EField::S,X) && NN_[Y]==space_->eIndB(EField::S,Y) ) {

			Ordinal i = space_->sIndB(EField::S,X);
			Ordinal j = space_->eIndB(EField::S,Y);

			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k ) {

				Ordinal I = getI( i, j, k );

				// set row zero
				for( Ordinal l=0; l<A->numCols(); ++ l ) (*A)(I,l) = 0.;

				(*A)(I,I) = 1.;
			}
		}
		if( space_->getBCLocal()->getBCU(X)>0 && space_->getBCLocal()->getBCL(Y)>0 &&
				NN_[X]==space_->eIndB(EField::S,X) && SS_[Y]==space_->sIndB(EField::S,Y) ) {

			Ordinal i = space_->eIndB(EField::S,X);
			Ordinal j = space_->sIndB(EField::S,Y);

			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k ) {

				Ordinal I = getI( i, j, k );

				// set row zero
				for( Ordinal l=0; l<A->numCols(); ++ l ) (*A)(I,l) = 0.;

				(*A)(I,I) = 1.;
			}
		}
		if( space_->getBCLocal()->getBCU(X)>0 && space_->getBCLocal()->getBCU(Y)>0 &&
				NN_[X]==space_->eIndB(EField::S,X) && NN_[Y]==space_->eIndB(EField::S,Y) ) {

			Ordinal i = space_->eIndB(EField::S,X);
			Ordinal j = space_->eIndB(EField::S,Y);

			for( Ordinal k=SS_[Z]; k<=NN_[Z]; ++k ) {

				Ordinal I = getI( i, j, k );

				// set row zero
				for( Ordinal l=0; l<A->numCols(); ++ l ) (*A)(I,l) = 0.;

				(*A)(I,I) = 1.;
			}
		}

		if( space_->getBCLocal()->getBCL(X)>0 && space_->getBCLocal()->getBCL(Z)>0 &&
				SS_[X]==space_->sIndB(EField::S,X) && SS_[Z]==space_->sIndB(EField::S,Z) ) {

			Ordinal i = space_->sIndB(EField::S,X);
			Ordinal k = space_->sIndB(EField::S,Z);

			for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j ) {

				Ordinal I = getI( i, j, k );

				// set row zero
				for( Ordinal l=0; l<A->numCols(); ++ l ) (*A)(I,l) = 0.;

				(*A)(I,I) = 1.;
			}
		}
		if( space_->getBCLocal()->getBCL(X)>0 && space_->getBCLocal()->getBCU(Z)>0 &&
				SS_[X]==space_->sIndB(EField::S,X) && NN_[Z]==space_->eIndB(EField::S,Z) ) {

			Ordinal i = space_->sIndB(EField::S,X);
			Ordinal k = space_->eIndB(EField::S,Z);

			for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j ) {

				Ordinal I = getI( i, j, k );

				// set row zero
				for( Ordinal l=0; l<A->numCols(); ++ l ) (*A)(I,l) = 0.;

				(*A)(I,I) = 1.;
			}
		}
		if( space_->getBCLocal()->getBCU(X)>0 && space_->getBCLocal()->getBCL(Z)>0 &&
				NN_[X]==space_->eIndB(EField::S,X) && SS_[Z]==space_->sIndB(EField::S,Z) ) {

			Ordinal i = space_->eIndB(EField::S,X);
			Ordinal k = space_->sIndB(EField::S,Z);

			for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j ) {

				Ordinal I = getI( i, j, k );

				// set row zero
				for( Ordinal l=0; l<A->numCols(); ++ l ) (*A)(I,l) = 0.;

				(*A)(I,I) = 1.;
			}
		}
		if( space_->getBCLocal()->getBCU(X)>0 && space_->getBCLocal()->getBCU(Z)>0 &&
				NN_[X]==space_->eIndB(EField::S,X) && NN_[Z]==space_->eIndB(EField::S,Z) ) {

			Ordinal i = space_->eIndB(EField::S,X);
			Ordinal k = space_->eIndB(EField::S,Z);

			for( Ordinal j=SS_[Y]; j<=NN_[Y]; ++j ) {

				Ordinal I = getI( i, j, k );

				// set row zero
				for( Ordinal l=0; l<A->numCols(); ++ l ) (*A)(I,l) = 0.;

				(*A)(I,I) = 1.;
			}
		}
		if( space_->getBCLocal()->getBCL(Y)>0 && space_->getBCLocal()->getBCL(Z)>0 &&
				SS_[Y]==space_->sIndB(EField::S,Y) && SS_[Z]==space_->sIndB(EField::S,Z) ) {

			Ordinal j = space_->sIndB(EField::S,Y);
			Ordinal k = space_->sIndB(EField::S,Z);

			for( Ordinal i=SS_[X]; i<=NN_[X]; ++i ) {

				Ordinal I = getI( i, j, k );

				// set row zero
				for( Ordinal l=0; l<A->numCols(); ++ l ) (*A)(I,l) = 0.;

				(*A)(I,I) = 1.;
			}
		}
		if( space_->getBCLocal()->getBCL(Y)>0 && space_->getBCLocal()->getBCU(Z)>0 &&
				SS_[Y]==space_->sIndB(EField::S,Y) && NN_[Z]==space_->eIndB(EField::S,Z) ) {

			Ordinal j = space_->sIndB(EField::S,Y);
			Ordinal k = space_->eIndB(EField::S,Z);

			for( Ordinal i=SS_[X]; i<=NN_[X]; ++i ) {

				Ordinal I = getI( i, j, k );

				// set row zero
				for( Ordinal l=0; l<A->numCols(); ++ l ) (*A)(I,l) = 0.;

				(*A)(I,I) = 1.;
			}
		}
		if( space_->getBCLocal()->getBCU(Y)>0 && space_->getBCLocal()->getBCL(Z)>0 &&
				SS_[Y]==space_->eIndB(EField::S,Y) && NN_[Z]==space_->sIndB(EField::S,Z) ) {

			Ordinal j = space_->eIndB(EField::S,Y);
			Ordinal k = space_->sIndB(EField::S,Z);

			for( Ordinal i=SS_[X]; i<=NN_[X]; ++i ) {

				Ordinal I = getI( i, j, k );

				// set row zero
				for( Ordinal l=0; l<A->numCols(); ++ l ) (*A)(I,l) = 0.;

				(*A)(I,I) = 1.;
			}
		}
		if( space_->getBCLocal()->getBCU(Y)>0 && space_->getBCLocal()->getBCU(Z)>0 &&
				NN_[Y]==space_->eIndB(EField::S,Y) && NN_[Z]==space_->eIndB(EField::S,Z) ) {

			Ordinal j = space_->eIndB(EField::S,Y);
			Ordinal k = space_->eIndB(EField::S,Z);

			for( Ordinal i=SS_[X]; i<=NN_[X]; ++i ) {

				Ordinal I = getI( i, j, k );

				// set row zero
				for( Ordinal l=0; l<A->numCols(); ++ l ) (*A)(I,l) = 0.;

				(*A)(I,I) = 1.;
			}
		}

	}



	void updateRHS( const Teuchos::RCP<const DivGradO2Op<SpaceT> > op,
			const DomainFieldT& x, Teuchos::RCP<VectorT> b ) const {


		TEUCHOS_TEST_FOR_EXCEPT( N_!=b->numRows()*b->numCols() );


		// boundary conditions in X
		if( space_->getBCLocal()->getBCL(X)<=0 || SS_[X]>space_->sIndB(EField::S,X) ) {

			Ordinal i = SS_[X];

			Ordinal SZ = SS_[Z];
			if( space_->getBCLocal()->getBCL(Z)>0 && SS_[Z]==space_->sIndB(EField::S,Z) )
				++SZ;
			Ordinal NZ = NN_[Z];
			if( space_->getBCLocal()->getBCU(Z)>0 && NN_[Z]==space_->eIndB(EField::S,Z) )
				--NZ;
			Ordinal SY = SS_[Y];
			if( space_->getBCLocal()->getBCL(Y)>0 && SS_[Y]==space_->sIndB(EField::S,Y) )
				++SY;
			Ordinal NY = NN_[Y];
			if( space_->getBCLocal()->getBCU(Y)>0 && NN_[Y]==space_->eIndB(EField::S,Y) )
				--NY;

			for( Ordinal k=SZ; k<=NZ; ++k )
				for( Ordinal j=SY; j<=NY; ++j )
					(*b)( getI(i,j,k) ) -= x.at(i-1,j,k)*op->getC( X, i, -1 );
		}

		if( space_->getBCLocal()->getBCU(X)<=0 || NN_[X]<space_->eIndB(EField::S,X) ) {

			Ordinal i = NN_[X];

			Ordinal SZ = SS_[Z];
			if( space_->getBCLocal()->getBCL(Z)>0 && SS_[Z]==space_->sIndB(EField::S,Z) )
				++SZ;
			Ordinal NZ = NN_[Z];
			if( space_->getBCLocal()->getBCU(Z)>0 && NN_[Z]==space_->eIndB(EField::S,Z) )
				--NZ;
			Ordinal SY = SS_[Y];
			if( space_->getBCLocal()->getBCL(Y)>0 && SS_[Y]==space_->sIndB(EField::S,Y) )
				++SY;
			Ordinal NY = NN_[Y];
			if( space_->getBCLocal()->getBCU(Y)>0 && NN_[Y]==space_->eIndB(EField::S,Y) )
				--NY;

			for( Ordinal k=SZ; k<=NZ; ++k )
				for( Ordinal j=SY; j<=NY; ++j )
					(*b)( getI(i,j,k) ) -= x.at(i+1,j,k)*op->getC( X, i, +1 );
		}

		// boundary conditions in Y
		if( space_->getBCLocal()->getBCL(Y)<=0 || SS_[Y]>space_->sIndB(EField::S,Y) ) {

			Ordinal j = SS_[Y];

			Ordinal SZ = SS_[Z];
			if( space_->getBCLocal()->getBCL(Z)>0 && SS_[Z]==space_->sIndB(EField::S,Z) )
				++SZ;
			Ordinal NZ = NN_[Z];
			if( space_->getBCLocal()->getBCU(Z)>0 && NN_[Z]==space_->eIndB(EField::S,Z) )
				--NZ;
			Ordinal SX = SS_[X];
			if( space_->getBCLocal()->getBCL(X)>0 && SS_[X]==space_->sIndB(EField::S,X) )
				++SX;
			Ordinal NX = NN_[X];
			if( space_->getBCLocal()->getBCU(X)>0 && NN_[X]==space_->eIndB(EField::S,X) )
				--NX;

			for( Ordinal k=SZ; k<=NZ; ++k )
				for( Ordinal i=SX; i<=NX; ++i )
					(*b)( getI(i,j,k) ) -= x.at(i,j-1,k)*op->getC( Y, j, -1 );
		}

		if( space_->getBCLocal()->getBCU(Y)<=0 || NN_[Y]<space_->eIndB(EField::S,Y) ) {

			Ordinal j = NN_[Y];

			Ordinal SZ = SS_[Z];
			if( space_->getBCLocal()->getBCL(Z)>0 && SS_[Z]==space_->sIndB(EField::S,Z) )
				++SZ;
			Ordinal NZ = NN_[Z];
			if( space_->getBCLocal()->getBCU(Z)>0 && NN_[Z]==space_->eIndB(EField::S,Z) )
				--NZ;
			Ordinal SX = SS_[X];
			if( space_->getBCLocal()->getBCL(X)>0 && SS_[X]==space_->sIndB(EField::S,X) )
				++SX;
			Ordinal NX = NN_[X];
			if( space_->getBCLocal()->getBCU(X)>0 && NN_[X]==space_->eIndB(EField::S,X) )
				--NX;

			for( Ordinal k=SZ; k<=NZ; ++k )
				for( Ordinal i=SX; i<=NX; ++i )
					(*b)( getI(i,j,k) ) -= x.at(i,j+1,k)*op->getC( Y, j, +1 );
		}

		// boundary conditions in Z
		if( space_->getBCLocal()->getBCL(Z)<=0 || SS_[Z]>space_->sIndB(EField::S,Z) ) {

			Ordinal k = SS_[Z];

			Ordinal SX = SS_[X];
			if( space_->getBCLocal()->getBCL(X)>0 && SS_[X]==space_->sIndB(EField::S,X) )
				++SX;
			Ordinal NX = NN_[X];
			if( space_->getBCLocal()->getBCU(X)>0 && NN_[X]==space_->eIndB(EField::S,X) )
				--NX;
			Ordinal SY = SS_[Y];
			if( space_->getBCLocal()->getBCL(Y)>0 && SS_[Y]==space_->sIndB(EField::S,Y) )
				++SY;
			Ordinal NY = NN_[Y];
			if( space_->getBCLocal()->getBCU(Y)>0 && NN_[Y]==space_->eIndB(EField::S,Y) )
				--NY;

			for( Ordinal j=SY; j<=NY; ++j )
				for( Ordinal i=SX; i<=NX; ++i )
					(*b)( getI(i,j,k) ) -= x.at(i,j,k-1)*op->getC( Z, k, -1 );

		}

		if( space_->getBCLocal()->getBCU(Z)<=0 || NN_[Z]<space_->eIndB(EField::S,Z) ) {

			Ordinal k = NN_[Z];

			Ordinal SX = SS_[X];
			if( space_->getBCLocal()->getBCL(X)>0 && SS_[X]==space_->sIndB(EField::S,X) )
				++SX;
			Ordinal NX = NN_[X];
			if( space_->getBCLocal()->getBCU(X)>0 && NN_[X]==space_->eIndB(EField::S,X) )
				--NX;
			Ordinal SY = SS_[Y];
			if( space_->getBCLocal()->getBCL(Y)>0 && SS_[Y]==space_->sIndB(EField::S,Y) )
				++SY;
			Ordinal NY = NN_[Y];
			if( space_->getBCLocal()->getBCU(Y)>0 && NN_[Y]==space_->eIndB(EField::S,Y) )
				--NY;

			for( Ordinal j=SY; j<=NY; ++j )
				for( Ordinal i=SX; i<=NX; ++i )
					(*b)( getI(i,j,k) ) -= x.at(i,j,k+1)*op->getC( Z, k, +1 );

		}

	}

	const std::string getLabel() const { return( "TeuchosTransfer" ); };

	void print( std::ostream& out=std::cout ) const {
		out << "--- " << getLabel() << " ---\n";
		out << "N: " << N_ << "\n";
		out << "cw: " << cw_ << "\n";
		out << "SS: (\t";
		for( int i=0; i<3; ++i ) out << SS_[i] << "\t";
		out << "\n";
		out << "NN: (\t";
		for( int i=0; i<3; ++i ) out << NN_[i] << "\t";
		out << "\n";


	}

}; // end of class TeuchosTransfer



} // end of namespace Pimpact



#ifdef COMPILE_ETI
extern template class Pimpact::TeuchosTransfer< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::TeuchosTransfer< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::TeuchosTransfer< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::TeuchosTransfer< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_TEUCHOSTRANSFER_HPP
