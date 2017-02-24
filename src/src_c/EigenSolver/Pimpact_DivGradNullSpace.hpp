#pragma once
#ifndef PIMPACT_DIVGRADNULLSPACE_HPP
#define PIMPACT_DIVGRADNULLSPACE_HPP


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Pimpact_Stencil.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Space.hpp" // just for createOstream<>




namespace Pimpact {



/// \todo make nullspace greate again, (store vectors) not hole field
template<class OperatorT>
class DivGradNullSpace {

protected:

	using SpaceT = typename OperatorT::SpaceT;

	using Scalar  = typename SpaceT::Scalar;
	using Ordinal = typename SpaceT::Ordinal;

	static const int dimNC = SpaceT::dimNC;
	static const int dim = SpaceT::dimension;

	using SW = StencilWidths<dim,dimNC>;

	using Stenc = Stencil< Scalar, Ordinal, SW::BL(0), SW::BL(0), SW::BU(0) >;

	using DomainFieldT= typename OperatorT::DomainFieldT;
	using RangeFieldT = typename OperatorT::RangeFieldT;

	using VectorT = Teuchos::SerialDenseVector<Ordinal,Scalar>;
	using MatrixT = Teuchos::SerialDenseMatrix<Ordinal,Scalar>;

	using SolverT = Teuchos::SerialQRDenseSolver<Ordinal,Scalar>;

public:

	/// \todo think about changing BC solver to proper one ( inner field, setting last component one and move it to the rhs)
	void computeNullSpace( const Teuchos::RCP<const OperatorT>& div, RangeFieldT& y, const bool& DJG_yes = true )  {

		Teuchos::RCP<const SpaceT> space = div->space();

		Teuchos::Tuple< Teuchos::RCP<VectorT>, SpaceT::sdim > x_ ;

		for( int dir=0; dir<SpaceT::sdim; ++dir ) {
			const Ordinal& N = space->nGlo(dir);
			x_[dir] = Teuchos::rcp( new VectorT(N) );

			if( -1==space->getBCGlobal()->getBCL(dir) ) {
				*x_[dir] = 1.;
			}
			else {
				// global stencil
				Ordinal nTempG = ( space->nGlo(dir) + space->bu(dir) - space->bl(dir) + 1 )
					*( space->bu(dir) - space->bl(dir) + 1);

				Stenc cG1( space->nGlo(dir) + space->bu(dir) );
				Stenc cG2( space->nGlo(dir) + space->bu(dir) );

				for( Ordinal i = space->begin(F::S,dir); i<=space->end(F::S,dir); ++i )
					for( Ordinal ii = space->dl(dir); ii<=space->du(dir); ++ii )
						cG1( i+space->getShift(dir), ii )= div->getC( static_cast<ECoord>(dir), i, ii );

				/// \todo Allgather would be better
				MPI_Allreduce(
						cG1.get(),    		                          // const void *sendbuf,
						cG2.get(),    		                          // void *recvbuf,
						nTempG,			                                // int count,
						MPI_REAL8,	                                // MPI_Datatype datatype,
						MPI_SUM,		                                // MPI_Op op,
						space->getProcGrid()->getCommSlice(dir) );	// MPI_Comm comm )


				// generate global Div Stencil(copy from transposed stencil)
				// generate Matrix
				Teuchos::RCP<MatrixT> D = Teuchos::rcp( new MatrixT( N + 1, N, true ) );
				for( Ordinal i=0; i<=N; ++i ) {
					for( int ii=space->gl(dir); ii<=space->gu(dir); ++ii ) {
						if( 0<=i+ii-1 && i+ii-1<N )
							(*D)(i,i+ii-1) = cG2(i+ii,-ii);
					}
				}
				// generate RHS
				Teuchos::RCP<VectorT> b = Teuchos::rcp( new VectorT( N + 1, true ) );

				if( space->getBCGlobal()->getBCL(dir)>0 ) {
					if( DJG_yes ) {
						(*b)(0) = -space->getInterpolateV2S()->getC( static_cast<Pimpact::ECoord>(dir),1,-1);
						(*b)(1) = -space->getInterpolateV2S()->getC( static_cast<Pimpact::ECoord>(dir),1, 0);
						(*b)(2) = -space->getInterpolateV2S()->getC( static_cast<Pimpact::ECoord>(dir),1, 1);
						(*b)(3) = -space->getInterpolateV2S()->getC( static_cast<Pimpact::ECoord>(dir),1, 2);
						MPI_Bcast(
								&(*b)(0),																		//void* data,
								4,																					//int count,
								MPI_DOUBLE,																	//MPI_Datatype datatype,
								0,																					//int root,
								space->getProcGrid()->getCommSlice(dir) );	// MPI_Comm comm )
					}
					else
						(*b)(0) = -1.;
				}

				if( space->getBCGlobal()->getBCU(dir)>0 ) {
					if( DJG_yes ) {
					(*b)(N-3) = space->getInterpolateV2S()->getC( static_cast<Pimpact::ECoord>(dir),space->end(F::S,dir),-3);
					(*b)(N-2) = space->getInterpolateV2S()->getC( static_cast<Pimpact::ECoord>(dir),space->end(F::S,dir),-2);
					(*b)(N-1) = space->getInterpolateV2S()->getC( static_cast<Pimpact::ECoord>(dir),space->end(F::S,dir),-1);
					(*b)(N  ) = space->getInterpolateV2S()->getC( static_cast<Pimpact::ECoord>(dir),space->end(F::S,dir), 0);
					MPI_Bcast(
							&(*b)(N-3),																	//void* data,
							4,																					//int count,
							MPI_DOUBLE,																	//MPI_Datatype datatype,
							space->getProcGrid()->getNP(dir)-1,					//int root,
							space->getProcGrid()->getCommSlice(dir) );	// MPI_Comm comm )
					}
					else
						(*b)(N) = 1.;
				}

				Teuchos::RCP<SolverT> solv_ = Teuchos::rcp( new SolverT() );
				solv_->factorWithEquilibration( true );
				solv_->setMatrix( D );
				solv_->factor();
				solv_->setVectors( x_[dir], b );
				solv_->solve();
				// solve with QRU
				if( 0==space->rankST() ) {
					std::cout << "dir: " << dir << "\n";
					std::cout << std::scientific;
					std::cout << std::setprecision( std::numeric_limits<long double>::digits10 + 1);
					std::cout << *x_[dir] << "\n";

					Teuchos::RCP<std::ostream> output = Pimpact::createOstream( "null.txt" );
					x_[dir]->print( *output );
				}
			}
		}
		if( 3==SpaceT::sdim ) {
			for( Ordinal k=space()->begin(F::S,Z); k<=space()->end(F::S,Z); ++k )
				for( Ordinal j=space()->begin(F::S,Y); j<=space()->end(F::S,Y); ++j )
					for( Ordinal i=space()->begin(F::S,X); i<=space()->end(F::S,X); ++i ) {
							y(i,j,k) = (*x_[0])(i-1+space->getShift(0))* (*x_[1])(j-1+space->getShift(1))* (*x_[2])(k-1+space->getShift(2));
					}
			//y(1,1,1) = 0.;
			//y(space()->end(F::S,X),1,1) = 0.;
			//y(1.,space()->end(F::S,Y),1) = 0.;
			//y(1.,1.,space()->end(F::S,Z)) = 0.;
			//y(1.,space()->end(F::S,Y),space()->end(F::S,Z)) = 0.;
			//y(space()->end(F::S,X),1.,space()->end(F::S,Z)) = 0.;
			//y(space()->end(F::S,X),space()->end(F::S,Y),1) = 0.;
			//y(space()->end(F::S,X),space()->end(F::S,Y),space()->end(F::S,Z)) = 0.;
		}
		else{
			for( Ordinal k=space()->begin(F::S,Z); k<=space()->end(F::S,Z); ++k )
				for( Ordinal j=space()->begin(F::S,Y); j<=space()->end(F::S,Y); ++j )
					for( Ordinal i=space()->begin(F::S,X); i<=space()->end(F::S,X); ++i )
						y(i,j,k) = (*x_[0])(i-1+space->getShift(0))* (*x_[1])(j-1+space->getShift(1));
		}
		Scalar blup = std::sqrt( 1./y.dot(y) );
		y.scale( blup );
		//y.write( SpaceT::sdim );
	}


}; // end of class SimpleVectorIteration


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DIVGRADNULLSPACE_HPP
