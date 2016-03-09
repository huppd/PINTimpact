#pragma once
#ifndef PIMPACT_DIVGRADO2OP_HPP
#define PIMPACT_DIVGRADO2OP_HPP


// for EV
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact{


extern "C" {

void Op_getCDG(
    const int& dimens,
    const int* const M,
    const int* const N,
    const int* const BL,
    const int* const BU,
    const int* const BCL,
    const int* const BCU,
    const double* const y1u,
    const double* const y2v,
    const double* const y3w,
    const double* const x1p,
    const double* const x2p,
    const double* const x3p,
    const double* const x1u,
    const double* const x2v,
    const double* const x3w,
    double* const cdg1,
    double* const cdg2,
    double* const cdg3 );

void OP_DivGradO2Op(
    const int& dimens,
    const int* const N,
    const int* const SR,
    const int* const ER,
    const int* const BL,
    const int* const BU,
    const int* const BCL,
    const int* const BCU,
    const double* const cdg1,
    const double* const cdg2,
    const double* const cdg3,
    const double* const phi,
          double* const Lap );

}



/// \brief "Laplace" for pressure 2nd Order.
///
/// independent of \c StencilWidths
/// \ingroup BaseOperator
template<class ST>
class DivGradO2Op {

public:

  using SpaceT = ST;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

  using TO = const Teuchos::Tuple<Ordinal,3>;

protected:

  using TS = const Teuchos::Tuple<Scalar*,3>;

  const Teuchos::RCP<const SpaceT> space_;

  TS c_;

	TO SR_;
	TO ER_;

public:


	DivGradO2Op( const Teuchos::RCP<const SpaceT>& space ): space_(space) {

		for( int i=0; i<3; ++i ) {
			// allocate stencil
			Ordinal nTemp = 3*( space_->nLoc(i) - 1 + 1 );
			c_[i] = new Scalar[ nTemp ];

			// inner field bounds
			SR_[i] = 1;
			ER_[i] = space_->nLoc(i) - 1;

			if( space_->getBCLocal()->getBCL(i) >  0 ) SR_[i] = 2;
			if( space_->getBCLocal()->getBCL(i) == 0 ) SR_[i] = 1;
			if( space_->getBCLocal()->getBCU(i) >  0 ) ER_[i] = space_->nLoc(i) - 1;
			if( space_->getBCLocal()->getBCU(i) == 0 ) ER_[i] = space_->nLoc(i);
		}

		Op_getCDG(
				space_->dim(),
				space_->nGlo(),
				space_->nLoc(),
				space_->bl(),
				space_->bu(),
				space_->getBCLocal()->getBCL(),
				space_->getBCLocal()->getBCU(),
				space_->getCoordinatesGlobal()->getX( ECoord::X, EField::U ),
				space_->getCoordinatesGlobal()->getX( ECoord::Y, EField::V ),
				space_->getCoordinatesGlobal()->getX( ECoord::Z, EField::W ),
				space_->getCoordinatesLocal()->getX( ECoord::X, EField::S ),
				space_->getCoordinatesLocal()->getX( ECoord::Y, EField::S ),
				space_->getCoordinatesLocal()->getX( ECoord::Z, EField::S ),
				space_->getCoordinatesLocal()->getX( ECoord::X, EField::U ),
				space_->getCoordinatesLocal()->getX( ECoord::Y, EField::V ),
				space_->getCoordinatesLocal()->getX( ECoord::Z, EField::W ),
				c_[0],
				c_[1],
				c_[2] );
	}


	void apply( const DomainFieldT& x, RangeFieldT& y, Belos::ETrans
			trans=Belos::NOTRANS ) const {

		x.exchange();

		// inner stencil
		if( 3==space()->dim() )
			for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
				for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
					for( Ordinal i=getSR(X); i<=getER(X); ++i ) {
						y.at(i,j,k) = innerStenc3D(x, i,j,k);
					}
		else
			for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
				for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
					for( Ordinal i=getSR(X); i<=getER(X); ++i ) {
						y.at(i,j,k) = innerStenc2D(x, i,j,k);
					}

		// boundaries
		for( int d=0; d<3; ++d ) {

			int d1 = ( d + 1 )%3;
			int d2 = ( d + 2 )%3;
			if( d2>d1 ) std::swap( d2, d1 );
			TO i;

			// lower boundaries
			if( space_->getBCLocal()->getBCL(d)>0 ) {
				i[d] = 1;
				for( i[d1]=getSR(d1); i[d1]<=getER(d1); ++i[d1] )
					for( i[d2]=getSR(d2); i[d2]<=getER(d2); ++i[d2] ) {
						TO ip = i;
						++ip[d];
						y.at(i[0],i[1],i[2]) =
							getC(d,i[d],0)*x.at(i[0],i[1],i[2]) + getC(d,i[d],+1)*x.at(ip[0],ip[1],ip[2]);
					}
			}

			// upper boundaries
			if( space_->getBCLocal()->getBCU(d)>0 ) {
				i[d] = space_->nLoc(d);
				for( i[d1]=getSR(d1); i[d1]<=getER(d1); ++i[d1] )
					for( i[d2]=getSR(d2); i[d2]<=getER(d2); ++i[d2] ) {
						TO ip = i;
						--ip[d];
						y.at(i[0],i[1],i[2]) =
							getC(d,i[d],0)*x.at(i[0],i[1],i[2]) + getC(d,i[d],-1)*x.at(ip[0],ip[1],ip[2]);
					}
			}

		}

		y.changed();

	}

	void computeResidual( const RangeFieldT& b, const DomainFieldT& x,
			RangeFieldT& res ) const {

		x.exchange();
		// inner stencil
		if( 3==space()->dim() )
			for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
				for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
					for( Ordinal i=getSR(X); i<=getER(X); ++i ) {
						res.at(i,j,k) = b.at(i,j,k) - innerStenc3D(x, i,j,k);
					}
		else
			for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
				for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
					for( Ordinal i=getSR(X); i<=getER(X); ++i ) {
						res.at(i,j,k) = b.at(i,j,k) - innerStenc2D(x, i,j,k);
					}

		// boundaries
		for( int d=0; d<3; ++d ) {

			int d1 = ( d + 1 )%3;
			int d2 = ( d + 2 )%3;
			if( d2>d1 ) std::swap( d2, d1 );
			TO i;

			// lower boundaries
			if( space_->getBCLocal()->getBCL(d)>0 ) {
				i[d] = 1;
				for( i[d1]=getSR(d1); i[d1]<=getER(d1); ++i[d1] )
					for( i[d2]=getSR(d2); i[d2]<=getER(d2); ++i[d2] ) {
						TO ip = i;
						++ip[d];
						res.at(i[0],i[1],i[2]) = b.at(i[0],i[1],i[2]) -
							getC(d,i[d],0 )*x.at(i[0], i[1], i[2] ) -
							getC(d,i[d],+1)*x.at(ip[0],ip[1],ip[2]);
					}
			}

			// upper boundaries
			if( space_->getBCLocal()->getBCU(d)>0 ) {
				i[d] = space_->nLoc(d);
				for( i[d1]=getSR(d1); i[d1]<=getER(d1); ++i[d1] )
					for( i[d2]=getSR(d2); i[d2]<=getER(d2); ++i[d2] ) {
						TO ip = i;
						--ip[d];
						res.at(i[0],i[1],i[2]) = b.at(i[0],i[1],i[2]) -
							getC(d,i[d],0 )*x.at(i[0], i[1], i[2] ) -
							getC(d,i[d],-1)*x.at(ip[0],ip[1],ip[2]);
					}
			}

		}

		res.changed();

	}

	void computeEV( const ECoord& dir, Scalar& evMax, Scalar& evMin ) const {

		using VectorT = Teuchos::SerialDenseVector<Ordinal,Scalar>;
		using MatrixT = Teuchos::SerialDenseMatrix<Ordinal,Scalar>;

		Ordinal n = space_->nLoc(dir)-1+1;

		Teuchos::RCP< MatrixT > A;
		//Teuchos::RCP< VectorT > X;
		//Teuchos::RCP< VectorT > B;

		A = Teuchos::rcp( new MatrixT( n, n, true ) );
		// construct A
		for( Ordinal i=1; i<=space_->nLoc(dir); ++i )
			for( int o=-1; o<=1; ++o ) {
				Ordinal I = i-1;
				if( (i+o)>=getSR(dir) && (i+o)<=getER(dir) ) 
					(*A)( I, I+o ) += getC( dir, i, o) ;
			}


		// lower boundaries
		if( space_->getBCLocal()->getBCL(dir)>0 ) {
			Ordinal i = 1;
			Ordinal I = i-1;
			// set row zero
			for( Ordinal l=0; l<A->numCols(); ++l ) 
				(*A)(I,l) = 0.;
			(*A)( I, I )   += getC( dir, i, 0) ;
			(*A)( I, I+1 ) += getC( dir, i,+1) ;

		}

		// upper boundaries
		if( space_->getBCLocal()->getBCU(dir)>0 ) {
			Ordinal i = space_->nLoc(dir);
			Ordinal I = i-1;
			// set row zero
			for( Ordinal l=0; l<A->numCols(); ++l ) 
				(*A)(I,l) = 0.;
			(*A)( I, I )   += getC( dir, i, 0) ;
			(*A)( I, I-1 ) += getC( dir, i,-1) ;
		}

		// compute EV, \todo move to TeuchosBla
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

		Teuchos::RCP<std::ostream> out = Pimpact::createOstream( "ev_"+toString(dir)+".txt" );
		for( Ordinal i=0; i<n; ++i ) {
			*out << (*evr)[i] << "\n";
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

	void computeEV( Scalar& evMax, Scalar& evMin ) const {

		evMax = 0.;
		evMin = 0.;
		
		for( int i=0; i<3; ++i ) {
			Scalar evMaxTemp, evMinTemp;
			computeEV( static_cast<Pimpact::ECoord>(i), evMaxTemp, evMinTemp );
			evMax += evMaxTemp;
			evMin += evMinTemp;
		}

	}

	/// \name setter
	/// @{ 

  void assignField ( const DomainFieldT& mv ) const {};

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

	///  @} 

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    out << " --- stencil: ---";
		out << " sr: " << SR_ << "\n";
		out << " er: " << ER_ << "\n";
    for( int dir=0; dir<3; ++dir ) {
      out << "\ndir: " << dir << "\n";
      for( int i=1; i<=space_->nLoc(dir); ++i ) {
        out << "\ni: " << i << "\t(";
        for( int k=-1; k<=1; ++k ) {
					out << getC(dir,i,k) << "\t" ;
        }
        out << ")\n";
      }
      out << "\n";
    }
  }

  void print2Mat(  ) const {

    for( int dir=0; dir<3; ++dir ) {
			std::string fn = "A_" + toString( static_cast<ECoord>(dir) ) + "_" + std::to_string(space_->nLoc(dir)) + ".txt";

			Teuchos::RCP<std::ostream> out = Pimpact::createOstream( fn );
      for( int i=1; i<=space_->nLoc(dir); ++i ) {
        for( int k=-1; k<=1; ++k ) {
					*out << getC(dir,i,k) << "\t" ;
        }
        *out << "\n";
      }
    }
  }

protected:

	inline Scalar innerStenc3D( const DomainFieldT& x, const Ordinal& i, const Ordinal& j,
			const Ordinal& k ) const {

		return( 
				getC(X,i,-1)*x.at(i-1,j  ,k  ) + getC(X,i,1)*x.at(i+1,j  ,k  ) +
				getC(Y,j,-1)*x.at(i  ,j-1,k  ) + getC(Y,j,1)*x.at(i  ,j+1,k  ) +
				getC(Z,k,-1)*x.at(i  ,j  ,k-1) + getC(Z,k,1)*x.at(i  ,j  ,k+1) +
				( getC(X,i,0) + getC(Y,j,0) + getC(Z,k,0) )*x.at(i,j,k)
				);
	}

	inline Scalar innerStenc2D( const DomainFieldT& x, const Ordinal& i, const Ordinal& j,
			const Ordinal& k ) const {

		return( 
				getC(X,i,-1)*x.at(i-1,j  ,k  ) + getC(X,i,1)*x.at(i+1,j  ,k  ) +
				getC(Y,j,-1)*x.at(i  ,j-1,k  ) + getC(Y,j,1)*x.at(i  ,j+1,k  ) +
				( getC(X,i,0) + getC(Y,j,0) )*x.at(i,j,k)
				);
	}

public:

	/// \name getters
	/// @{ 

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(space_); };

	inline const Scalar* getC( const ECoord& dir) const  {
		return( getC( static_cast<const int&>(dir) ) );
  }

  inline const Scalar* getC( const int& dir) const  {
		return( c_[dir] );
  }

	inline const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const  {
		return( getC( static_cast<const int&>(dir), i, off ) );
  }

	inline const Scalar& getC( const int& dir, Ordinal i, Ordinal off ) const  {
		return( c_[dir][ off + 1 + (i-1)*3 ] );
  }

	inline const Ordinal* getSR() const { return( SR_.getRawPtr() ); }
	inline const Ordinal* getER() const { return( ER_.getRawPtr() ); }

	inline const Ordinal& getSR( const Ordinal& coord ) const { return( SR_[coord] ); }
	inline const Ordinal& getER( const Ordinal& coord ) const { return( ER_[coord] ); }

	inline const Ordinal& getSR( const ECoord& coord ) const { return( getSR( static_cast<Ordinal>(coord) ) ); }
	inline const Ordinal& getER( const ECoord& coord ) const { return( getER( static_cast<Ordinal>(coord) ) ); }

	const std::string getLabel() const { return( "DivGradO2" ); };

	///  @} 


}; // end of class DivGradO2Op





} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_DIVGRADO2OP_HPP
