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

	using VectorT = Teuchos::SerialDenseVector<Ordinal,Scalar>;
	using MatrixT = Teuchos::SerialDenseMatrix<Ordinal,Scalar>;

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
			//SR_[i] = 1;
			//ER_[i] = space_->nLoc(i) - 1;
			SR_[i] = 1;
			ER_[i] = space_->nLoc(i);

			//if( space_->getBCLocal()->getBCL(i) >  0 ) SR_[i] = 2;
			//if( space_->getBCLocal()->getBCL(i) == 0 ) SR_[i] = 1;
			//if( space_->getBCLocal()->getBCU(i) >  0 ) ER_[i] = space_->nLoc(i) - 1;
			//if( space_->getBCLocal()->getBCU(i) == 0 ) ER_[i] = space_->nLoc(i);
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

		//// boundaries
		//for( int d=0; d<3; ++d ) {

			//int d1 = ( d + 1 )%3;
			//int d2 = ( d + 2 )%3;
			//if( d2>d1 ) std::swap( d2, d1 );
			//TO i;

			//// lower boundaries
			//if( space_->getBCLocal()->getBCL(d)>0 ) {
				//i[d] = 1;
				//for( i[d1]=getSR(d1); i[d1]<=getER(d1); ++i[d1] )
					//for( i[d2]=getSR(d2); i[d2]<=getER(d2); ++i[d2] ) {
						//TO ip = i;
						//++ip[d];
						//y.at(i[0],i[1],i[2]) =
							//getC(d,i[d],0)*x.at(i[0],i[1],i[2]) + getC(d,i[d],+1)*x.at(ip[0],ip[1],ip[2]);
					//}
			//}

			//// upper boundaries
			//if( space_->getBCLocal()->getBCU(d)>0 ) {
				//i[d] = space_->nLoc(d);
				//for( i[d1]=getSR(d1); i[d1]<=getER(d1); ++i[d1] )
					//for( i[d2]=getSR(d2); i[d2]<=getER(d2); ++i[d2] ) {
						//TO ip = i;
						//--ip[d];
						//y.at(i[0],i[1],i[2]) =
							//getC(d,i[d],0)*x.at(i[0],i[1],i[2]) + getC(d,i[d],-1)*x.at(ip[0],ip[1],ip[2]);
					//}
			//}

		//}

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
