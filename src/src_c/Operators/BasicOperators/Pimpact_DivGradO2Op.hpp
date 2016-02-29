#pragma once
#ifndef PIMPACT_DIVGRADO2OP_HPP
#define PIMPACT_DIVGRADO2OP_HPP


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



/// \brief "laplace" for pressure 2nd Order.
///
/// independent of \c StencilWidths
/// \ingroup BaseOperator
/// \todo instead of hardcode 2nd Order it would be pretty to use new space with \c StencilWidths<3,2>
/// \todo handle corner
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


	DivGradO2Op( const Teuchos::RCP<const SpaceT>& space ):
		space_(space) {

			for( int i=0; i<3; ++i ) {
				// alocate stencilt
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


	void apply(const DomainFieldT& x, RangeFieldT& y,
			Belos::ETrans trans=Belos::NOTRANS ) const {

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

		// boundary conditions in X
		if( space_->getBCLocal()->getBCL(X)>0 ) {

			Ordinal i = 1;

			for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
				for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
					y.at(i,j,k) = getC(X,i,0)*x.at(i,j,k) + getC(X,i,+1)*x.at(i+1,j,k);

		}
		if( space_->getBCLocal()->getBCU(X)>0 ) {

			Ordinal i = space_->nLoc(X);

			for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
				for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
					y.at(i,j,k) = getC(X,i,-1)*x.at(i-1,j,k) + getC(X,i,0)*x.at(i,j,k);

		}

		// boundary conditions in Y
		if( space_->getBCLocal()->getBCL(Y)>0 ) {

			Ordinal j = 1;

			for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
				for( Ordinal i=getSR(X); i<=getER(X); ++i )
					y.at(i,j,k) = getC(Y,j,0)*x.at(i,j,k) + getC(Y,j,+1)*x.at(i,j+1,k);

		}
		if( space_->getBCLocal()->getBCU(Y)>0 ) {

			Ordinal j = space_->nLoc(Y);

			for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
				for( Ordinal i=getSR(X); i<=getER(X); ++i )
					y.at(i,j,k) = getC(Y,j,-1)*x.at(i,j-1,k) + getC(Y,j,0)*x.at(i,j,k);

		}

		// boundary conditions in Z
		if( space_->getBCLocal()->getBCL(Z)>0 ) {

			Ordinal k = 1;

			for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
				for( Ordinal i=getSR(X); i<=getER(X); ++i )
					y.at(i,j,k) = getC(Z,k,0)*x.at(i,j,k) + getC(Z,k,+1)*x.at(i,j,k+1);

		}
		if( space_->getBCLocal()->getBCU(Z)>0 ) {

			Ordinal k = space_->nLoc(Z);

			for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
				for( Ordinal i=getSR(X); i<=getER(X); ++i )
					y.at(i,j,k) = getC(Z,k,-1)*x.at(i,j,k-1) + getC(Z,k,0)*x.at(i,j,k);

		}

		//y.setCornersZero(); // ???

		y.changed();

	}

	void computeResidual( const RangeFieldT& b, const DomainFieldT& x, RangeFieldT& res ) const {

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

		// boundary conditions in X
		if( space_->getBCLocal()->getBCL(X)>0 ) {

			Ordinal i = 1;

			for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
				for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
					res.at(i,j,k) = b.at(i,j,k) - getC(X,i,0)*x.at(i,j,k) - getC(X,i,+1)*x.at(i+1,j,k);

		}
		if( space_->getBCLocal()->getBCU(X)>0 ) {

			Ordinal i = space_->nLoc(X);

			for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
				for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
					res.at(i,j,k) = b.at(i,j,k) - getC(X,i,-1)*x.at(i-1,j,k) - getC(X,i,0)*x.at(i,j,k);

		}

		// boundarres.conditions in Y
		if( space_->getBCLocal()->getBCL(Y)>0 ) {

			Ordinal j = 1;

			for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
				for( Ordinal i=getSR(X); i<=getER(X); ++i )
					res.at(i,j,k) = b.at(i,j,k) - getC(Y,j,0)*x.at(i,j,k) - getC(Y,j,+1)*x.at(i,j+1,k);

		}
		if( space_->getBCLocal()->getBCU(Y)>0 ) {

			Ordinal j = space_->nLoc(Y);

			for( Ordinal k=getSR(Z); k<=getER(Z); ++k )
				for( Ordinal i=getSR(X); i<=getER(X); ++i )
					res.at(i,j,k) = b.at(i,j,k) - getC(Y,j,-1)*x.at(i,j-1,k) - getC(Y,j,0)*x.at(i,j,k);

		}

		// boundary conditions in Z
		if( space_->getBCLocal()->getBCL(Z)>0 ) {

			Ordinal k = 1;

			for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
				for( Ordinal i=getSR(X); i<=getER(X); ++i )
					res.at(i,j,k) = b.at(i,j,k) - getC(Z,k,0)*x.at(i,j,k) - getC(Z,k,+1)*x.at(i,j,k+1);

		}
		if( space_->getBCLocal()->getBCU(Z)>0 ) {

			Ordinal k = space_->nLoc(Z);

			for( Ordinal j=getSR(Y); j<=getER(Y); ++j )
				for( Ordinal i=getSR(X); i<=getER(X); ++i )
					res.at(i,j,k) = b.at(i,j,k) - getC(Z,k,-1)*x.at(i,j,k-1) - getC(Z,k,0)*x.at(i,j,k);

		}

		//res.setCornersZero(); // ???

		res.changed();
		//apply( x, res );
		//res.add( 1., b, -1., res );
	}

  void assignField ( const DomainFieldT& mv ) const {};

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(space_); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

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

	Scalar innerStenc3D( const DomainFieldT& x, const Ordinal& i, const Ordinal& j,
			const Ordinal& k ) const {

		return( 
				getC(X,i,-1)*x.at(i-1,j  ,k  ) + getC(X,i,1)*x.at(i+1,j  ,k  ) +
				getC(Y,j,-1)*x.at(i  ,j-1,k  ) + getC(Y,j,1)*x.at(i  ,j+1,k  ) +
				getC(Z,k,-1)*x.at(i  ,j  ,k-1) + getC(Z,k,1)*x.at(i  ,j  ,k+1) +
				( getC(X,i,0) + getC(Y,j,0) + getC(Z,k,0) )*x.at(i,j,k)
				);
	}

	Scalar innerStenc2D( const DomainFieldT& x, const Ordinal& i, const Ordinal& j,
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

	const Scalar* getC( const ECoord& dir) const  {
		return( getC( static_cast<const int&>(dir) ) );
  }

  const Scalar* getC( const int& dir) const  {
		return( c_[dir] );
  }

	const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const  {
		return( getC( static_cast<const int&>(dir), i, off ) );
  }

	const Scalar& getC( const int& dir, Ordinal i, Ordinal off ) const  {
		return( c_[dir][ off + 1 + (i-1)*3 ] );
  }

	const Ordinal* getSR() const { return( SR_.getRawPtr() ); }
	const Ordinal* getER() const { return( ER_.getRawPtr() ); }

	const Ordinal& getSR( const Ordinal& coord ) const { return( SR_[coord] ); }
	const Ordinal& getER( const Ordinal& coord ) const { return( ER_[coord] ); }

	const Ordinal& getSR( const ECoord& coord ) const { return( getSR( static_cast<Ordinal>(coord) ) ); }
	const Ordinal& getER( const ECoord& coord ) const { return( getER( static_cast<Ordinal>(coord) ) ); }

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
