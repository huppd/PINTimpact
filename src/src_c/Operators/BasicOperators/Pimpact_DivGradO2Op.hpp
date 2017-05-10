#pragma once
#ifndef PIMPACT_DIVGRADO2OP_HPP
#define PIMPACT_DIVGRADO2OP_HPP


// for EV
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact{



/// \brief "Laplace" for pressure 2nd Order.
///
/// independent of \c StencilWidths
/// \ingroup BaseOperator
template<class SpT>
class DivGradO2Op {

public:

  using SpaceT = SpT;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

	static constexpr int epsI = 1e2;

protected:

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

	using Stenc = Stencil< ST, OT, 1, -1, 1 >;
  using TS = const Teuchos::Tuple< Stenc, SpT::sdim >;

  const Teuchos::RCP<const SpaceT> space_;

  TS c_;

public:

	DivGradO2Op( const Teuchos::RCP<const SpaceT>& space ): space_(space) {

		using StencG = Stencil< ST, OT, 0,  0, 1 >;
		using StencD = Stencil< ST, OT, 1, -1, 0 >;

		const auto& x = *space_->getCoordinatesLocal();
		const auto& y = *space_->getCoordinatesGlobal();

		for( int dir=0; dir<SpT::sdim; ++dir ) {
			F f = static_cast<F>( dir );
			// allocate stencils
			StencG cG( space_->nLoc(dir) );
			StencD cD( space_->nLoc(dir) );

			c_[dir] = Stenc( space_->nLoc(dir) );

			for( OT i=0; i<=space_->nLoc(dir); ++i ) {
				ST h = 1./( x(F::S,dir,i+1) - x(F::S,dir,i) );
				cG(i,0) = - h;
				cG(i,1) =   h;
			}

			for( OT i=1; i<=space_->nLoc(dir); ++i ) {
				ST h = 1./( x(f,dir,i) - x(f,dir,i-1) );
				cD(i,-1) = - h;
				cD(i, 0) =   h;
			}

			if( 0<space_->bcl(dir) ) {
				cG(0,0) = 0.;
				cG(0,1) = 0.;
				cD(1,-1) = 0.;
				//cD(1, 0) = 0.;
				cD(1, 0) = 1./( y(f,dir,1) - y(f,dir,0) );
				//cD(1, 0) = 1./( x(f,dir,1) - x(f,dir,0) );
				//cD(1, 0) = 1./( x(F::S,dir,1) - x(F::S,dir,0) );
			}
			if( 0<space_->bcu(dir) ) {
				cG( space_->nLoc(dir), 0 ) = 0.;
				cG( space_->nLoc(dir), 1 ) = 0.;
				cD( space_->nLoc(dir),-1 ) = 0.;
				cD( space_->nLoc(dir), 0 ) = 0.;
				cD( space_->nLoc(dir),-1 ) = -1./( y(f,dir,space_->nGlo(dir)) - y(f,dir,space_->nGlo(dir)-1) );
				//cD( space_->nLoc(dir),-1 ) = -1./( x(f,dir,space_->nLoc(dir)) - x(f,dir,space_->nLoc(dir)-1) );
				//cD( space_->nLoc(dir),-1 ) = -1./( x(F::S,dir,space_->nLoc(dir)) - x(F::S,dir,space_->nLoc(dir)-1) );
			}

			for( OT i=1; i<=space_->nLoc(dir); ++i ) {
				c_[dir](i,-1) = cD(i,-1)*cG(i-1,0);
				c_[dir](i, 0) = cD(i,-1)*cG(i-1,1) + cD(i,0)*cG(i,0);
				c_[dir](i, 1) =                      cD(i,0)*cG(i,1);
			}

			if( 0<space_->bcl(dir) ) {
				c_[dir](1,-1)  = 0.;
				c_[dir](1, 0) *= 2.;
				c_[dir](1, 1) *= 2.;
			}
			if( 0<space_->bcu(dir) ) {
				c_[dir](space_->nLoc(dir),-1) *= 2.;
				c_[dir](space_->nLoc(dir), 0) *= 2.;
				c_[dir](space_->nLoc(dir), 1)  = 0.;
			}

			if( BC::Neumann==space_->bcl(dir) ) {
				c_[dir](1,-1) = 0.;
				c_[dir](1, 0) = 0.;
				c_[dir](1, 1) = 0.;
			}
			if( BC::Neumann==space_->bcu(dir) ) {
				c_[dir](space_->nLoc(dir),-1) = 0.;
				c_[dir](space_->nLoc(dir), 0) = 0.;
				c_[dir](space_->nLoc(dir), 1) = 0.;
			}

			if( BC::Symmetry==space_->bcl(dir) ) {
				c_[dir](1, 1) += c_[dir](1,-1);
				c_[dir](1,-1)  = 0.;
			}
			if( BC::Symmetry==space_->bcu(dir) ) {
				c_[dir](space_->nLoc(dir),-1) += c_[dir](space_->nLoc(dir),1);
				c_[dir](space_->nLoc(dir), 1)  = 0.;
			}
		}
	}


	void apply( const DomainFieldT& x, RangeFieldT& y, const Add& add=Add::N ) const {

		x.exchange();

		if( 3==SpaceT::sdim ) {
			for( OT k=space()->si(F::S,Z); k<=space()->ei(F::S,Z); ++k )
				for( OT j=space()->si(F::S,Y); j<=space()->ei(F::S,Y); ++j )
					for( OT i=space()->si(F::S,X); i<=space()->ei(F::S,X); ++i ) {
						if( Add::N==add ) y(i,j,k) = 0.;
						y(i,j,k) += innerStenc3D(x, i,j,k);
					}
		}
		else {
			for( OT k=space()->si(F::S,Z); k<=space()->ei(F::S,Z); ++k )
				for( OT j=space()->si(F::S,Y); j<=space()->ei(F::S,Y); ++j )
					for( OT i=space()->si(F::S,X); i<=space()->ei(F::S,X); ++i ) {
						if( Add::N==add ) y(i,j,k) = 0.;
						y(i,j,k) += innerStenc2D(x, i,j,k);
					}
		}
	}


	void computeResidual( const RangeFieldT& b, const DomainFieldT& x,
			RangeFieldT& res ) const {

		x.exchange();
		// inner stencil
		if( 3==SpaceT::sdim ) {
			for( OT k=space()->si(F::S,Z); k<=space()->ei(F::S,Z); ++k )
				for( OT j=space()->si(F::S,Y); j<=space()->ei(F::S,Y); ++j )
					for( OT i=space()->si(F::S,X); i<=space()->ei(F::S,X); ++i )
						res(i,j,k) = b(i,j,k) - innerStenc3D(x, i,j,k);
		}
		else {
			for( OT k=space()->si(F::S,Z); k<=space()->ei(F::S,Z); ++k )
				for( OT j=space()->si(F::S,Y); j<=space()->ei(F::S,Y); ++j )
					for( OT i=space()->si(F::S,X); i<=space()->ei(F::S,X); ++i )
						res(i,j,k) = b(i,j,k) - innerStenc2D(x, i,j,k);
		}

		res.changed();
	}


	void applyInvDiag( const DomainFieldT& x, RangeFieldT& y ) const {

		//const ST& eps = 0.1;
		const ST& eps = 1./static_cast<ST>(epsI);

		if( 3==SpaceT::sdim ) {
			for( OT k=space()->si(F::S,Z); k<=space()->ei(F::S,Z); ++k )
				for( OT j=space()->si(F::S,Y); j<=space()->ei(F::S,Y); ++j )
					for( OT i=space()->si(F::S,X); i<=space()->ei(F::S,X); ++i ) {

						const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->si(F::S,X) ) ||
							(               space()->getBCLocal()->getBCU(X) > 0 && i==space()->ei(F::S,X) ) ;
						const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->si(F::S,Y) ) ||
							(               space()->getBCLocal()->getBCU(Y) > 0 && j==space()->ei(F::S,Y) ) ;
						const bool bcZ = (space()->getBCLocal()->getBCL(Z) > 0 && k==space()->si(F::S,Z) ) ||
							(               space()->getBCLocal()->getBCU(Z) > 0 && k==space()->ei(F::S,Z) ) ;

						const ST epsX = ( (bcY||bcZ)?eps:1. );
						const ST epsY = ( (bcX||bcZ)?eps:1. );
						const ST epsZ = ( (bcX||bcY)?eps:1. );

						ST diag = std::fabs( epsX*getC(X,i,0) + epsY*getC(Y,j,0) + epsZ*getC(Z,k,0) );
						assert( diag!=0. );
						y(i,j,k) = x(i,j,k)/diag;
					}
		}
		else {
			for( OT k=space()->si(F::S,Z); k<=space()->ei(F::S,Z); ++k )
				for( OT j=space()->si(F::S,Y); j<=space()->ei(F::S,Y); ++j )
					for( OT i=space()->si(F::S,X); i<=space()->ei(F::S,X); ++i ) {
						const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->si(F::S,X) ) ||
							(               space()->getBCLocal()->getBCU(X) > 0 && i==space()->ei(F::S,X) ) ;
						const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->si(F::S,Y) ) ||
							(               space()->getBCLocal()->getBCU(Y) > 0 && j==space()->ei(F::S,Y) ) ;

						const ST epsX = ( bcY?eps:1. );
						const ST epsY = ( bcX?eps:1. );

						y(i,j,k) = x(i,j,k)/std::fabs( epsX*getC(X,i,0) + epsY*getC(Y,j,0) );
					}
		}

		y.changed();
	}

	/// \name setter
	/// @{ 

  void assignField ( const DomainFieldT& mv ) const {};

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

	///  @} 

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    out << " --- stencil: ---";
    for( int dir=0; dir<SpT::sdim; ++dir ) {
      out << "\ndir: " << dir << "\n";
			c_[dir].print( out );
    }
  }

  void print2Mat(  ) const {

    for( int dir=0; dir<SpT::sdim; ++dir ) {
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


	constexpr ST innerStenc3D( const DomainFieldT& x, const OT& i, const OT& j,
			const OT& k ) const {

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->si(F::S,X) ) ||
			(               space()->getBCLocal()->getBCU(X) > 0 && i==space()->ei(F::S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->si(F::S,Y) ) ||
			(               space()->getBCLocal()->getBCU(Y) > 0 && j==space()->ei(F::S,Y) ) ;
		const bool bcZ = (space()->getBCLocal()->getBCL(Z) > 0 && k==space()->si(F::S,Z) ) ||
			(               space()->getBCLocal()->getBCU(Z) > 0 && k==space()->ei(F::S,Z) ) ;

		const ST& eps = 1./static_cast<ST>(epsI);

		const ST epsX = ( (bcY||bcZ)?eps:1. );
		const ST epsY = ( (bcX||bcZ)?eps:1. );
		const ST epsZ = ( (bcX||bcY)?eps:1. );

		return( 
				epsX*getC(X,i,-1)*x(i-1,j  ,k  ) + epsX*getC(X,i,1)*x(i+1,j  ,k  ) +
				epsY*getC(Y,j,-1)*x(i  ,j-1,k  ) + epsY*getC(Y,j,1)*x(i  ,j+1,k  ) +
				epsZ*getC(Z,k,-1)*x(i  ,j  ,k-1) + epsZ*getC(Z,k,1)*x(i  ,j  ,k+1) +
				( epsX*getC(X,i,0) + epsY*getC(Y,j,0) + epsZ*getC(Z,k,0) )*x(i,j,k)
				);
	}

	constexpr ST innerStenc2D( const DomainFieldT& x, const OT& i, const OT& j,
			const OT& k ) const {

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->si(F::S,X) ) ||
		           (space()->getBCLocal()->getBCU(X) > 0 && i==space()->ei(F::S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->si(F::S,Y) ) ||
		           (space()->getBCLocal()->getBCU(Y) > 0 && j==space()->ei(F::S,Y) ) ;

		const ST& eps = 1./static_cast<ST>(epsI);

		const ST epsX = (bcY)?eps:1.;
		const ST epsY = (bcX)?eps:1.;

		return( 
				epsX*getC(X,i,-1)*x(i-1,j  ,k  ) + epsX*getC(X,i,1)*x(i+1,j  ,k  ) +
				epsY*getC(Y,j,-1)*x(i  ,j-1,k  ) + epsY*getC(Y,j,1)*x(i  ,j+1,k  ) +
				( epsX*getC(X,i,0) + epsY*getC(Y,j,0) )*x(i,j,k)
				);
	}

	constexpr ST innerDiag3D( const OT& i, const OT& j, const OT& k ) {

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->si(F::S,X) ) ||
			(               space()->getBCLocal()->getBCU(X) > 0 && i==space()->ei(F::S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->si(F::S,Y) ) ||
			(               space()->getBCLocal()->getBCU(Y) > 0 && j==space()->ei(F::S,Y) ) ;
		const bool bcZ = (space()->getBCLocal()->getBCL(Z) > 0 && k==space()->si(F::S,Z) ) ||
			(               space()->getBCLocal()->getBCU(Z) > 0 && k==space()->ei(F::S,Z) ) ;

		const ST& eps = 1./static_cast<ST>(epsI);

		const ST epsX = ( (bcY||bcZ)?eps:1. );
		const ST epsY = ( (bcX||bcZ)?eps:1. );
		const ST epsZ = ( (bcX||bcY)?eps:1. );

		return( epsX*getC(X,i,0) + epsY*getC(Y,j,0) + epsZ*getC(Z,k,0) );
	}

	constexpr ST innerDiag2D( const OT& i, const OT& j,
			const OT& k ) {

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->si(F::S,X) ) ||
			(               space()->getBCLocal()->getBCU(X) > 0 && i==space()->ei(F::S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->si(F::S,Y) ) ||
			(               space()->getBCLocal()->getBCU(Y) > 0 && j==space()->ei(F::S,Y) ) ;

		const ST& eps = 1./static_cast<ST>(epsI);

		const ST epsX = ( bcY?eps:1. );
		const ST epsY = ( bcX?eps:1. );

		return( epsX*getC(X,i,0) + epsY*getC(Y,j,0) );
	}

	/// \name getters
	/// @{ 

  bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(space_); };

	constexpr const ST* getC( const ECoord& dir) const  {
		return( getC( static_cast<const int&>(dir) ) );
  }

  constexpr const ST* getC( const int& dir) const  {
		return( c_[dir].get() );
  }

	constexpr const ST& getC( const ECoord& dir, OT i, OT off ) const  {
		return( getC( static_cast<const int&>(dir), i, off ) );
  }

	constexpr const ST& getC( const int& dir, OT i, OT off ) const  {
		return( c_[dir](i,off) );
  }

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
