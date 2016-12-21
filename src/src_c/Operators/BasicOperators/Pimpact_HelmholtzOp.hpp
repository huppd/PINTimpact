#pragma once
#ifndef PIMPACT_HELMHOLTZOP_HPP
#define PIMPACT_HELMHOLTDOP_HPP


#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_Stencil.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{


/// \brief Helmholtz operator
/// \ingroup BaseOperator
///
/// computes \f$ y = ( mulI_ I - mulL_ \Delta) x \f$
template<class ST>
class HelmholtzOp {

public:

  using SpaceT = ST;

  using DomainFieldT = VectorField<SpaceT>;
  using RangeFieldT = VectorField<SpaceT>;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

	static const int dimNC = ST::dimNC;
	static const int dim = ST::dimension;

	using SW = StencilWidths<dim,dimNC>;

	using Stenc = Stencil< Scalar, Ordinal, 0, SW::BL(0), SW::BU(0) >;

  using TO = const Teuchos::Tuple<Stenc,ST::sdim>; 

  const Teuchos::RCP<const SpaceT> space_;

	Scalar mulI_;
	Scalar mulL_;
	
  TO cS_;
  TO cV_;

public:

	/// \todo change stencil at B for Neumann BC
	HelmholtzOp(
			const Teuchos::RCP<const SpaceT>& space ):
		space_( space ),
		mulI_( static_cast<Scalar>(0.) ),
		mulL_( 1./space_->getDomainSize()->getRe() ) {

    for( int dir=0; dir<SpaceT::sdim; ++dir ) {
			cS_[dir] = Stenc( space_->nLoc(dir) );
			FD_getDiffCoeff(
					space_->nLoc(dir),
					space_->bl(dir),
					space_->bu(dir),
					space_->bl(dir),
					space_->bu(dir),
					space_->getBCLocal()->getBCL(dir),
					space_->getBCLocal()->getBCU(dir),
					space_->getShift(dir),
					int(EField::S)+1,
					dir+1,
					2,
					0,
					//true,
					false, // mapping
					space_->getStencilWidths()->getDimNcbC(dir),
					space_->getStencilWidths()->getNcbC(dir),
					space_->getCoordinatesLocal()->getX( dir, EField::S ),
					space_->getCoordinatesLocal()->getX( dir, EField::S ),
					cS_[dir].get() );

			//if( DirichletBC==space_->bcl(dir) ) {
				//Ordinal i = space_->begin(S,dir,With::B);
				//for( int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii )
					//cS_[dir](i,ii) = 0.;
				//cS_[dir](i,0) = -1./mulL_;
			//}
			//if( DirichletBC==space_->bcu(dir) ) {
				//Ordinal i = space_->end(S,dir,With::B);
				//for( int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii )
					//cS_[dir](i,ii) = 0.;
				//cS_[dir](i,0) = -1./mulL_;
			//}

		
			cV_[dir] = Stenc( space_->nLoc(dir) );
			FD_getDiffCoeff(
					space_->nLoc(dir),
					space_->bl(dir),
					space_->bu(dir),
					space_->bl(dir),
					space_->bu(dir),
					space_->getBCLocal()->getBCL(dir),
					space_->getBCLocal()->getBCU(dir),
					space_->getShift(dir),
					1,
					dir+1,
					2,
					0,
					true,
					//false,
					space_->getStencilWidths()->getDimNcbC(dir),
					space_->getStencilWidths()->getNcbC(dir),
					space_->getCoordinatesLocal()->getX( dir, dir ),
					space_->getCoordinatesLocal()->getX( dir, dir ),
					cV_[dir].get() );

			//if( DirichletBC==space_->bcl(dir) ) {
				//Ordinal i = space_->begin(dir,dir,With::B);
				//for( int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii )
					//cV_[dir](i,ii) = 0.;
				//for( int ii=SW::DL(0); ii<=SW::DU(0); ++ii )
					//cV_[dir](i,ii) = -space_->getInterpolateV2S()->getC(dir,i+1,ii)/mulL_;
			//}
			//if( DirichletBC==space_->bcu(dir) ) {
				//Ordinal i = space_->end(dir,dir,With::B);
				//for( int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii )
					//cV_[dir](i,ii) = 0.;
				//for( int ii=SW::DL(0); ii<=SW::DU(0); ++ii )
					//cV_[dir](i,ii) = -space_->getInterpolateV2S()->getC(dir,i,ii)/mulL_;
			//}
    }
  };



	/// \todo change with or without bc, everywhere
  void apply( const DomainFieldT& x, RangeFieldT& y, const Add& add=Add::No  ) const {

		const With& wB = ( (Add::No==add) ? With::B : With::noB );
		const With& wnB = With::noB;

		for( int dir=0; dir<SpaceT::sdim; ++dir ) {

			EField fType = static_cast<EField>(dir);

			x.getField(fType).exchange();

			if( 3==SpaceT::sdim ) {
				for( Ordinal k=space()->begin(fType,Z,wnB); k<=space()->end(fType,Z,wnB); ++k )
					for( Ordinal j=space()->begin(fType,Y,wnB); j<=space()->end(fType,Y,wnB); ++j )
						for( Ordinal i=space()->begin(fType,X,wnB); i<=space()->end(fType,X,wnB); ++i ) {
							if( Add::No==add ) y.getField(fType)(i,j,k) = 0.;
							y.getField(fType)(i,j,k) +=
								mulI_*x.getConstField(fType)(i,j,k) - mulL_*innerStenc3D( x.getConstField(fType), fType, i, j, k);
						}
			}
			else{
				for( Ordinal k=space()->begin(fType,Z,wnB); k<=space()->end(fType,Z,wnB); ++k )
					for( Ordinal j=space()->begin(fType,Y,wnB); j<=space()->end(fType,Y,wnB); ++j )
						for( Ordinal i=space()->begin(fType,X,wnB); i<=space()->end(fType,X,wnB); ++i ) {
							if( Add::No==add ) y.getField(fType)(i,j,k) = 0.;
							y.getField(fType)(i,j,k) +=
								mulI_*x.getConstField(fType)(i,j,k) - mulL_*innerStenc2D( x.getConstField(fType), fType, i, j, k);
						}
			}
			if( wB==With::B ) applyBC( x.getConstField(fType), y.getField(fType) );
		}
    y.changed();
  }


	/// \brief implements Dirichlet boundary conditions as identity in tangential
	/// velocity direction or interpolation in wand normal direction
	void applyBC( const ScalarField<SpaceT>& x, ScalarField<SpaceT>& y	) const {

		assert( x.getType()==y.getType() );

		// U-field
		if( U==x.getType() ) {
			if( DirichletBC==space()->bcl(Y) ) {
				Ordinal j = space()->begin(U,Y,With::B);
				for( Ordinal k=space()->begin(U,Z,With::B); k<=space()->end(U,Z,With::B); ++k )
					for( Ordinal i=space()->begin(U,X,With::B); i<=space()->end(U,X,With::B); ++i ) {
						y(i,j,k) = x(i,j,k);
					}
			}
			if( DirichletBC==space()->bcu(Y) ) {
				Ordinal j = space()->end(U,Y,With::B);
				for( Ordinal k=space()->begin(U,Z,With::B); k<=space()->end(U,Z,With::B); ++k )
					for( Ordinal i=space()->begin(U,X,With::B); i<=space()->end(U,X,With::B); ++i ) {
						y(i,j,k) = x(i,j,k);
					}
			}

			if( DirichletBC==space()->bcl(Z) ) {
				Ordinal k = space()->begin(U,Z,With::B);
				for( Ordinal j=space()->begin(U,Y,With::B); j<=space()->end(U,Y,With::B); ++j )
					for( Ordinal i=space()->begin(U,X,With::B); i<=space()->end(U,X,With::B); ++i ) {
						y(i,j,k) = x(i,j,k);
					}
			}
			if( DirichletBC==space()->bcl(Z) ) {
				Ordinal k = space()->end(U,Z,With::B);
				for( Ordinal j=space()->begin(U,Y,With::B); j<=space()->end(U,Y,With::B); ++j )
					for( Ordinal i=space()->begin(U,X,With::B); i<=space()->end(U,X,With::B); ++i ) {
						y(i,j,k) = x(i,j,k);
					}
			}

			if( DirichletBC==space()->bcl(X) ) {
				Ordinal i = space()->begin(U,X,With::B);
				for( Ordinal k=space()->begin(U,Z,With::B); k<=space()->end(U,Z,With::B); ++k )
					for( Ordinal j=space()->begin(U,Y,With::B); j<=space()->end(U,Y,With::B); ++j ) {
						y(i,j,k) = 0.;
						for( Ordinal ii=space()->dl(X); ii<=space()->du(X); ++ii )
							y(i,j,k) += space()->getInterpolateV2S()->getC( X, i+1, ii )*x(1+i+ii,j,k);
					}
			}
			if( DirichletBC==space()->bcu(X) ) {
				Ordinal i = space()->end(U,X,With::B);
				for( Ordinal k=space()->begin(U,Z,With::B); k<=space()->end(U,Z,With::B); ++k )
					for( Ordinal j=space()->begin(U,Y,With::B); j<=space()->end(U,Y,With::B); ++j ) {
						y(i,j,k) = 0.;
						for( Ordinal ii=space()->dl(X); ii<=space()->du(X); ++ii )
							y(i,j,k) += space()->getInterpolateV2S()->getC( X, i, ii )*x(i+ii,j,k);
					}
			}
		}

		// V-field
		if( V==x.getType() ) {
			if( DirichletBC==space()->bcl(X) ) {
				Ordinal i = space()->begin(V,X,With::B);
				for( Ordinal k=space()->begin(V,Z,With::noB); k<=space()->end(V,Z,With::noB); ++k )
					for( Ordinal j=space()->begin(V,Y,With::B); j<=space()->end(V,Y,With::B); ++j )
						y(i,j,k) = x(i,j,k);
			}
			if( DirichletBC==space()->bcu(X) ) {
				Ordinal i = space()->end(V,X,With::B);
				for( Ordinal k=space()->begin(V,Z,With::noB); k<=space()->end(V,Z,With::noB); ++k )
					for( Ordinal j=space()->begin(V,Y,With::B); j<=space()->end(V,Y,With::B); ++j )
						y(i,j,k) = x(i,j,k);
			}

			if( DirichletBC==space()->bcl(Z) ) {
				Ordinal k = space()->begin(V,Z,With::B);
				for( Ordinal j=space()->begin(V,Y,With::noB); j<=space()->end(V,Y,With::noB); ++j )
					for( Ordinal i=space()->begin(V,X,With::B); i<=space()->end(V,X,With::B); ++i ) {
						y(i,j,k) = x(i,j,k);
					}
			}
			if( DirichletBC==space()->bcu(Z) ) {
				Ordinal k = space()->end(V,Z,With::B);
				for( Ordinal j=space()->begin(V,Y,With::noB); j<=space()->end(V,Y,With::noB); ++j )
					for( Ordinal i=space()->begin(V,X,With::B); i<=space()->end(V,X,With::B); ++i ) {
						y(i,j,k) = x(i,j,k);
					}
			}

			if( DirichletBC==space()->bcl(Y) ) {
				Ordinal j = space()->begin(V,Y,With::B);
				for( Ordinal k=space()->begin(V,Z,With::B); k<=space()->end(V,Z,With::B); ++k )
					for( Ordinal i=space()->begin(V,X,With::B); i<=space()->end(V,X,With::B); ++i ) {
						y(i,j,k) = 0.;
						for( Ordinal jj=space()->dl(Y); jj<=space()->du(Y); ++jj )
							y(i,j,k) += space()->getInterpolateV2S()->getC( Y, j+1, jj )*x(i,1+j+jj,k);
					}
			}
			if( DirichletBC==space()->bcu(Y) ) {
				Ordinal j = space()->end(V,Y,With::B);
				for( Ordinal k=space()->begin(V,Z,With::B); k<=space()->end(V,Z,With::B); ++k )
					for( Ordinal i=space()->begin(V,X,With::B); i<=space()->end(V,X,With::B); ++i ) {
						y(i,j,k) = 0.;
						for( Ordinal jj=space()->dl(Y); jj<=space()->du(Y); ++jj )
							y(i,j,k) += space()->getInterpolateV2S()->getC( Y, j, jj )*x(i,j+jj,k);
					}
			}
		}

		// W-field
		if( W==x.getType() ) {
			if( DirichletBC==space()->bcl(X) ) {
				Ordinal i = space()->begin(W,X,With::B);
				for( Ordinal k=space()->begin(W,Z,With::noB); k<=space()->end(W,Z,With::noB); ++k )
					for( Ordinal j=space()->begin(W,Y,With::B); j<=space()->end(W,Y,With::B); ++j )
						y(i,j,k) = x(i,j,k);
			}
			if( DirichletBC==space()->bcu(X) ) {
				Ordinal i = space()->end(W,X,With::B);
				for( Ordinal k=space()->begin(W,Z,With::noB); k<=space()->end(W,Z,With::noB); ++k )
					for( Ordinal j=space()->begin(W,Y,With::B); j<=space()->end(W,Y,With::B); ++j )
						y(i,j,k) = x(i,j,k);
			}

			if( DirichletBC==space()->bcl(Y) ) {
				Ordinal j = space()->begin(W,Y,With::B);
				for( Ordinal k=space()->begin(W,Z,With::noB); k<=space()->end(W,Z,With::noB); ++k )
					for( Ordinal i=space()->begin(W,X,With::B); i<=space()->end(W,X,With::B); ++i ) {
						y(i,j,k) = x(i,j,k);
					}
			}
			if( DirichletBC==space()->bcu(Y) ) {
				Ordinal j = space()->end(W,Y,With::B);
				for( Ordinal k=space()->begin(W,Z,With::noB); k<=space()->end(W,Z,With::noB); ++k )
					for( Ordinal i=space()->begin(W,X,With::B); i<=space()->end(W,X,With::B); ++i ) {
						y(i,j,k) = x(i,j,k);
					}
			}

			if( DirichletBC==space()->bcl(Z) ) {
				Ordinal k = space()->begin(W,Z,With::B);
				for( Ordinal j=space()->begin(W,Y,With::B); j<=space()->end(W,Y,With::B); ++j )
					for( Ordinal i=space()->begin(W,X,With::B); i<=space()->end(W,X,With::B); ++i ) {
						y(i,j,k) = 0.;
						for( Ordinal kk=space()->dl(Z); kk<=space()->du(Z); ++kk )
							y(i,j,k) += space()->getInterpolateV2S()->getC( Z, k+1, kk )*x(i,j,1+k+kk);
					}
			}
			if( DirichletBC==space()->bcu(Z) ) {
				Ordinal k = space()->end(W,Z,With::B);
				for( Ordinal j=space()->begin(W,Y,With::B); j<=space()->end(W,Y,With::B); ++j )
					for( Ordinal i=space()->begin(W,X,With::B); i<=space()->end(W,X,With::B); ++i ) {
						y(i,j,k) = 0.;
						for( Ordinal kk=space()->dl(Z); kk<=space()->du(Z); ++kk )
							y(i,j,k) += space()->getInterpolateV2S()->getC( Z, k, kk )*x(i,j,k+kk);
					}
			}
		}
	}


  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

  void print( std::ostream& out=std::cout ) const {

    out << "--- " << getLabel() << " ---\n";
    out << " --- scalar stencil: ---";

    for( int dir=0; dir<SpaceT::sdim; ++dir ) {

			out << "\ncoord: " << toString( static_cast<ECoord>(dir) ) << "\n";

			cS_[dir].print( out );
    }
    out << " --- velocity stencil: ---";

    for( int dir=0; dir<SpaceT::sdim; ++dir ) {

			out << "\ncoord: " << toString( static_cast<ECoord>(dir) ) << "\n";

			cV_[dir].print( out );
    }
  }


	constexpr const Teuchos::RCP<const SpaceT>& space() const { return( space_ ); };

  constexpr const Scalar* getC( const int& dir, const int& ftype ) const {
		return( (dir==ftype)?cV_[dir].get():cS_[dir].get() );
  }

	constexpr const Scalar& getC( const ECoord& dir, const EField& ftype, Ordinal i, Ordinal ii ) const {
		return(
				(static_cast<int>(dir)==static_cast<int>(ftype))?
					cV_[dir](i,ii):
					cS_[dir](i,ii) );
	}

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		mulI_ = para->get<Scalar>( "mulI", 0. );
		mulL_ = para->get<Scalar>( "mulL", 1./space_->getDomainSize()->getRe() );
	}

	const std::string getLabel() const { return( "Helmholtz" ); };


	inline constexpr Scalar innerStenc3D( const ScalarField<SpaceT>& x, const EField& fType,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar lap = 0.;

		for( int ii=space_->bl(X); ii<=space_->bu(X); ++ii ) 
			lap += getC(X,fType,i,ii)*x(i+ii,j,k);

		for( int jj=space_->bl(Y); jj<=space_->bu(Y); ++jj ) 
			lap += getC(Y,fType,j,jj)*x(i,j+jj,k);

		for( int kk=space_->bl(Z); kk<=space_->bu(Z); ++kk ) 
			lap += getC(Z,fType,k,kk)*x(i,j,k+kk);

		return( lap );
	}

	inline constexpr Scalar innerStenc2D( const ScalarField<SpaceT>& x, const EField& fType,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar lap = 0.;

		for( int ii=space_->bl(X); ii<=space_->bu(X); ++ii ) 
			lap += getC(X,fType,i,ii)*x(i+ii,j,k);

		for( int jj=space_->bl(Y); jj<=space_->bu(Y); ++jj ) 
			lap += getC(Y,fType,j,jj)*x(i,j+jj,k);

		return( lap );
	}

	inline constexpr Scalar innerDiag3D( const EField& fType,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		return( getC(X,fType,i,0) + getC(Y,fType,j,0) + getC(Z,fType,k,0) );

	}

	inline constexpr Scalar innerDiag2D( const EField& fType,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		return( getC(X,fType,i,0) + getC(Y,fType,j,0) );

	}

}; // end of class HelmholtzOp



} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::HelmholtzOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::HelmholtzOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::HelmholtzOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::HelmholtzOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_HELMHOLTZOP_HPP
