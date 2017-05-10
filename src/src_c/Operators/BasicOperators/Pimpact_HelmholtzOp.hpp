#pragma once
#ifndef PIMPACT_HELMHOLTZOP_HPP
#define PIMPACT_HELMHOLTDOP_HPP


#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_Stencil.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{


/// \brief Helmholtz operator
/// \ingroup BaseOperator
///
/// computes \f$ y = ( mulI_ I - mulL_ \Delta) x \f$
template<class SpT>
class HelmholtzOp {

public:

  using SpaceT = SpT;

  using DomainFieldT = VectorField<SpaceT>;
  using RangeFieldT = VectorField<SpaceT>;

protected:

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

	static const int dimNC = SpT::dimNC;
	static const int dim = SpT::dimension;

	using SW = typename SpT::SW;

	using Stenc = Stencil< ST, OT, 0, SW::BL(0), SW::BU(0) >;

  using TO = const Teuchos::Tuple<Stenc,SpT::sdim>; 

  const Teuchos::RCP<const SpaceT> space_;

	ST mulI_;
	ST mulL_;

	const ST& mulBC_ = 1.0e+00; /// good idea, but has to be considered also in Force( what happens in Schurcomplement?
	
  TO cS_;
  TO cV_;

public:

	HelmholtzOp(
			const Teuchos::RCP<const SpaceT>& space ):
		space_( space ),
		mulI_( static_cast<ST>(0.) ),
		mulL_( 1./space_->getDomainSize()->getRe() ) {

			//const bool mapping = true; // order ~2
			const bool mapping = false; // order ~6

			for( int dir=0; dir<SpaceT::sdim; ++dir ) {

				F fdir = static_cast<F>( dir );

				// scalar stencil
				cS_[dir] = Stenc( space_->nLoc(dir) );
				FD_getDiffCoeff(
						1,
						space_->nLoc(dir),
						space_->bl(dir),
						space_->bu(dir),
						space_->bl(dir),
						space_->bu(dir),
						space_->getBCLocal()->getBCL(dir),
						space_->getBCLocal()->getBCU(dir),
						space_->getShift(dir),
						5,
						dir+1,
						2,
						0,
						mapping, // mapping
						space_->getStencilWidths()->getDimNcbC(dir),
						space_->getStencilWidths()->getNcbC(dir),
						space_->getCoordinatesLocal()->getX( F::S, dir ),
						space_->getCoordinatesLocal()->getX( F::S, dir ),
						cS_[dir].get() );

				// lower BC
				if( BC::Dirichlet==space_->bcl(dir) ) {
					for( int ii=SW::BL(X); ii<=SW::BU(X); ++ii )
						cS_[dir](1,ii) = 0.;
					cS_[dir](1,0) = mulBC_;
				}
				else if( BC::Neumann==space_->bcl(dir) ) {
					FD_getDiffCoeff(
							1,
							1,
							space_->bl(dir),
							space_->bu(dir),
							space_->bl(dir),
							space_->bu(dir),
							space_->getBCLocal()->getBCL(dir),
							0,
							space_->getShift(dir),
							5,
							dir+1,
							1,
							0,
							mapping, // mapping
							space_->getStencilWidths()->getDimNcbC(dir),
							space_->getStencilWidths()->getNcbC(dir),
							space_->getCoordinatesLocal()->getX( F::S, dir ),
							space_->getCoordinatesLocal()->getX( F::S, dir ),
							cS_[dir].get() );
				}
				// upper BC
				if( BC::Dirichlet==space_->bcu(dir) ) {
					for( int ii=SW::BL(X); ii<=SW::BU(X); ++ii )
						cS_[dir](space_->nLoc(dir),ii) = 0.;
					cS_[dir](space_->nLoc(dir),0) = mulBC_;
				}
				if( BC::Neumann==space_->bcu(dir) ) {
					FD_getDiffCoeff(
							space_->nLoc(dir),
							space_->nLoc(dir),
							space_->bl(dir),
							space_->bu(dir),
							space_->bl(dir),
							space_->bu(dir),
							0,
							space_->getBCLocal()->getBCU(dir),
							space_->getShift(dir),
							5,
							dir+1,
							1,
							0,
							mapping, // mapping
							space_->getStencilWidths()->getDimNcbC(dir),
							space_->getStencilWidths()->getNcbC(dir),
							space_->getCoordinatesLocal()->getX( F::S, dir ),
							space_->getCoordinatesLocal()->getX( F::S, dir ),
							cS_[dir].get() );
				}


				// velocity stencil
				cV_[dir] = Stenc( space_->nLoc(dir) );
				FD_getDiffCoeff(
						0,
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
						mapping,
						space_->getStencilWidths()->getDimNcbC(dir),
						space_->getStencilWidths()->getNcbC(dir),
						space_->getCoordinatesLocal()->getX( fdir, dir ),
						space_->getCoordinatesLocal()->getX( fdir, dir ),
						cV_[dir].get() );

				if( BC::Dirichlet==space_->bcl(dir) ) {
					for( OT ii=SW::BL(dir); ii<=SW::BU(dir); ++ii )
						cV_[dir](0,ii) = 0.;

					for( OT ii=SW::DL(dir); ii<=SW::DU(dir); ++ii )
						cV_[dir](0,ii+1) = space()->getInterpolateV2S()->getC( dir, 1, ii );
				}
				else if( BC::Neumann==space_->bcl(dir) ) {

					using StencD = Stencil< ST, OT, 0, SW::DL(0), SW::DU(0) >;

					StencD c_( space_->nLoc(dir) );

					FD_getDiffCoeff(
							1,
							space_->nLoc(dir),
							space_->bl(dir),
							space_->bu(dir),
							space_->dl(dir),
							space_->du(dir),
							space_->getBCLocal()->getBCL(dir),
							space_->getBCLocal()->getBCU(dir),
							space_->getShift(dir),
							3,
							dir+1,
							1,
							0,
							mapping, // mapping
							space_->getStencilWidths()->getDimNcbD(dir),
							space_->getStencilWidths()->getNcbD(dir),
							space_->getCoordinatesLocal()->getX( fdir, dir ),
							space_->getCoordinatesLocal()->getX( F::S, dir ),
							c_.get() );

					for( OT ii=SW::BL(dir); ii<=SW::BU(dir); ++ii )
						cV_[dir](0,ii) = 0.;

					for( OT ii=SW::DL(dir); ii<=SW::DU(dir); ++ii )
						cV_[dir](0,ii+1) = c_(1,ii);
				}


				if( BC::Dirichlet==space_->bcu(dir) ) {
					for( OT ii=SW::BL(dir); ii<=SW::BU(dir); ++ii )
						cV_[dir]( space()->ei(fdir,dir,B::Y), ii ) = 0.;

					for( OT ii=SW::DL(dir); ii<=SW::DU(dir); ++ii )
						cV_[dir](space()->ei(fdir,dir,B::Y), ii) =
							space()->getInterpolateV2S()->getC( dir, space()->ei(F::S,dir,B::Y), ii );
				}
				else if( BC::Neumann==space_->bcu(dir) ) {
					using StencD = Stencil< ST, OT, 0, SW::DL(0), SW::DU(0) >;

					StencD c_( space_->nLoc(dir) );

					FD_getDiffCoeff(
							1,
							space_->nLoc(dir),
							space_->bl(dir),
							space_->bu(dir),
							space_->dl(dir),
							space_->du(dir),
							space_->getBCLocal()->getBCL(dir),
							space_->getBCLocal()->getBCU(dir),
							space_->getShift(dir),
							3,
							dir+1,
							1,
							0,
							mapping, // mapping
							space_->getStencilWidths()->getDimNcbD(dir),
							space_->getStencilWidths()->getNcbD(dir),
							space_->getCoordinatesLocal()->getX( fdir, dir ),
							space_->getCoordinatesLocal()->getX( F::S, dir ),
							c_.get() );

					for( OT ii=SW::BL(dir); ii<=SW::BU(dir); ++ii )
						cV_[dir]( space()->ei(fdir,dir,B::Y), ii ) = 0.;

					for( OT ii=SW::DL(dir); ii<=SW::DU(dir); ++ii )
						cV_[dir](space()->ei(fdir,dir,B::Y), ii) = c_(space()->ei(F::S,dir,B::Y),ii);
				}
			}
		};



	/// \todo change with or without bc, everywhere
  void apply( const DomainFieldT& x, RangeFieldT& y, const Add& add=Add::N  ) const {

		const B& nb = B::N;

		for( int dir=0; dir<SpaceT::sdim; ++dir ) {

			F f = static_cast<F>(dir);

			x(f).exchange();

			if( 3==SpaceT::sdim ) {
				for( OT k=space()->si(f,Z,nb); k<=space()->ei(f,Z,nb); ++k )
					for( OT j=space()->si(f,Y,nb); j<=space()->ei(f,Y,nb); ++j )
						for( OT i=space()->si(f,X,nb); i<=space()->ei(f,X,nb); ++i ) {
							if( Add::N==add ) y(f)(i,j,k) = 0.;
							y(f)(i,j,k) +=
								mulI_*x(f)(i,j,k) - mulL_*innerStenc3D( x(f), f, i, j, k);
						}
			}
			else{
				for( OT k=space()->si(f,Z,nb); k<=space()->ei(f,Z,nb); ++k )
					for( OT j=space()->si(f,Y,nb); j<=space()->ei(f,Y,nb); ++j )
						for( OT i=space()->si(f,X,nb); i<=space()->ei(f,X,nb); ++i ) {
							if( Add::N==add ) y(f)(i,j,k) = 0.;
							y(f)(i,j,k) +=
								mulI_*x(f)(i,j,k) - mulL_*innerStenc2D( x(f), f, i, j, k);
						}
			}
		}
		if( Add::N==add ) applyBC( x, y );
    y.changed();
  }


	void applyBC( const VectorField<SpaceT>& x, VectorField<SpaceT>& y ) const {
		for( F field=F::U; field<SpaceT::sdim; ++field )
			applyBC( x(field), y(field) );
	}


	/// \brief implements Dirichlet boundary conditions as identity in tangential
	/// velocity direction or interpolation in wand normal direction
	void applyBC( const ScalarField<SpaceT>& x, ScalarField<SpaceT>& y	) const {

		assert( x.getType()==y.getType() );

		const F& f = x.getType();

		// U-field
		if( F::U==f ) {

			// tangential direction: Y
			if( 0<space()->bcl(Y) ) {
				OT j = space()->si(f,Y,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
						y(i,j,k) = 0.;
						for( OT jj=0; jj<=SW::BU(Y); ++jj )
							y(i,j,k) += getC(Y,f,j,jj)*x(i,j+jj,k);
					}
			}
			if( 0<space()->bcu(Y) ) {
				OT j = space()->ei(f,Y,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
						y(i,j,k) = 0.;
						for( OT jj=SW::BL(Y); jj<=0; ++jj )
							y(i,j,k) += getC(Y,f,j,jj)*x(i,j+jj,k);
					}
			}

			// tangential direction: Z
			if( 0<space()->bcl(Z) ) {
				OT k = space()->si(f,Z,B::Y);
				for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
					for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
						y(i,j,k) = 0.;
						for( OT kk=0; kk<=SW::BU(Z); ++kk )
							y(i,j,k) += getC(Z,f,k,kk)*x(i,j,k+kk);
					}
			}
			if( 0<space()->bcu(Z) ) {
				OT k = space()->ei(f,Z,B::Y);
				for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
					for( OT i=space()->si(f,X,B::N); i<=space()->ei(f,X,B::N); ++i ) {
						y(i,j,k) = 0.;
						for( OT kk=SW::BL(Z); kk<=0; ++kk )
							y(i,j,k) += getC(Z,f,k,kk)*x(i,j,k+kk);
					}
			}

			// normal direction: X
			if( 0<space()->bcl(X) ) {
				OT i = space()->si(f,X,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
						y(i,j,k) = 0.;
						for( OT ii=0; ii<=SW::BU(X); ++ii )
							y(i,j,k) += getC(X,f,i,ii)*x(i+ii,j,k);
					}
			}
			if( 0<space()->bcu(X) ) {
				OT i = space()->ei(f,X,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
						y(i,j,k) = 0.;
						for( OT ii=SW::BL(X); ii<=0; ++ii )
							y(i,j,k) += getC(X,f,i,ii)*x(i+ii,j,k);
					}
			}
		}

		// V-field
		if( F::V==f ) {

			// tangential direction: X
			if( 0<space()->bcl(X) ) {
				OT i = space()->si(f,X,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j ) {
						y(i,j,k) = 0.;
						for( OT ii=0; ii<=SW::BU(X); ++ii )
							y(i,j,k) += getC(X,f,i,ii)*x(i+ii,j,k);
					}
			}
			if( 0<space()->bcu(X) ) {
				OT i = space()->ei(f,X,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j ) {
						y(i,j,k) = 0.;
						for( OT ii=SW::BL(X); ii<=0; ++ii )
							y(i,j,k) += getC(X,f,i,ii)*x(i+ii,j,k);
					}
			}

			// tangential direction: Z
			if( 0<space()->bcl(Z) ) {
				OT k = space()->si(f,Z,B::Y);
				for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						y(i,j,k) = 0.;
						for( OT kk=0; kk<=SW::BU(Z); ++kk )
							y(i,j,k) += getC(Z,f,k,kk)*x(i,j,k+kk);
					}
			}
			if( 0<space()->bcu(Z) ) {
				OT k = space()->ei(f,Z,B::Y);
				for( OT j=space()->si(f,Y,B::N); j<=space()->ei(f,Y,B::N); ++j )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						y(i,j,k) = 0.;
						for( OT kk=SW::BL(Z); kk<=0; ++kk )
							y(i,j,k) += getC(Z,f,k,kk)*x(i,j,k+kk);
					}
			}

			// normal direction: Y
			if( 0<space()->bcl(Y) ) {
				OT j = space()->si(f,Y,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						y(i,j,k) = 0.;
						for( OT jj=0; jj<=SW::BU(Y); ++jj )
							y(i,j,k) += getC(Y,f,j,jj)*x(i,j+jj,k);
					}
			}
			if( 0<space()->bcu(Y) ) {
				OT j = space()->ei(f,Y,B::Y);
				for( OT k=space()->si(f,Z,B::Y); k<=space()->ei(f,Z,B::Y); ++k )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						y(i,j,k) = 0.;
						for( OT jj=SW::BL(Y); jj<=0; ++jj )
							y(i,j,k) += getC(Y,f,j,jj)*x(i,j+jj,k);
					}
			}
			
		}

		// W-field
		if( F::W==f ) {

			// tangential direction: X
			if( 0<space()->bcl(X) ) {
				OT i = space()->si(f,X,B::Y);
				for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
					for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
						y(i,j,k) = 0.;
						for( OT ii=0; ii<=SW::BU(X); ++ii )
							y(i,j,k) += getC(X,f,i,ii)*x(i+ii,j,k);
					}
			}
			if( 0<space()->bcu(X) ) {
				OT i = space()->ei(f,X,B::Y);
				for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
					for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j ) {
						y(i,j,k) = 0.;
						for( OT ii=SW::BL(X); ii<=0; ++ii )
							y(i,j,k) += getC(X,f,i,ii)*x(i+ii,j,k);
					}
			}

			// tangential direction: Y
			if( 0<space()->bcl(Y) ) {
				OT j = space()->si(f,Y,B::Y);
				for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						y(i,j,k) = 0.;
						for( OT jj=0; jj<=SW::BU(Y); ++jj )
							y(i,j,k) += getC(Y,f,j,jj)*x(i,j+jj,k);
					}
			}
			if( 0<space()->bcu(Y) ) {
				OT j = space()->ei(f,Y,B::Y);
				for( OT k=space()->si(f,Z,B::N); k<=space()->ei(f,Z,B::N); ++k )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						y(i,j,k) = 0.;
						for( OT jj=SW::BL(Y); jj<=0; ++jj )
							y(i,j,k) += getC(Y,f,j,jj)*x(i,j+jj,k);
					}
			}

			// normal direction: Z
			if( 0<space()->bcl(Z) ) {
				OT k = space()->si(f,Z,B::Y);
				for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						y(i,j,k) = 0.;
						for( OT kk=0; kk<=SW::BU(Z); ++kk )
							y(i,j,k) += getC(Z,f,k,kk)*x(i,j,k+kk);
					}
			}
			if( 0<space()->bcu(Z) ) {
				OT k = space()->ei(f,Z,B::Y);
				for( OT j=space()->si(f,Y,B::Y); j<=space()->ei(f,Y,B::Y); ++j )
					for( OT i=space()->si(f,X,B::Y); i<=space()->ei(f,X,B::Y); ++i ) {
						y(i,j,k) = 0.;
						for( OT kk=SW::BL(Z); kk<=0; ++kk )
							y(i,j,k) += getC(Z,f,k,kk)*x(i,j,k+kk);
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

			out << "\ncoord: " << static_cast<ECoord>(dir) << "\n";

			cS_[dir].print( out );
    }
    out << " --- velocity stencil: ---";

    for( int dir=0; dir<SpaceT::sdim; ++dir ) {

			out << "\ncoord: " << static_cast<ECoord>(dir) << "\n";

			cV_[dir].print( out );
    }
  }


	constexpr const Teuchos::RCP<const SpaceT>& space() { return( space_ ); };

  constexpr const ST* getC( const int& dir, const F& ftype ) {
		return( (dir==ftype)?cV_[dir].get():cS_[dir].get() );
  }

	constexpr const ST& getC( const ECoord& dir, const F& ftype, OT i, OT ii ) {
		return(
				(static_cast<int>(dir)==static_cast<int>(ftype))?
					cV_[dir](i,ii):
					cS_[dir](i,ii) );
	}

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		if( para->name()!="Linear Solver" ) {
			mulI_ = para->get<ST>( "mulI", 0. );
			mulL_ = para->get<ST>( "mulL", 1./space_->getDomainSize()->getRe() );
		}
	}

	const std::string getLabel() const { return( "Helmholtz" ); };


	constexpr ST innerStenc3D( const ScalarField<SpaceT>& x, const F& fType,
			const OT& i, const OT& j, const OT& k ) const {

		ST lap = 0.;

		for( int ii=SW::BL(X); ii<=SW::BU(X); ++ii ) 
			lap += getC(X,fType,i,ii)*x(i+ii,j,k);

		for( int jj=SW::BL(Y); jj<=SW::BU(Y); ++jj ) 
			lap += getC(Y,fType,j,jj)*x(i,j+jj,k);

		for( int kk=SW::BL(Z); kk<=SW::BU(Z); ++kk ) 
			lap += getC(Z,fType,k,kk)*x(i,j,k+kk);

		return( lap );
	}

	constexpr ST innerStenc2D( const ScalarField<SpaceT>& x, const F& fType,
			const OT& i, const OT& j, const OT& k ) const {

		ST lap = 0.;

		for( int ii=SW::BL(X); ii<=SW::BU(X); ++ii ) 
			lap += getC(X,fType,i,ii)*x(i+ii,j,k);

		for( int jj=SW::BL(Y); jj<=SW::BU(Y); ++jj ) 
			lap += getC(Y,fType,j,jj)*x(i,j+jj,k);

		return( lap );
	}

	constexpr ST innerDiag3D( const F& fType,
			const OT& i, const OT& j, const OT& k ) const {

		return( getC(X,fType,i,0) + getC(Y,fType,j,0) + getC(Z,fType,k,0) );

	}

	constexpr ST innerDiag2D( const F& fType,
			const OT& i, const OT& j, const OT& k ) const {

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
