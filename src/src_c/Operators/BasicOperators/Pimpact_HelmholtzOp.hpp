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
					true,
					//false, // mapping
					space_->getStencilWidths()->getDimNcbC(dir),
					space_->getStencilWidths()->getNcbC(dir),
					space_->getCoordinatesLocal()->getX( dir, EField::S ),
					space_->getCoordinatesLocal()->getX( dir, EField::S ),
					cS_[dir].get() );

			if( DirichletBC==space_->bcl(dir) ) {
				Ordinal i = space_->begin(S,dir,With::B);
				for( int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii )
					cS_[dir].at(i,ii) = 0.;
				cS_[dir].at(i,0) = 1.;
			}
			if( DirichletBC==space_->bcu(dir) ) {
				Ordinal i = space_->end(S,dir,With::B);
				for( int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii )
					cS_[dir].at(i,ii) = 0.;
				cS_[dir].at(i,0) = 1.;
			}

		
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

			if( DirichletBC==space_->bcl(dir) ) {
				Ordinal i = space_->begin(dir,dir,With::B);
				for( int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii )
					cV_[dir].at(i,ii) = 0.;
				for( int ii=SW::DL(0); ii<=SW::DU(0); ++ii )
					cV_[dir].at(i,ii) = space_->getInterpolateV2S()->getC(dir,i+1,ii);
			}
			if( DirichletBC==space_->bcu(dir) ) {
				Ordinal i = space_->end(dir,dir,With::B);
				for( int ii=Stenc::bl(); ii<=Stenc::bu(); ++ii )
					cV_[dir].at(i,ii) = 0.;
				for( int ii=SW::DL(0); ii<=SW::DU(0); ++ii )
					cV_[dir].at(i,ii) = space_->getInterpolateV2S()->getC(dir,i,ii);
			}
    }
  };



	/// \todo change with or without bc, everywhere
  void apply(const DomainFieldT& x, RangeFieldT& y ) const {

		for( int dir=0; dir<SpaceT::sdim; ++dir ) {

			EField fType = static_cast<EField>(dir);

			x.getField(fType).exchange();

			if( 3==SpaceT::sdim ) {
				for( Ordinal k=space()->begin(fType,Z); k<=space()->end(fType,Z); ++k )
					for( Ordinal j=space()->begin(fType,Y); j<=space()->end(fType,Y); ++j )
						for( Ordinal i=space()->begin(fType,X); i<=space()->end(fType,X); ++i )
							y.getField(fType).at(i,j,k) =
								mulI_*x.getConstField(fType).at(i,j,k) - mulL_*innerStenc3D( x.getConstField(fType), fType, i, j, k);
			}
			else{
				for( Ordinal k=space()->begin(fType,Z); k<=space()->end(fType,Z); ++k )
					for( Ordinal j=space()->begin(fType,Y); j<=space()->end(fType,Y); ++j )
						for( Ordinal i=space()->begin(fType,X); i<=space()->end(fType,X); ++i )
							y.getField(fType).at(i,j,k) =
								mulI_*x.getConstField(fType).at(i,j,k) - mulL_*innerStenc2D( x.getConstField(fType), fType, i, j, k);
			}
		}

    y.changed();
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
					cV_[dir].at(i,ii):
					cS_[dir].at(i,ii) );
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
			lap += getC(X,fType,i,ii)*x.at(i+ii,j,k);

		for( int jj=space_->bl(Y); jj<=space_->bu(Y); ++jj ) 
			lap += getC(Y,fType,j,jj)*x.at(i,j+jj,k);

		for( int kk=space_->bl(Z); kk<=space_->bu(Z); ++kk ) 
			lap += getC(Z,fType,k,kk)*x.at(i,j,k+kk);

		return( lap );
	}

	inline constexpr Scalar innerStenc2D( const ScalarField<SpaceT>& x, const EField& fType,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar lap = 0.;

		for( int ii=space_->bl(X); ii<=space_->bu(X); ++ii ) 
			lap += getC(X,fType,i,ii)*x.at(i+ii,j,k);

		for( int jj=space_->bl(Y); jj<=space_->bu(Y); ++jj ) 
			lap += getC(Y,fType,j,jj)*x.at(i,j+jj,k);

		return( lap );
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
