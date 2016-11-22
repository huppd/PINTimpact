#pragma once
#ifndef PIMPACT_CONVECTIONSOP_HPP
#define PIMPACT_CONVECTIONSOP_HPP


#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Stencil.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact {


/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup NonliearOperator
template<class ST>
class ConvectionSOp {

public:

  using SpaceT = ST;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using FluxFieldT = Teuchos::Tuple< Teuchos::RCP< ScalarField<SpaceT> >, 3 >;
  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

	static const int dimNC = ST::dimNC;
	static const int dim = ST::dimension;
	static const int sdim = ST::sdim;

	using SW = StencilWidths<dim,dimNC>;

	using Stenc = Stencil< Scalar, Ordinal, 0, SW::NL(0), SW::NU(0) >;
	using TO = const Teuchos::Tuple< Stenc*, sdim >;

  const Teuchos::RCP<const SpaceT> space_;

  TO cSD_;
  TO cVD_;

  TO cSU_;
  TO cVU_;

public:

	ConvectionSOp( const Teuchos::RCP<const SpaceT>& space  ):
		space_(space) {

			for( int i=0; i<sdim; ++i ) {

				cSD_[i] = new Stenc( space_->nLoc(i) );

				FD_getDiffCoeff(
						space_->nLoc(i),
						space_->bl(i),
						space_->bu(i),
						space_->nl(i),
						space_->nu(i),
						space_->getBCLocal()->getBCL(i),
						space_->getBCLocal()->getBCU(i),
						space_->getShift(i),
						int(EField::S)+1,
						i+1,
						1,
						-1,
						true, // mapping
						//false, // mapping
						space_->getStencilWidths()->getDimNcbC(i),
						space_->getStencilWidths()->getNcbC(i),
						space_->getCoordinatesLocal()->getX( i, EField::S ),
						space_->getCoordinatesLocal()->getX( i, EField::S ),
						cSD_[i]->get() );

				cSU_[i] = new Stenc( space_->nLoc(i) );

				FD_getDiffCoeff(
						space_->nLoc(i),
						space_->bl(i),
						space_->bu(i),
						space_->nl(i),
						space_->nu(i),
						space_->getBCLocal()->getBCL(i),
						space_->getBCLocal()->getBCU(i),
						space_->getShift(i),
						int(EField::S)+1,
						i+1,
						1,
						+1,
						true, // mapping
						//false, // mapping
						space_->getStencilWidths()->getDimNcbC(i),
						space_->getStencilWidths()->getNcbC(i),
						space_->getCoordinatesLocal()->getX( i, EField::S ),
						space_->getCoordinatesLocal()->getX( i, EField::S ),
						cSU_[i]->get() );

				cVD_[i] = new Stenc( space_->nLoc(i) );

				FD_getDiffCoeff(
						space_->nLoc(i),
						space_->bl(i),
						space_->bu(i),
						space_->nl(i),
						space_->nu(i),
						space_->getBCLocal()->getBCL(i),
						space_->getBCLocal()->getBCU(i),
						space_->getShift(i),
						1,
						i+1,
						1,
						-1,
						true, // mapping
						//false, // mapping
						space_->getStencilWidths()->getDimNcbC(i),
						space_->getStencilWidths()->getNcbC(i),
						space_->getCoordinatesLocal()->getX( i, i ),
						space_->getCoordinatesLocal()->getX( i, i ),
						cVD_[i]->get() );

				cVU_[i] = new Stenc( space_->nLoc(i) );

				FD_getDiffCoeff(
						space_->nLoc(i),
						space_->bl(i),
						space_->bu(i),
						space_->nl(i),
						space_->nu(i),
						space_->getBCLocal()->getBCL(i),
						space_->getBCLocal()->getBCU(i),
						space_->getShift(i),
						1,
						i+1,
						1,
						+1,
						true, // mapping
						//false, // mapping
						space_->getStencilWidths()->getDimNcbC(i),
						space_->getStencilWidths()->getNcbC(i),
						space_->getCoordinatesLocal()->getX( i, i ),
						space_->getCoordinatesLocal()->getX( i, i ),
						cVU_[i]->get() );
			}
		};


  void assignField( const RangeFieldT& mv ) {};


	void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z,
			Scalar mul, Scalar mulI, Scalar mulC, Scalar mulL ) const {
		apply( x, y, z, mul, mulC );
	}


	void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, Scalar mul=0., Scalar mulC=1. ) const {

		int m = static_cast<int>(z.getType());

		assert( z.getType() == y.getType() );
		for( int i=0; i<SpaceT::sdim; ++i ) 
			assert( x[i]->getType()==y.getType() );
#ifndef NDEBUG
		TEUCHOS_TEST_FOR_EXCEPTION(
				z.getType() != y.getType(),
				std::logic_error,
				"Pimpact::ConvectionSOP can only be applied to same fieldType !!!\n");


		for( int i=0; i<SpaceT::sdim; ++i ) {
			TEUCHOS_TEST_FOR_EXCEPTION(
					x[i]->getType() != y.getType(),
					std::logic_error,
					"Pimpact::ConvectionSOP can only be applied to same fieldType !!!\n");
		}
#endif

		for( int vel_dir=0; vel_dir<SpaceT::sdim; ++vel_dir )
			x[vel_dir]->exchange();

		y.exchange();

		if( 3==SpaceT::sdim )
			for( Ordinal k=space()->begin(m,Z); k<=space()->end(m,Z); ++k )
				for( Ordinal j=space()->begin(m,Y); j<=space()->end(m,Y); ++j )
					for( Ordinal i=space()->begin(m,X); i<=space()->end(m,X); ++i )
						z.at(i,j,k) = mul*z.at(i,j,k) + mulC*innerStenc3D( x[0]->at(i,j,k), x[1]->at(i,j,k), x[2]->at(i,j,k), y, i, j, k );
		else
			for( Ordinal k=space()->begin(m,Z); k<=space()->end(m,Z); ++k )
				for( Ordinal j=space()->begin(m,Y); j<=space()->end(m,Y); ++j )
					for( Ordinal i=space()->begin(m,X); i<=space()->end(m,X); ++i )
						z.at(i,j,k) = mul*z.at(i,j,k) + mulC*innerStenc2D( x[0]->at(i,j,k), x[1]->at(i,j,k), y, i,j,k);

		z.changed();
	}

	void print( std::ostream& out=std::cout ) const {
		out << " --- ConvectioSOp ---\n";
		for( int i=0; i<SpaceT::sdim; ++i ) {
			out << "dir: " << toString( static_cast<ECoord>(i) ) << "\n ";
			out << "cSD:\n";
			cSD_[i]->print( out );
			out << "cSU:\n";
			cSU_[i]->print( out );
			out << "cVD:\n";
			cVD_[i]->print( out );
			out << "cVU:\n";
			cVU_[i]->print( out );
		}
	}


  bool hasApplyTranspose() const { return( false ); }

  constexpr const Teuchos::RCP<const SpaceT>&  space() const { return( space_ ); }

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}


  constexpr const Scalar* getCU( const ECoord& dir, const EField& ftype ) const  {
    return( ( ((int)dir)==((int)ftype) )?cVU_[dir]->get():cSU_[dir]->get() );
  }

  constexpr const Scalar* getCD( const ECoord& dir, const EField& ftype ) const  {
    return( ( ((int)dir)==((int)ftype) )?cVD_[dir]->get():cSD_[dir]->get() );
  }

	inline constexpr const Scalar& getC( const Scalar& wind, const ECoord& dir, const EField& ftype, const int& i, const int ii ) const {
		return( 
				( static_cast<int>(dir)==static_cast<int>(ftype) )?
					(wind>=0? cVU_[dir]->at(i,ii):cVD_[dir]->at(i,ii))
					:
					(wind>=0? cSU_[dir]->at(i,ii):cSD_[dir]->at(i,ii))
				);
	}


	inline constexpr Scalar innerStenc2D( const Scalar& u, const Scalar& v, const
			RangeFieldT& x, const Ordinal& i, const Ordinal& j, const Ordinal& k )
		const {

		Scalar dx = 0.;
		for( int ii=space_->nl(X); ii<=space_->nu(X); ++ii ) 
			dx += getC( u,X, x.getType(),i,ii)*x.at(i+ii,j,k);

		Scalar dy = 0.;
		for( int jj=space_->nl(Y); jj<=space_->nu(Y); ++jj ) 
			dy += getC( v,Y, x.getType(),j,jj)*x.at(i,j+jj,k);

		return( u*dx+v*dy );
	}

	inline constexpr Scalar innerStenc3D( const Scalar& u, const Scalar& v,
			const Scalar& w, const RangeFieldT& x, const Ordinal& i, const Ordinal&
			j, const Ordinal& k ) const {

		Scalar dx = 0.;
		for( int ii=space_->nl(X); ii<=space_->nu(X); ++ii ) 
			dx += getC( u,X, x.getType(),i,ii)*x.at(i+ii,j,k);

		Scalar dy = 0.;
		for( int jj=space_->nl(Y); jj<=space_->nu(Y); ++jj ) 
			dy += getC( v,Y, x.getType(),j,jj)*x.at(i,j+jj,k);

		Scalar dz = 0.;
		for( int kk=space_->nl(Z); kk<=space_->nu(Z); ++kk ) 
			dz += getC( w,Z, x.getType(),k,kk)*x.at(i,j,k+kk);

		return( u*dx+v*dy+w*dz );
	}

	constexpr const std::string getLabel() const { return( "Convection" ); };


}; // end of class ConvectionSOp


} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::ConvectionSOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::ConvectionSOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::ConvectionSOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::ConvectionSOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_CONVECTIONSOP_HPP
