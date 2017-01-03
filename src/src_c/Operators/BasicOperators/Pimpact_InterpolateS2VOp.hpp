#pragma once
#ifndef PIMPACT_INTERPOLATES2VDOP_HPP
#define PIMPACT_INTERPOLATES2VDOP_HPP


#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact{



/// \ingroup BaseOperator
template<class ST>
class InterpolateS2V {

public:

	using SpaceT = ST;

	using DomainFieldT = ScalarField<SpaceT>;
	using RangeFieldT  = ScalarField<SpaceT>;

protected:

	using Scalar = typename SpaceT::Scalar;
	using Ordinal = typename SpaceT::Ordinal;

	static const int dimNC = ST::dimNC;
	static const int dim = ST::dimension;

	using SW = StencilWidths<dim,dimNC>;

	using Stenc = Stencil< Scalar, Ordinal, 0, SW::GL(0), SW::GU(0) >;
	using TO = const Teuchos::Tuple<Stenc,ST::sdim>;

	Teuchos::RCP< const SpaceT> space_;

	TO c_;

public:

	InterpolateS2V( const Teuchos::RCP<const SpaceT>& space ):
		space_(space) {

			for( int i=0; i<SpaceT::sdim; ++i ) {

				c_[i] = Stenc( space_->nLoc(i) );

				FD_getDiffCoeff(
						space_->nLoc(i),
						space_->bl(i),
						space_->bu(i),
						space_->gl(i),
						space_->gu(i),
						space_->getBCLocal()->getBCL(i),
						space_->getBCLocal()->getBCU(i),
						space_->getShift(i),
						2,    // grid_type ???
						i+1,  // direction
						0,    // derivative
						0,    // central
						//true, // not working with grid stretching. mapping
						false, // mapping
						space_->getStencilWidths()->getDimNcbG(i),
						space_->getStencilWidths()->getNcbG(i),
						space_->getCoordinatesLocal()->getX( i, EField::S ),
						space_->getCoordinatesLocal()->getX( i, i ),
						c_[i].get() );
			}
		};



	void apply( const DomainFieldT& x, RangeFieldT& y, const Add& add=Add::N ) const {

		assert( x.getType() == S );
		assert( y.getType() != S );
		assert( !(y.getType() == W && SpaceT::sdim==2) );

		int m = static_cast<int>( y.getType() );
		const EField& field = y.getType();


		x.exchange(m);
		//
		for( Ordinal k=space()->begin(field,Z,B::Y); k<=space()->end(field,Z,B::Y); ++k )
			for( Ordinal j=space()->begin(field,Y,B::Y); j<=space()->end(field,Y,B::Y); ++j )
				for( Ordinal i=space()->begin(field,X,B::Y); i<=space()->end(field,X,B::Y); ++i ) {
					if( Add::N==add ) y(i,j,k) = 0.;
					for( int ii = space_->gl(m); ii<=space_->gu(m); ++ii ) {
						if( U==field ) {
							y(i,j,k) += getC(static_cast<ECoord>(m),i,ii)*x(i+ii,j,k) ;
						}
						else if( V==field ) {
							y(i,j,k) += getC(static_cast<ECoord>(m),j,ii)*x(i,j+ii,k) ;
						}
						else if( W==field ) {
							y(i,j,k) += getC(static_cast<ECoord>(m),k,ii)*x(i,j,k+ii) ;
						}
					}
				}

		y.changed();
	}

	void assignField( const RangeFieldT& mv ) {};

	bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(space_); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

	void print( std::ostream& out=std::cout ) const {

		out << "\n--- " << getLabel() << " ---\n";
		out << "--- stencil: ---";

		for( int dir=0; dir<SpaceT::sdim; ++dir ) {
			out << "\ncoord: " << toString( static_cast<ECoord>(dir) ) << "\n";
			c_[dir].print( out );
		}
	}

	const std::string getLabel() const { return( "InterpolateS2V" ); };

	constexpr const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const {
		return( c_[dir]( i, off ) );
	}


}; // end of class InterpolateS2V


} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_INTERPOLATES2VDOP_HPP
