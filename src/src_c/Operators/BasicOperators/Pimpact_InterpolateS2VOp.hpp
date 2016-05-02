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

	using TO = const Teuchos::Tuple<Scalar*,3>;

	Teuchos::RCP< const SpaceT> space_;

	TO c_;

public:

	InterpolateS2V( const Teuchos::RCP<const SpaceT>& space ):
		space_(space) {

			for( int i=0; i<3; ++i ) {

				Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->gu(i) - space_->gl(i) + 1);
				c_[i] = new Scalar[ nTemp ];

				if( i<space_->dim() )
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
							c_[i] );
			}
		};


	~InterpolateS2V() {
		for( int i=0; i<3; ++i )
			delete[] c_[i];
	}



	void apply( const DomainFieldT& x, RangeFieldT& y ) const {

#ifndef NDEBUG
		TEUCHOS_TEST_FOR_EXCEPT( x.getType() != S );
		TEUCHOS_TEST_FOR_EXCEPT( y.getType() == S );
#endif

		int m = static_cast<int>( y.getType() );
		const EField& field = y.getType();


		x.exchange(m);
		for( Ordinal k=space()->sInd(field,Z); k<=space()->eInd(field,Z); ++k )
			for( Ordinal j=space()->sInd(field,Y); j<=space()->eInd(field,Y); ++j )
				for( Ordinal i=space()->sInd(field,X); i<=space()->eInd(field,X); ++i ) {
					y.at(i,j,k) = 0.;
					for( int ii = space_->gl(m); ii<=space_->gu(m); ++ii ) {
						if( U==field ) {
							y.at(i,j,k) += getC(static_cast<ECoord>(m),i,ii)*x.at(i+ii,j,k) ;
						}
						else if( V==field ) {
							y.at(i,j,k) += getC(static_cast<ECoord>(m),j,ii)*x.at(i,j+ii,k) ;
						}
						else if( W==field ) {
							y.at(i,j,k) += getC(static_cast<ECoord>(m),k,ii)*x.at(i,j,k+ii) ;
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

		for( int dir=0; dir<3; ++dir ) {

			out << "\ncoord: " << toString( static_cast<ECoord>(dir) ) << "\n";

			Ordinal nTemp2 = space_->gu(dir) - space_->gl(dir) + 1;

			for( Ordinal i=0; i<=space_->nLoc(dir); ++i ) {
				out << "\ni: " << i << "\t(";
				for( Ordinal k=space_->gl(dir); k<=space_->gu(dir); ++k ) {
					out << getC(static_cast<ECoord>(dir),i,k) <<", ";
				}
				out << ")\n";
			}
			out << "\n";
		}
	}

	const std::string getLabel() const { return( "InterpolateS2V" ); };

	constexpr const Scalar* getC( const ECoord& dir ) const {
		return( c_[dir] );
	}

	constexpr const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const {
		return( c_[dir][ off - space_->gl(dir) + i*( space_->gu(dir) - space_->gl(dir) + 1) ] );
	}


}; // end of class InterpolateS2V


} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::InterpolateS2V< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::InterpolateS2V< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::InterpolateS2V< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::InterpolateS2V< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_INTERPOLATES2VDOP_HPP
