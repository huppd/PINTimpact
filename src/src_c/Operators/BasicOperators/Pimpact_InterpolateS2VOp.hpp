#pragma once
#ifndef PIMPACT_INTERPOLATES2VDOP_HPP
#define PIMPACT_INTERPOLATES2VDOP_HPP


#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact{



extern "C"
void OP_S2VOp(
		const int& dir,
		const int* const N,
		const int* const bl,
		const int* const bu,
		const int* const gl,
		const int* const gu,
		const int* const ss,
		const int* const nn,
		const double* const c,
		const double* const phi,
		double* const grad );



/// \ingroup BaseOperator
template<class ST>
class InterpolateS2V {

public:

	using SpaceT = ST;

	using Scalar = typename SpaceT::Scalar;
	using Ordinal = typename SpaceT::Ordinal;


	using DomainFieldT = ScalarField<SpaceT>;
	using RangeFieldT  = ScalarField<SpaceT>;

protected:

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

		TEUCHOS_TEST_FOR_EXCEPT( x.getType() != S );
		TEUCHOS_TEST_FOR_EXCEPT( y.getType() == S );

		int m = static_cast<int>( y.getType() );

		x.exchange( m );

		OP_S2VOp(
				m+1,
				space_->nLoc(),
				space_->bl(),
				space_->bu(),
				space_->gl(),
				space_->gu(),
				space_->sInd(m),
				space_->eInd(m),
				c_[m],
				x.getConstRawPtr(),
				y.getRawPtr() );

		y.changed();

	}

	void assignField( const RangeFieldT& mv ) {};

	bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(space_); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

	void print( std::ostream& out=std::cout ) const {

		out << "\n--- " << getLabel() << " ---\n";
		out << "--- stencil: ---";

		for( int i=0; i<3; ++i ) {

			out << "\ncoord: " << toString( static_cast<ECoord>(i) ) << "\n";

			Ordinal nTemp1 = space_->nLoc(i) + 1;
			Ordinal nTemp2 = space_->gu(i) - space_->gl(i) + 1;

			for( Ordinal j=0; j<nTemp1; ++j ) {
				out << "\ni: " << j << "\t(";
				for( Ordinal k=0; k<nTemp2; ++k ) {
					out << c_[i][k+nTemp2*j] <<", ";
				}

				out << ")\n";
			}
			out << "\n";
		}
	}

	const std::string getLabel() const { return( "InterpolateS2V" ); };

}; // end of class InterpolateS2V


} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::InterpolateS2V< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::InterpolateS2V< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::InterpolateS2V< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::InterpolateS2V< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_INTERPOLATES2VDOP_HPP
