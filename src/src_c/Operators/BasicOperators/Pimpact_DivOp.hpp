#pragma once
#ifndef PIMPACT_DIVOP_HPP
#define PIMPACT_DIVOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{



/// \brief Divergence operator.
/// \ingroup BaseOperator
template<class ST>
class DivOp {

public:

  using SpaceT = ST;

  using DomainFieldT = VectorField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using TO = const Teuchos::Tuple<Scalar*,3>;

  Teuchos::RCP<const SpaceT> space_;

  TO c_;

public:

	DivOp( const Teuchos::RCP<const SpaceT>& space ):
		space_(space) {

			for( int i=0; i<3; ++i ) {
				Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->du(i) - space_->dl(i) + 1);

				c_[i] = new Scalar[ nTemp ];
				if( i<space_->dim() )
					FD_getDiffCoeff(
							space_->nLoc(i),
							space_->bl(i),
							space_->bu(i),
							space_->dl(i),
							space_->du(i),
							space_->getBCLocal()->getBCL(i),
							space_->getBCLocal()->getBCU(i),
							space_->getShift(i),
							3,
							i+1,
							1,
							0,
							//true, // mapping
							false,
							space_->getStencilWidths()->getDimNcbD(i),
							space_->getStencilWidths()->getNcbD(i),
							space_->getCoordinatesLocal()->getX( i, i ),
							space_->getCoordinatesLocal()->getX( i, EField::S ),
							c_[i] );
			}
  };


	~DivOp() {
		for( int i=0; i<3; ++i )
			delete[] c_[i];
	}


	void apply( const DomainFieldT& x, RangeFieldT& y,
			Belos::ETrans trans=Belos::NOTRANS ) const {

		for( int dir=0; dir<space_->dim(); ++dir )
			x.exchange( dir, dir );

		if( 3==space_->dim() )  {

			for( Ordinal k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				for( Ordinal j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j )
					for( Ordinal i=space()->sInd(S,Y); i<=space()->eInd(S,X); ++i )
						y.at(i,j,k) = innerStenc3D( x, i, j, k );
		}
		else{

			for( Ordinal k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				for( Ordinal j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j )
					for( Ordinal i=space()->sInd(S,Y); i<=space()->eInd(S,X); ++i )
						y.at(i,j,k) = innerStenc2D( x, i, j, k );
		}

		y.changed();
	}

  void assignField( const RangeFieldT& mv ) const {};
  void assignField( const DomainFieldT& mv ) const {};

  bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(space_); };

	constexpr const Scalar* getC( const ECoord& dir ) const {
		return( c_[dir] );
	}

	constexpr const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const  {
		return( c_[dir][ off - space_->dl(dir) + i*( space_->du(dir) - space_->dl(dir) + 1) ] );
	}

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

	void print( std::ostream& out=std::cout ) const {
		out << "--- " << getLabel() << " ---\n";
		out << " --- stencil: ---";
		for( int dir=0; dir<3; ++dir ) {
			out << "\ndir: " << toString(static_cast<ECoord>(dir)) << "\n";

			for( int i=0; i<=space_->nLoc(dir); ++i ) {
				out << "\ni: " << i << "\n";
				for( int k=space_->dl(dir); k<=space_->du(dir); ++k ) {
					out << getC(static_cast<ECoord>(dir),i,k ) << ", ";
				}
				out << "\n";
			}
			out << "\n";
		}
	}

	const std::string getLabel() const { return( "Div" ); };

protected:

	inline constexpr Scalar innerStenc3D( const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar div = 0.;

		for( int ii=space_->dl(X); ii<=space_->du(X); ++ii ) 
			div += getC(X,i,ii)*x.getField(U).at(i+ii,j,k);

		for( int jj=space_->dl(Y); jj<=space_->du(Y); ++jj ) 
			div += getC(Y,j,jj)*x.getField(V).at(i,j+jj,k);

		for( int kk=space_->dl(Z); kk<=space_->du(Z); ++kk ) 
			div += getC(Z,k,kk)*x.getField(W).at(i,j,k+kk);

		return( div );
	}

	inline constexpr Scalar innerStenc2D( const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar div = 0.;

		for( int ii=space_->dl(X); ii<=space_->du(X); ++ii ) 
			div += getC(X,i,ii)*x.getField(U).at(i+ii,j,k);

		for( int jj=space_->dl(Y); jj<=space_->du(Y); ++jj ) 
			div += getC(Y,j,jj)*x.getField(V).at(i,j+jj,k);

		return( div );
	}

};


} // end of namespace Pimpact



#ifdef COMPILE_ETI
extern template class Pimpact::DivOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::DivOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::DivOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::DivOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_DIVOP_HPP
