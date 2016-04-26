#pragma once
#ifndef PIMPACT_GRADOP_HPP
#define PIMPACT_GRADOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{



extern "C" void OP_extrapolateBC(
		const int& m,         
    const int* const N,         
    const int* const bL,
		const int* const bU,     
    const int& dL,
		const int& dU,     
		const int* const BC_L,
		const int* const BC_U, 
		const int* const SB,
		const int* const NB,
		const double* const c,    
		const double*       phi );




/// \ingroup BaseOperator
template<class ST>
class GradOp {

public:

  using SpaceT = ST;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using TO = const Teuchos::Tuple<Scalar*,3>;

  Teuchos::RCP<const SpaceT> space_;

  TO c_;

public:

  using DomainFieldT =  ScalarField<SpaceT>;
  using RangeFieldT =  VectorField<SpaceT>;

  GradOp( const Teuchos::RCP< const SpaceT>& space):
    space_(space) {

    for( int i=0; i<3; ++i ) {
      Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->gu(i) - space_->gl(i) + 1);

      c_[i] = new Scalar[ nTemp ];

      if( i<space_->dim() )
        FD_getDiffCoeff(
//            space_->rankST(),
            space_->nLoc(i),
            space_->bl(i),
            space_->bu(i),
            space_->gl(i),
            space_->gu(i),
            space_->getBCLocal()->getBCL(i),
            space_->getBCLocal()->getBCU(i),
            space_->getShift(i),
            2,
            i+1,
            1,
            0,
						//true, // mapping
						false,
            space_->getStencilWidths()->getDimNcbG(i),
            space_->getStencilWidths()->getNcbG(i),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            space_->getCoordinatesLocal()->getX( i, i ),
            c_[i] );
    }
  };


  ~GradOp() {
    for( int i=0; i<3; ++i )
      delete[] c_[i];
  }



  void apply(const DomainFieldT& x, RangeFieldT& y) const {

		x.exchange(X);
		for( Ordinal k=space()->sIndB(U,Z); k<=space()->eIndB(U,Z); ++k )
			for( Ordinal j=space()->sIndB(U,Y); j<=space()->eIndB(U,Y); ++j )
				for( Ordinal i=space()->sIndB(U,X); i<=space()->eIndB(U,X); ++i )
					y.getField(U).at(i,j,k) = innerStencU( x, i, j, k );

		x.exchange(Y);
		for( Ordinal k=space()->sIndB(V,Z); k<=space()->eIndB(V,Z); ++k )
			for( Ordinal j=space()->sIndB(V,Y); j<=space()->eIndB(V,Y); ++j )
				for( Ordinal i=space()->sIndB(V,X); i<=space()->eIndB(V,X); ++i )
					y.getField(V).at(i,j,k) = innerStencV( x, i, j, k );

		if( 3==space_->dim() )  {

			x.exchange(Z);
			for( Ordinal k=space()->sIndB(W,Z); k<=space()->eIndB(W,Z); ++k )
				for( Ordinal j=space()->sIndB(W,Y); j<=space()->eIndB(W,Y); ++j )
					for( Ordinal i=space()->sIndB(W,X); i<=space()->eIndB(W,X); ++i )
						y.getField(W).at(i,j,k) = innerStencW( x, i, j, k );
		}

		for( int i=0; i<space()->dim(); ++i ) {
			OP_extrapolateBC(
					i+1,
					space_->nLoc(),
					space_->bl(),
					space_->bu(),
					space_->dl(i),
					space_->du(i),
					space_->getBCLocal()->getBCL(),
					space_->getBCLocal()->getBCU(),
					space_->sIndB(i),
					space_->eIndB(i),
					space_->getInterpolateV2S()->getCM( static_cast<ECoord>(i) ),
					y.getRawPtr(i) );
		}

    y.changed();
  }

  void assignField( const RangeFieldT& mv ) {};
  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(space_); };

	constexpr const Scalar* getC( const ECoord& dir ) const {
		return( c_[dir] );
	}

	constexpr const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const {
		return( c_[dir][ off - space_->gl(dir) + i*( space_->gu(dir) - space_->gl(dir) + 1) ] );
	}

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    out << " --- stencil: ---";
    for( int dir=0; dir<3; ++dir ) {
			out << "\ndir: " << toString(static_cast<ECoord>(dir)) << "\n";

      for( int i=0; i<=space_->nLoc(dir); ++i ) {
        out << "\ni: " << i << "\t(";
        for( int ii=space_->gl(dir); ii<=space_->gu(dir); ++ii ) {
          out << getC( static_cast<ECoord>(dir), i, ii ) <<", ";
        }
        out << ")\n";
      }
      out << "\n";
    }
  }

	const std::string getLabel() const { return( "Grad" ); };

protected:

	inline constexpr Scalar innerStencU( const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar grad = 0.;

		for( int ii=space_->gl(X); ii<=space_->gu(X); ++ii ) 
			grad += getC(X,i,ii)*x.at(i+ii,j,k);

		return( grad );
	}

	inline constexpr Scalar innerStencV( const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar grad = 0.;

		for( int jj=space_->gl(Y); jj<=space_->gu(Y); ++jj ) 
			grad += getC(Y,j,jj)*x.at(i,j+jj,k);

		return( grad );
	}

	inline constexpr Scalar innerStencW( const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar grad = 0.;

		for( int kk=space_->gl(Z); kk<=space_->gu(Z); ++kk ) 
			grad += getC(Z,k,kk)*x.at(i,j,k+kk);

		return( grad );
	}


}; // end of class GradOp


} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::GradOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::GradOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::GradOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::GradOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_GRADOP_HPP
