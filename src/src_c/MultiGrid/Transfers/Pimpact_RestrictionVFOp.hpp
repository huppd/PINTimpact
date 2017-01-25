#pragma once
#ifndef PIMPACT_RESTRICTIONVFOP_HPP
#define PIMPACT_RESTRICTIONVFOP_HPP
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Pimpact_RestrictionBaseOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Space.hpp"
#include "Pimpact_Stencil.hpp"




namespace Pimpact {




/// \brief Opetartor that restricts from a fine space to a coarse space
///
/// \tparam ST type of the \c Space
template<class ST>
class RestrictionVFOp : private RestrictionBaseOp<ST>{

  static const int dimension = ST::dimension;

	using Scalar = typename ST::Scalar;
  using Ordinal = typename ST::Ordinal;

public:

	using SpaceT = ST;

	using FSpaceT = SpaceT;
	using CSpaceT = SpaceT;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

	using StencS = Stencil< Scalar, Ordinal, 1, -1, 1 >;
	using StencV = Stencil< Scalar, Ordinal, 0,  1, 2 >;

protected:

  Teuchos::Tuple<StencS,3> cRS_;
  Teuchos::Tuple<StencV,3> cRV_;

	void initVF() {

			// ------------------------- CRS, CRV
			for( int i=0; i<3; ++i ) {
				F fi = static_cast<F>( i );

				cRS_[i] = StencS( this->iimax_[i] );

				MG_getCRVS(
						this->iimax_[i],
						(this->nGather_[i]>1)?
						spaceF()->getBCLocal()->getBCL(i):
						spaceC()->getBCLocal()->getBCL(i),
						(this->nGather_[i]>1)?
						spaceF()->getBCLocal()->getBCU(i):
						spaceC()->getBCLocal()->getBCU(i),
						this->dd_[i],
						spaceF()->getGridSizeLocal()->get(i),
						spaceF()->bl(i),
						spaceF()->bu(i),
						spaceF()->getCoordinatesLocal()->getX( F::S, i ),
						cRS_[i].get() );

				cRV_[i] = StencV( this->iimax_[i] );

				MG_getCRV(
						spaceF()->getGridSizeLocal()->get(i),
						spaceC()->bl(i),
						spaceC()->bu(i),
						this->iimax_[i],
						spaceC()->getBCLocal()->getBCL(i),
						spaceC()->getBCLocal()->getBCU(i),
						spaceF()->getCoordinatesLocal()->getX( fi, i ),
						spaceF()->getCoordinatesLocal()->getX( F::S, i ),
						this->dd_[i],
						cRV_[i].get() );
			}
	}

public:

	RestrictionVFOp(
			const Teuchos::RCP<const SpaceT>& spaceF,
			const Teuchos::RCP<const SpaceT>& spaceC ):
		RestrictionBaseOp<ST>( spaceF, spaceC ) {

			initVF();
  }


	RestrictionVFOp(
			const Teuchos::RCP<const SpaceT>& spaceF,
			const Teuchos::RCP<const SpaceT>& spaceC,
		  const Teuchos::Tuple<int,dimension>& np ):
		RestrictionBaseOp<ST>( spaceF, spaceC, np ) {

			initVF();
  }



	void apply( const DomainFieldT& x, RangeFieldT& y ) const {

		assert( x.getType()==y.getType() );
		assert( x.getType()!=F::S );

		F fType  = x.getType();
		int dir = static_cast<int>( fType );
		x.exchange( );

		const int sdimens = ST::sdim;
		MG_restrictFWV(
				sdimens,
				dir+1,
				spaceF()->nLoc(),
				spaceF()->bl(),
				spaceF()->bu(),
				spaceF()->sIndB(fType),
				spaceF()->eIndB(fType),
				spaceC()->nLoc(),
				spaceC()->bl(),
				spaceC()->bu(),
				spaceC()->sIndB(fType),
				spaceC()->eIndB(fType),
				this->iimax_.getRawPtr(),
				this->dd_.getRawPtr(),
				cRV_[dir].get(),
				cRS_[0].get(),
				cRS_[1].get(),
				cRS_[2].get(),
				x.getConstRawPtr(),
				y.getRawPtr() );

		this->gather( y.getRawPtr() );

		y.changed();
	}


	void print(  std::ostream& out=std::cout ) const {

		out << "=== Restriction OP ===\n";
		out << "nGather:\t" << this->nGather_ << "\n";
		out << "rankc2:\t" << this->rankc2_ << "\n";
		out << "comm2:\t" << this->comm2_ << "\n";

		out << " --- scalar stencil: ---";
		for( int j=0; j<3; ++j ) {
			out << "\ndir: " << j << "\n";
			cRS_[j].print( out );
		}

		out << " --- velocity stencil: ---";
		for( int j=0; j<3; ++j ) {
			out << "\ndir: " << j << "\n";
			cRV_[j].print( out );
		}
	}

	
	Teuchos::Tuple<Ordinal,dimension> getDD() const { return( this->dd_ ); };

	Teuchos::RCP<const SpaceT> spaceC() const { return( this->spaceC_ ); };
	Teuchos::RCP<const SpaceT> spaceF() const { return( this->spaceF_ ); };

	const std::string getLabel() const { return( "Restriction VF" ); };


}; // end of class RestrictionVFOp



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_RESTRICTIONVFOP_HPP
