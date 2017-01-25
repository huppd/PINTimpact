#pragma once
#ifndef PIMPACT_RESTRICTIONSFOP_HPP
#define PIMPACT_RESTRICTIONSFOP_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Pimpact_RestrictionBaseOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Space.hpp"
#include "Pimpact_Stencil.hpp"




namespace Pimpact {


/// \brief Opetartor that restricts from a fine space to a coarse space
/// \todo c++fy
/// \tparam ST type of the \c Space
template<class ST>
class RestrictionSFOp : private RestrictionBaseOp<ST> {

  static const int sdim = ST::sdim;
  static const int dimension = ST::dimension;

	using Scalar = typename ST::Scalar;
  using Ordinal = typename ST::Ordinal;

public:

	using SpaceT = ST;

	using FSpaceT = SpaceT;
	using CSpaceT = SpaceT;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

	using Stenc = Stencil< Scalar, Ordinal, 1, -1, 1 >;

protected:

  Teuchos::Tuple<Stenc,3> cRS_;

	void initSF() {

		for( int i=0; i<3; ++i ) {

			cRS_[i] = Stenc( this->iimax_[i] );

			MG_getCRS(
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
		}
	}

public:

	RestrictionSFOp(
			const Teuchos::RCP<const SpaceT>& spaceF,
			const Teuchos::RCP<const SpaceT>& spaceC ):
		RestrictionBaseOp<ST>( spaceF, spaceC ) {

			initSF();
	}


	RestrictionSFOp(
			const Teuchos::RCP<const SpaceT>& spaceF,
			const Teuchos::RCP<const SpaceT>& spaceC,
			const Teuchos::Tuple<int,dimension>& np ):
		RestrictionBaseOp<ST>( spaceF, spaceC, np ) {

			initSF();
	}


	void apply( const DomainFieldT& x, RangeFieldT& y ) const {

		assert( x.getType()==y.getType() );
		assert( x.getType()==F::S );

		x.exchange();

		const int sdimens = sdim; // so ugly
		MG_restrictFW(
				sdimens,
				spaceF()->nLoc(),
				spaceF()->bl(),
				spaceF()->bu(),
				spaceC()->nLoc(),
				spaceC()->bl(),
				spaceC()->bu(),
				this->iimax_.getRawPtr(),
				this->dd_.getRawPtr(),
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
	}

	
	Teuchos::Tuple<Ordinal,dimension> getDD() const { return( this->dd_ ); };

	Teuchos::RCP<const SpaceT> spaceC() const { return( this->spaceC_ ); };
	Teuchos::RCP<const SpaceT> spaceF() const { return( this->spaceF_ ); };

	const std::string getLabel() const { return( "Restriction SF" ); };

}; // end of class RestrictionSFOp


} // end of namespace Pimpact




#endif // end of #ifndef PIMPACT_RESTRICTIONSFOP_HPP
