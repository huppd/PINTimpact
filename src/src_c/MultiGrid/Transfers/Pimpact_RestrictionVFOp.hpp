#pragma once
#ifndef PIMPACT_RESTRICTIONVFOP_HPP
#define PIMPACT_RESTRICTIONVFOP_HPP
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Pimpact_RestrictionBaseOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Space.hpp"




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

protected:

  Teuchos::Tuple<Scalar*,3> cRS_;
  Teuchos::Tuple<Scalar*,3> cRV_;

	void initVF() {

			// ------------------------- CRS, CRV
			for( int i=0; i<3; ++i ) {

				cRS_[i] = new Scalar[ 3*this->iimax_[i]  ];

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
						spaceF()->getCoordinatesLocal()->getX( i, EField::S ),
						cRS_[i] );

				cRV_[i] = new Scalar[ 2*( this->iimax_[i]-0+1 ) ];

				MG_getCRV(
						spaceF()->getGridSizeLocal()->get(i),
						spaceC()->bl(i),
						spaceC()->bu(i),
						this->iimax_[i],
						spaceC()->getBCLocal()->getBCL(i),
						spaceC()->getBCLocal()->getBCU(i),
						spaceF()->getCoordinatesLocal()->getX( i, i ),
						spaceF()->getCoordinatesLocal()->getX( i, EField::S ),
						this->dd_[i],
						cRV_[i] );
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


  ~RestrictionVFOp() {
    for( int i=0; i<3; ++i ) {
      delete[] cRS_[i];
      delete[] cRV_[i];
    }
  }



	void apply( const DomainFieldT& x, RangeFieldT& y ) const {

		TEUCHOS_TEST_FOR_EXCEPT( x.getType()!=y.getType() );
		TEUCHOS_TEST_FOR_EXCEPT( x.getType()==EField::S );

		EField fType  = x.getType();
		int dir = fType;
		//x.exchange( dir );
		x.exchange( );

		MG_restrictFWV(
				spaceF()->dim(),
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
				cRV_[dir],
				cRS_[0],
				cRS_[1],
				cRS_[2],
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

			Ordinal nTemp1 = this->iimax_[j];
			Ordinal nTemp2 = 3;

			for( int i=0; i<nTemp1; ++i ) {
				out << "\ni: " << i+1 << "\t(";
				for( int k=0; k<nTemp2; ++k ) {
					out << cRS_[j][k+nTemp2*i] << ", ";
				}
				out << ")\n";
			}
			out << "\n";
		}

		out << " --- velocity stencil: ---";
		for( int j=0; j<3; ++j ) {

			out << "\ndir: " << j << "\n";

			Ordinal nTemp1 = spaceC()->eIndB(j)[j]-spaceC()->sIndB(j)[j]+1 ;
			Ordinal nTemp2 = 2;

			for( int i=0; i<nTemp1; ++i ) {
				out << "\ni: " << i << "\t(";
				for( int k=0; k<nTemp2; ++k ) {
					out << cRV_[j][k+nTemp2*i] << ", ";
				}
				out << ")\n";
			}
			out << "\n";
		}
	}

	
	Teuchos::Tuple<Ordinal,dimension> getDD() const { return( this->dd_ ); };

	Teuchos::RCP<const SpaceT> spaceC() const { return( this->spaceC_ ); };
	Teuchos::RCP<const SpaceT> spaceF() const { return( this->spaceF_ ); };

	const std::string getLabel() const { return( "Restriction VF" ); };


}; // end of class RestrictionVFOp



} // end of namespace Pimpact




#endif // end of #ifndef PIMPACT_RESTRICTIONVFOP_HPP
