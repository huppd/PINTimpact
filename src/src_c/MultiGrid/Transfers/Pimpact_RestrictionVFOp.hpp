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

	/// \todo mv MG_getCRVS to Base class deal with BC here
	void initVF() {

			// ------------------------- CRS, CRV
			for( int dir=0; dir<3; ++dir ) {

				const Ordinal& iimax = this->iimax_[dir];

				cRS_[dir] = StencS( iimax );

				MG_getCRVS(
						this->iimax_[dir],
						(this->nGather_[dir]>1)?
						spaceF()->getBCLocal()->getBCL(dir):
						spaceC()->getBCLocal()->getBCL(dir),
						(this->nGather_[dir]>1)?
						spaceF()->getBCLocal()->getBCU(dir):
						spaceC()->getBCLocal()->getBCU(dir),
						this->dd_[dir],
						spaceF()->getGridSizeLocal()->get(dir),
						spaceF()->bl(dir),
						spaceF()->bu(dir),
						spaceF()->getCoordinatesLocal()->getX( F::S, dir ),
						cRS_[dir].get() );


				cRV_[dir] = StencV( iimax );

				const auto& xv = spaceF()->getCoordinatesLocal()->getV(dir);
				const auto& xs = spaceF()->getCoordinatesLocal()->getS(dir);

				// Restriktion, linienweise, 1d 
				// fine
				//     xv(i)        xs(i+1)        xv(i+1)
				//  ----->-------------o------------->-----
				//       |-----------Dx12------------|
				//                  
				// coarse
				for( Ordinal ii =0; ii<=iimax; ++ii ) {

					Ordinal i = this->dd_[dir]*( ii-1 ) + 1;

					Scalar dx12 = xv[i+1] - xv[i];

					cRV_[dir](ii,1) = ( xv[i+1]-xs[i+1] )/dx12;
					cRV_[dir](ii,2) = ( xs[i+1]-xv[i  ] )/dx12;
				}

				// a little bit shaky, please verify this when IO is ready
				if( 0<spaceC()->getBCLocal()->getBCL(dir) ) {
					Scalar x0 = spaceF()->getCoordinatesLocal()->getV(dir)[0];
					cRV_[dir](0,1) = 0.;
					cRV_[dir](0,2) = 1.;
				}
				if( BC::Symmetry==spaceC()->getBCLocal()->getBCL(dir) ) {
					cRV_[dir](0,1) = 0.;
					cRV_[dir](0,2) = 0.;
				}

				if( 0<spaceC()->getBCLocal()->getBCU(dir) ) {
					cRV_[dir](iimax,1) = 1.;
					cRV_[dir](iimax,2) = 0.;
				}
				if( BC::Symmetry==spaceC()->getBCLocal()->getBCU(dir) ) {
					cRV_[dir](iimax,1) = 0.;
					cRV_[dir](iimax,2) = 0.;
				}
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
		x.exchange( );

		switch( fType ) {
			case( F::U ) :
				{ 
					if( 2==ST::sdim ) {
						for( Ordinal kk=spaceC()->begin(fType,Z,B::Y); kk<=this->iimax_[Z]; ++kk ) {
							Ordinal k = kk;
							for( Ordinal jj=spaceC()->begin(fType,Y,B::Y); jj<=this->iimax_[Y]; ++jj ) {
								Ordinal j = this->dd_[Y]*( jj - 1 ) + 1;
								for( Ordinal ii=spaceC()->begin(fType,X,B::Y); ii<=this->iimax_[X]; ++ii ) {
									Ordinal i = this->dd_[X]*( ii - 1 ) + 1;

									y(ii,jj,kk) = 0.;

									if( 1==this->dd_[X] )
										for( int jjj=-1; jjj<=1; ++jjj )
											y(ii,jj,kk) +=
												cRS_[Y](jj,jjj)*x(i,j+jjj,k);
											
									else 
										for( int jjj=-1; jjj<=1; ++jjj )
											for( int iii=1; iii<=2; ++iii )
												y(ii,jj,kk) +=
													cRV_[X](ii,iii)*cRS_[Y](jj,jjj)*x(i+iii-1,j+jjj,k) ;
									
								}
							}
						}
					}
					else {
						for( Ordinal kk=spaceC()->begin(fType,Z,B::Y); kk<=this->iimax_[Z]; ++kk ) {
							Ordinal k = this->dd_[Z]* ( kk - 1 ) + 1;
							for( Ordinal jj=spaceC()->begin(fType,Y,B::Y); jj<=this->iimax_[Y]; ++jj ) {
								Ordinal j = this->dd_[Y]*( jj - 1 ) + 1;
								for( Ordinal ii=spaceC()->begin(fType,X,B::Y); ii<=this->iimax_[X]; ++ii ) {
									Ordinal i = this->dd_[X]*( ii - 1 ) + 1;

									y(ii,jj,kk) = 0.;

									if( 1==this->dd_[X] )
										for( int kkk=-1; kkk<=1; ++kkk )
											for( int jjj=-1; jjj<=1; ++jjj )
												y(ii,jj,kk) +=
													cRS_[Y](jj,jjj)*cRS_[Z](kk,kkk)*x(i,j+jjj,k+kkk);
									else
										for( int kkk=-1; kkk<=1; ++kkk )
											for( int jjj=-1; jjj<=1; ++jjj )
												for( int iii=1; iii<=2; ++iii )
													y(ii,jj,kk) +=
														cRV_[X](ii,iii)*cRS_[Y](jj,jjj)*cRS_[Z](kk,kkk)*x(i+iii-1,j+jjj,k+kkk) ;
								}
							}
						}
					}
					break;
				}
			case( F::V ) :
				{

					if( 2==ST::sdim ) {
						for( Ordinal kk=spaceC()->begin(fType,Z,B::Y); kk<=this->iimax_[Z]; ++kk ) {
							Ordinal k = kk;
							for( Ordinal jj=spaceC()->begin(fType,Y,B::Y); jj<=this->iimax_[Y]; ++jj ) {
								Ordinal j = this->dd_[Y]*( jj - 1 ) + 1;
								for( Ordinal ii=spaceC()->begin(fType,X,B::Y); ii<=this->iimax_[X]; ++ii ) {
									Ordinal i = this->dd_[X]*( ii - 1 ) + 1;

									y(ii,jj,kk) = 0.;

									if( 1==this->dd_[Y] )
										for( int iii=-1; iii<=1; ++iii )
											y(ii,jj,kk) +=
												cRS_[X](ii,iii)*x(i+iii,j,k);
											
									else 
										for( int jjj=1; jjj<=2; ++jjj )
											for( int iii=-1; iii<=1; ++iii )
												y(ii,jj,kk) += 
													cRS_[X](ii,iii)*cRV_[Y](jj,jjj)*x(i+iii,j+jjj-1,k) ;
									
								}
							}
						}
					}
					else {
						for( Ordinal kk=spaceC()->begin(fType,Z,B::Y); kk<=this->iimax_[Z]; ++kk ) {
							Ordinal k = this->dd_[Z]* ( kk - 1 ) + 1;
							for( Ordinal jj=spaceC()->begin(fType,Y,B::Y); jj<=this->iimax_[Y]; ++jj ) {
								Ordinal j = this->dd_[Y]*( jj - 1 ) + 1;
								for( Ordinal ii=spaceC()->begin(fType,X,B::Y); ii<=this->iimax_[X]; ++ii ) {
									Ordinal i = this->dd_[X]*( ii - 1 ) + 1;

									y(ii,jj,kk) = 0.;

									if( 1==this->dd_[Y] )
										for( int kkk=-1; kkk<=1; ++kkk )
											for( int iii=-1; iii<=1; ++iii )
												y(ii,jj,kk) +=
													cRS_[X](ii,iii)*cRS_[Z](kk,kkk)*x(i+iii,j,k+kkk);
									else
										for( int kkk=-1; kkk<=1; ++kkk )
											for( int jjj=1; jjj<=2; ++jjj )
												for( int iii=-1; iii<=1; ++iii )
													y(ii,jj,kk) +=
														cRS_[X](ii,iii)*cRV_[Y](jj,jjj)*cRS_[Z](kk,kkk)*x(i+iii,j+jjj-1,k+kkk) ;
								}
							}
						}
					}
					break;
				}
			case( F::W ) :
				{ 
					
					for( Ordinal kk=spaceC()->begin(fType,Z,B::Y); kk<=this->iimax_[Z]; ++kk ) {
						Ordinal k = this->dd_[Z]* ( kk - 1 ) + 1;
						for( Ordinal jj=spaceC()->begin(fType,Y,B::Y); jj<=this->iimax_[Y]; ++jj ) {
							Ordinal j = this->dd_[Y]*( jj - 1 ) + 1;
							for( Ordinal ii=spaceC()->begin(fType,X,B::Y); ii<=this->iimax_[X]; ++ii ) {
								Ordinal i = this->dd_[X]*( ii - 1 ) + 1;

								y(ii,jj,kk) = 0.;

								if( 1==this->dd_[Z] )
									for( int jjj=-1; jjj<=1; ++jjj )
										for( int iii=-1; iii<=1; ++iii )
											y(ii,jj,kk) +=
												cRS_[X](ii,iii)*cRS_[Y](jj,jjj)*x(i+iii,j+jjj,k);
								else
									for( int kkk=1; kkk<=2; ++kkk )
										for( int jjj=-1; jjj<=1; ++jjj )
											for( int iii=-1; iii<=1; ++iii )
												y(ii,jj,kk) +=
													cRS_[X](ii,iii)*cRS_[Y](jj,jjj)*cRV_[Z](kk,kkk)*x(i+iii,j+jjj,k+kkk-1) ;
							}
						}
					}
				break;
				}
			case( F::S ) :
				{ // todo: throw exception
					break;
				}
		}

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
