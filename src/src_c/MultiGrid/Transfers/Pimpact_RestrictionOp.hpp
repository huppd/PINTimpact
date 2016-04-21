#pragma once
#ifndef PIMPACT_RESTRICTIONOP_HPP
#define PIMPACT_RESTRICTIONOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "Teuchos_TestForException.hpp"

#include "Pimpact_Space.hpp"
#include "Pimpact_ScalarField.hpp"




namespace Pimpact {


extern "C" {

void MG_getCRS(
    const int& N,
    const int& BC_L,
    const int& BC_U,
    double* const cR );

void MG_getCRV(
    const int& N,
    const int& bL,
    const int& bU,
    const int& SS,
    const int& NN,
    const int& BC_L,
    const int& BC_U,
    const double* const xs,
    const double* const xv,
    double* const cRV );

void MG_restrict(
		const int& dimens,
    const int* const Nf,
    const int* const bLf,
    const int* const bUf,
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
		const int* const iimax,          
		const int* const dd,             
    const double* const cR1,
    const double* const cR2,
    const double* const cR3,
    const double* const phif,
    double* const phic );

void MG_restrictV(
    const int& dimens,
    const int& dir,
    const int* const Nf,
    const int* const bLf,
    const int* const bUf,
    const int* const SSf,
    const int* const NNf,
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
    const int* const SSc,
    const int* const NNc,
		const int* const iimax,          
		const int* const dd,             
//    const int& BC_L,
//    const int& BC_U,
    const double* const cRV,
    const double* const phif,
    double* const phic );


void MG_RestrictGather(
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
		const int* const iimax,          
		const int* const n_gather,       
		const bool& participate_yes,
		const int& rankc2,         
		const int& comm2,          
		const int* const recvR,          
		const int* const dispR,          
		const int* const sizsR,          
		const int* const offsR,          
    double* const phic );
}



template<class SpaceT>
class RestrictionOp {

public:

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef SpaceT FSpaceT;
  typedef SpaceT CSpaceT;

  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;

protected:

  Teuchos::RCP<const SpaceT> spaceF_;
  Teuchos::RCP<const SpaceT> spaceC_;

	int rankc2_;

	MPI_Comm comm2_;

	Teuchos::Tuple<int,3> nGather_;

	Teuchos::Tuple<Ordinal,3> iimax_;
	Teuchos::Tuple<Ordinal,3> dd_;

  Ordinal* offsR_;
  Ordinal* sizsR_;
  Ordinal* recvR_;
  Ordinal* dispR_;

  Teuchos::Tuple<Scalar*,3> cRS_;
  Teuchos::Tuple<Scalar*,3> cRV_;

	void init( const Teuchos::Tuple<int,SpaceT::dimension>& nb ) {
		
			// ------------- nGather_, iimax_
			Teuchos::Tuple<int,SpaceT::dimension> periodic = spaceF_->getDomain()->getBCGlobal()->periodic();
			Teuchos::Tuple<Ordinal,3> iiShift;

			//		bool gather_yes = false;
			int nGatherTotal = 1; 
			for( int i=0; i<3; ++i ) {
				nGather_[i] = spaceF_->getNProc(i)/spaceC_->getNProc(i);
				nGatherTotal *= nGather_[i];

				iimax_[i] = (spaceC_->nLoc(i) - 1)/nGather_[i] + 1;
				dd_[i] = (spaceF_->nLoc(i) - 1)/( iimax_[i] -1 );

				if( spaceF_->getStencilWidths()->getLS(i)==0 && (spaceF_->getDomain()->getBCGlobal()->getBCL()[i]==0 || spaceF_->getDomain()->getBCGlobal()->getBCL()[i]==-1) )
					iimax_[i] = iimax_[i]-1;
				if( spaceF_->getStencilWidths()->getLS(i)==-1 && (spaceF_->getDomain()->getBCGlobal()->getBCU()[i]==0 || spaceF_->getDomain()->getBCGlobal()->getBCU()[i]==-1) )
					iimax_[i] = iimax_[i]-1;

				iiShift[i] = ( iimax_[i] - 1 )*( ( spaceF_->procCoordinate()[i] -1 )%nGather_[i] );
				if( spaceF_->getStencilWidths()->getLS(i)==0 && (spaceF_->getDomain()->getBCGlobal()->getBCL()[i]==0 || spaceF_->getDomain()->getBCGlobal()->getBCL()[i]==-1) )
					iiShift[i] = iiShift[i] - 1;
			}
			offsR_ = new Ordinal[3*nGatherTotal];
			sizsR_ = new Ordinal[3*nGatherTotal];
			recvR_ = new Ordinal[  nGatherTotal];
			dispR_ = new Ordinal[  nGatherTotal];

			//			std::cout << "ngather: " << nGather_ << "\n";

			// ------------- rank2_, comm2_
			if( nGatherTotal>1 ) {
				int * newRanks = new int[nGatherTotal];

				MPI_Comm commWorld = spaceF_->getProcGrid()->getCommWorld();
				MPI_Comm commTemp;
				MPI_Group baseGroup, newGroup;
				MPI_Comm_group( commWorld, &baseGroup );

				for( int kk=0; kk<spaceC_->getNProc(2); ++kk )
					for( int jj=0; jj<spaceC_->getNProc(1); ++jj )
						for( int ii=0; ii<spaceC_->getNProc(0); ++ii ) {
							bool member_yes = false;

							for( int k=0; k<nGather_[2]; ++k )
								for( int j=0; j<nGather_[1]; ++j )
									for( int i=0; i<nGather_[0]; ++i ) {

										int coord[3] = {
											((ii*nGather_[0]+i)*nb[0]/spaceF_->getNProc(0))%nb[0],
											((jj*nGather_[1]+j)*nb[1]/spaceF_->getNProc(1))%nb[1],
											((kk*nGather_[2]+k)*nb[2]/spaceF_->getNProc(2))%nb[2] };

										MPI_Cart_rank( commWorld, coord, &newRanks[i+nGather_[0]*j+nGather_[0]*nGather_[1]*k]  );
										if( newRanks[i+nGather_[0]*j+nGather_[0]*nGather_[1]*k]==spaceF_->rankST() )
											member_yes = true;
									}

							MPI_Group_incl( baseGroup, nGatherTotal, newRanks, &newGroup );
							MPI_Comm_create( commWorld, newGroup, &commTemp );
							MPI_Group_free( &newGroup );

							if( member_yes ) {
								MPI_Cart_create( commTemp, 3, nGather_.getRawPtr(), periodic.getRawPtr(), false, &comm2_ );
								MPI_Comm_free( &commTemp );
								int rank_comm2 = 0;
								if( spaceC_->getProcGrid()->participating() ) 
									MPI_Comm_rank( comm2_, &rank_comm2 );
								MPI_Allreduce( &rank_comm2, &rankc2_, 1, MPI_INTEGER, MPI_SUM, comm2_ );

							}
						}

				MPI_Group_free( &baseGroup );
				delete[] newRanks;
				// ------------------------- offsR_, sizsR_

//				{
//					int size=0;
//					if( comm2_!=MPI_COMM_NULL ) {
//						MPI_Comm_size( comm2_, &size );
//						std::cout << " rank: " << spaceF_->rankST() << " comm2_sze: " << size << "\n";
//					}
//					else
//						std::cout << " rank: " << spaceF_->rankST() << " no comm\n";
//				}
				if( spaceF_->getProcGrid()->participating() )  {
//					std::cout << " rank: " << spaceF_->rankST() << " comm2_: " << comm2_ << "\n";
					int rank_comm2;
//					if( comm2_!=MPI_COMM_NULL )
					MPI_Comm_rank( comm2_, &rank_comm2 );

					std::vector<Ordinal> offs_global(3*nGatherTotal);
					std::vector<Ordinal> sizs_global(3*nGatherTotal);

					for( Ordinal i=0; i<3*nGatherTotal; ++i ) {
						offs_global[i] = 0;
						sizs_global[i] = 0;
					}


					for( Ordinal i=0; i<3; ++i ) {
						offs_global[ i + rank_comm2*3 ] = iiShift[i];
						sizs_global[ i + rank_comm2*3 ] = iimax_[i];
					}

					//				std::cout << "rank_comm2: " << rank_comm2 << "\n" << "nGather: " << nGather_ << "\n"<< "nGatherTotal: " << nGatherTotal << "\n";

					MPI_Allreduce( offs_global.data(), offsR_, 3*nGatherTotal, MPI_INTEGER, MPI_SUM, comm2_ );
					MPI_Allreduce( sizs_global.data(), sizsR_, 3*nGatherTotal, MPI_INTEGER, MPI_SUM, comm2_ );

					Ordinal counter = 0;
					for( int k=0; k<nGather_[2]; ++k )
						for( int j=0; j<nGather_[1]; ++j )
							for( int i=0; i<nGather_[0]; ++i ) {
								recvR_[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ]
									= sizsR_[ 0 + 3*( i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ) ]
									* sizsR_[ 1 + 3*( i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ) ]
									* sizsR_[ 2 + 3*( i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ) ];

								dispR_[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ] = counter;
								counter += recvR_[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ];
							}
				}
			}

			// ------------------------- CRS, CRV
			for( int i=0; i<3; ++i ) {

				cRS_[i] = new Scalar[ 3*iimax_[i]  ];
				MG_getCRS(
						iimax_[i],
						(nGather_[i]>1)?
						spaceF_->getDomain()->getBCLocal()->getBCL(i):
						spaceC_->getDomain()->getBCLocal()->getBCL(i),
						(nGather_[i]>1)?
						spaceF_->getDomain()->getBCLocal()->getBCU(i):
						spaceC_->getDomain()->getBCLocal()->getBCU(i),
						cRS_[i] );

				cRV_[i] = new Scalar[ 2*( iimax_[i]-0+1 ) ];
				MG_getCRV(
						spaceC_->getGridSizeLocal()->get(i),
						spaceC_->bl(i),
						spaceC_->bu(i),
						spaceC_->sIndB(i)[i],
						iimax_[i],
						//          spaceC_->eIndB(i)[i],
						spaceC_->getDomain()->getBCLocal()->getBCL(i),
						spaceC_->getDomain()->getBCLocal()->getBCU(i),
						spaceC_->getCoordinatesLocal()->getX( i, EField::S ),
						spaceC_->getCoordinatesLocal()->getX( i, i ),
						cRV_[i] );
			}
	}

public:

	RestrictionOp(
			const Teuchos::RCP<const SpaceT>& spaceF,
			const Teuchos::RCP<const SpaceT>& spaceC ):
		spaceF_(spaceF),
		spaceC_(spaceC),
		comm2_(MPI_COMM_NULL) {

			init( spaceF_->getProcGridSize()->getTuple() );
  }


	RestrictionOp(
			const Teuchos::RCP<const SpaceT>& spaceF,
			const Teuchos::RCP<const SpaceT>& spaceC,
		  const Teuchos::Tuple<int,SpaceT::dimension>& nb ):
		spaceF_(spaceF),
		spaceC_(spaceC),
		comm2_(MPI_COMM_NULL) {

			init( nb );
  }


  ~RestrictionOp() {
    for( int i=0; i<3; ++i ) {
      delete[] cRS_[i];
      delete[] cRV_[i];
    }
		delete[] offsR_;
		delete[] sizsR_;
		delete[] recvR_;
		delete[] dispR_;
  }



	void apply( const DomainFieldT& x, RangeFieldT& y ) const {

		TEUCHOS_TEST_FOR_EXCEPT( x.getType()!=y.getType() );

		EField fType = x.getType();


		if( EField::S==fType ) {
			x.exchange();

			MG_restrict(
					spaceF_->dim(),
					spaceF_->nLoc(),
					spaceF_->bl(),
					spaceF_->bu(),
					spaceC_->nLoc(),
					spaceC_->bl(),
					spaceC_->bu(),
					iimax_.getRawPtr(),
					dd_.getRawPtr(),
					cRS_[0],
					cRS_[1],
					cRS_[2],
					x.getConstRawPtr(),
					y.getRawPtr() );

		}
		else {
			int dir = fType;
			x.exchange( dir );

			MG_restrictV(
					spaceF_->dim(),
					dir+1,
					spaceF_->nLoc(),
					spaceF_->bl(),
					spaceF_->bu(),
					spaceF_->sIndB(fType),
					spaceF_->eIndB(fType),
					spaceC_->nLoc(),
					spaceC_->bl(),
					spaceC_->bu(),
					spaceC_->sIndB(fType),
					spaceC_->eIndB(fType),
					iimax_.getRawPtr(),
					dd_.getRawPtr(),
					cRV_[dir],
					x.getConstRawPtr(),
					y.getRawPtr() );

		}

		if( nGather_[0]*nGather_[1]*nGather_[2]>1 )
			MG_RestrictGather(
					spaceC_->nLoc(),
					spaceC_->bl(),
					spaceC_->bu(),
					iimax_.getRawPtr(),
					nGather_.getRawPtr(),
					spaceC_->getProcGrid()->participating(),
					rankc2_,
					MPI_Comm_c2f(comm2_),
					recvR_,          
					dispR_,          
					sizsR_,          
					offsR_,          
					y.getRawPtr() );

		y.changed();
	}


	void print(  std::ostream& out=std::cout ) const {

		out << "=== Restriction OP ===\n";
		out << "nGather:\t" << nGather_ << "\n";
		out << "rankc2:\t" << rankc2_ << "\n";
		out << "comm2:\t" << comm2_ << "\n";

		out << " --- scalar stencil: ---";
		for( int j=0; j<3; ++j ) {

			out << "\ndir: " << j << "\n";

			Ordinal nTemp1 = iimax_[j];
			Ordinal nTemp2 = 3;

			//for( int i=0; i<3*( spaceC_->eInd(EField::S)[j]-spaceC_->sInd(EField::S)[j]+1 ); ++i)
			//for( int k=0; k<
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

			Ordinal nTemp1 = spaceC_->eIndB(j)[j]-spaceC_->sIndB(j)[j]+1 ;
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


//	int getRankC2() const { return( rankc2_ ); }
//	MPI_Comm getComm2() const { return( comm2_ ); }

    Teuchos::RCP<const SpaceT> get_spaceC_() const { return(spaceC_); };
    Teuchos::RCP<const SpaceT> get_spaceF_() const { return(spaceF_); };

}; // end of class RestrictionOp



/// \todo colect all create methods in one file
template<template<class> class OpT, class SpaceT>
Teuchos::RCP<const OpT<SpaceT> > create(
    const Teuchos::RCP<const SpaceT>& spaceF,
    const Teuchos::RCP<const SpaceT>& spaceC ) {

  return( Teuchos::rcp( new OpT<SpaceT>(spaceF,spaceC) ) );
}


} // end of namespace Pimpact



#ifdef COMPILE_ETI
extern template class Pimpact::RestrictionOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::RestrictionOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::RestrictionOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::RestrictionOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_RESTRICTIONOP_HPP
