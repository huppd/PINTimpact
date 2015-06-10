#pragma once
#ifndef PIMPACT_INTERPOLATIONOP_HPP
#define PIMPACT_INTERPOLATIONOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

#include "Teuchos_TestForException.hpp"

#include "Pimpact_Space.hpp"
#include "Pimpact_ScalarField.hpp"




namespace Pimpact {


extern "C" {

void MG_getCIS(
    const int& N,
    const int& bL,
    const int& bU,
    const double* const xs,
    double* const cI );

void MG_getCIV(
    const int& Nc,
    const int& bLc,
    const int& bUc,
    const int& SSc,
    const int& NNc,
    const int& BC_L,
    const int& BC_U,
    const int& Nf,
    const int& bLf,
    const int& bUf,
    const int& SSf,
    const double* const xc,
    const double* const xf,
    const int& dd,
    double* const cIV );

void MG_InterpolateCorners(
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
    const int* const BCL,
    const int* const BCU,
    const double* const phic );

void MG_InterpolateScatter(
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
    const int* const iimax,
    const int* const n_gather,
    const bool& participating,
    const int& rank2,
    const int& comm2,
//    const int* const recvI,
    const int* const dispI,
//    const int* const sizsI,
    const int* const offsI,
    const double* const phic );

void MG_interpolate(
    const int& dimens,
    const int* const Nc,
    const int* const bLc,
    const int* const bUc,
//    const int* const BCL,
//    const int* const BCU,
    const int* const Nf,
    const int* const bLf,
    const int* const bUf,
    const int* const iimax,
    const int* const dd,
    const double* const cI1,
    const double* const cI2,
    const double* const cI3,
    const double* const phic,
    double* const phif );

void MG_interpolateV(
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
//    const int* const BCL,
//    const int* const BCU,
    const int* const iimax,
    const int* const dd,
    const double* const cIV,
    const double* const cI1,
    const double* const cI2,
    const double* const cI3,
    const double* const phif,
    double* const phic );

}



template<class SpaceT>
class InterpolationOp {

	typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;


//  typedef Space<Scalar,Ordinal,dimension> SpaceT;

  Teuchos::RCP<const SpaceT> spaceC_;
  Teuchos::RCP<const SpaceT> spaceF_;

	int rankc2_;

	MPI_Comm comm2_;

	Teuchos::Tuple<int,3> nGather_;

	Teuchos::Tuple<Ordinal,3> iimax_;
	Teuchos::Tuple<Ordinal,3> iiShift_;
	Teuchos::Tuple<Ordinal,3> dd_;

  Ordinal* offsI_;
  Ordinal* sizsI_;
  Ordinal* recvI_;
  Ordinal* dispI_;

  Teuchos::Tuple<Scalar*,3> cIS_;
  Teuchos::Tuple<Scalar*,3> cIV_;

	void init( const Teuchos::Tuple<int,SpaceT::dimension> nb ) {

			// ------------- nGather_, iimax_
			Teuchos::Tuple<int,SpaceT::dimension> periodic = spaceF_->getDomain()->getBCGlobal()->periodic();

			int nGatherTotal = 1; 
			for( int i=0; i<3; ++i ) {
				nGather_[i] = spaceF_->getNProc(i)/spaceC_->getNProc(i); // check
				nGatherTotal *= nGather_[i]; // check

				iimax_[i] = (spaceC_->nLoc(i) - 1)/nGather_[i] + 1; // check
				dd_[i] = (spaceF_->nLoc(i) - 1)/( iimax_[i] -1 ); // check

				iiShift_[i] = ( iimax_[i] - 1 )*( ( spaceF_->procCoordinate()[i] -1 )%nGather_[i] ); // check
			}
			offsI_ = new Ordinal[3*nGatherTotal];
			sizsI_ = new Ordinal[3*nGatherTotal];
			recvI_ = new Ordinal[  nGatherTotal];
			dispI_ = new Ordinal[  nGatherTotal];

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
								MPI_Cart_create( commTemp, 3, nGather_.getRawPtr(), periodic.getRawPtr(), true, &comm2_ );
								MPI_Comm_free( &commTemp );
								int rank_comm2 = 0;
								if( spaceC_->getProcGrid()->participating() ) 
									MPI_Comm_rank( comm2_, &rank_comm2 );
								MPI_Allreduce( &rank_comm2, &rankc2_, 1, MPI_INTEGER, MPI_SUM, comm2_ );

							}

						}

				MPI_Group_free( &baseGroup );
				delete[] newRanks;
				// ------------------------- offsI_, sizsI_
				if( spaceF_->getProcGrid()->participating() )  {
					int rank_comm2;
					MPI_Comm_rank( comm2_, &rank_comm2 );

					std::vector<Ordinal> offs_global(3*nGatherTotal);
					std::vector<Ordinal> sizs_global(3*nGatherTotal);

					for( Ordinal i=0; i<3*nGatherTotal; ++i ) {
						offs_global[i] = 0;
						sizs_global[i] = 0;
					}

					for( Ordinal i=0; i<3; ++i ) {
						offs_global[ i + rank_comm2*3 ] = iiShift_[i];
						sizs_global[ i + rank_comm2*3 ] = iimax_[i]+1;
					}

					MPI_Allreduce( offs_global.data(), offsI_, 3*nGatherTotal, MPI_INTEGER, MPI_SUM, comm2_ );
					MPI_Allreduce( sizs_global.data(), sizsI_, 3*nGatherTotal, MPI_INTEGER, MPI_SUM, comm2_ );

					Ordinal counter = 0;
					for( int k=0; k<nGather_[2]; ++k )
						for( int j=0; j<nGather_[1]; ++j )
							for( int i=0; i<nGather_[0]; ++i ) {
								recvI_[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ]
									= sizsI_[ 0 + 3*( i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ) ]
									* sizsI_[ 1 + 3*( i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ) ]
									* sizsI_[ 2 + 3*( i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ) ];

								dispI_[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ] = counter;
								counter += recvI_[ i + j*nGather_[0] + k*nGather_[0]*nGather_[1] ];
							}

				}
			}
			// ------------------ cIS, cIV
			for( int i=0; i<3; ++i ) {

				cIS_[i] = new Scalar[ 2*( spaceC_->nLoc(i)-1+1 ) ];
				MG_getCIS(
						spaceC_->nLoc(i),
						spaceC_->bl(i),
						spaceC_->bu(i),
						spaceC_->getCoordinatesLocal()->getX( i, EField::S ),
						cIS_[i] );

				cIV_[i] = new Scalar[ 2*( spaceF_->nLoc(i)-0+1 ) ];
				//      if( i<spaceC_->dim() )
				
				Ordinal offset = 0;
				if( 1!=nGather_[i] )
					offset = (iimax_[i]-1)*( spaceF_->procCoordinate()[i]-1 );

				MG_getCIV(
						spaceC_->nLoc(i),
						spaceC_->bl(i),
						spaceC_->bu(i),
						spaceC_->sInd(i)[i],
						spaceC_->eInd(i)[i],
						spaceF_->getDomain()->getBCLocal()->getBCL(i),
						spaceF_->getDomain()->getBCLocal()->getBCU(i),
						spaceF_->nLoc(i),
						spaceF_->bl(i),
						spaceF_->bu(i),
						offset,
						spaceC_->getCoordinatesLocal()->getX( i, i ),
						spaceF_->getCoordinatesLocal()->getX( i, i ),
						dd_[i],
						cIV_[i] );
			}

	}

public:

	typedef SpaceT FSpaceT;
  typedef SpaceT CSpaceT;

	InterpolationOp( const Teuchos::RCP<const SpaceT>& spaceC, const Teuchos::RCP<const SpaceT>& spaceF ):
		spaceC_(spaceC),
		spaceF_(spaceF),
		comm2_(MPI_COMM_NULL) {

			auto nb = spaceF_->getProcGridSize()->getTuple();
			init( nb );

	}

	InterpolationOp(
			const Teuchos::RCP<const SpaceT>& spaceC,
			const Teuchos::RCP<const SpaceT>& spaceF,
			Teuchos::Tuple<int,SpaceT::dimension> nb ):
		spaceC_(spaceC),
		spaceF_(spaceF),
		comm2_(MPI_COMM_NULL) {

			init( nb );

	}


	~InterpolationOp() {
		for( int i=0; i<3; ++i ) {
			delete[] cIS_[i];
			delete[] cIV_[i];
		}
		delete[] offsI_;
		delete[] sizsI_;
		delete[] recvI_;
		delete[] dispI_;
	}

	void apply( const DomainFieldT& x, RangeFieldT& y ) const {

		y.initField();

		EField fType = x.getType();

		TEUCHOS_TEST_FOR_EXCEPT( x.getType()!=y.getType() );

		if( EField::S==fType ) {

			MG_InterpolateCorners(
					spaceC_->nLoc(),
					spaceC_->bl(),
					spaceC_->bu(),
					spaceC_->getDomain()->getBCLocal()->getBCL(),
					spaceC_->getDomain()->getBCLocal()->getBCU(),
					x.getConstRawPtr() );

			if( spaceC_->getProcGrid()->participating() )
				x.exchange();

			if( nGather_[0]*nGather_[1]*nGather_[2]>1 ) {
				MG_InterpolateScatter(
						spaceC_->nLoc(),
						spaceC_->bl(),
						spaceC_->bu(),
						iimax_.getRawPtr(),
						nGather_.getRawPtr(),
						spaceC_->getProcGrid()->participating(),
						rankc2_,
						MPI_Comm_c2f(comm2_),
						dispI_,          
						offsI_,          
						x.getConstRawPtr() );
			}

			MG_interpolate(
					spaceC_->dim(),
					spaceC_->nLoc(),
					spaceC_->bl(),
					spaceC_->bu(),
					spaceF_->nLoc(),
					spaceF_->bl(),
					spaceF_->bu(),
					iimax_.getRawPtr(),
					dd_.getRawPtr(),
					cIS_[0],
					cIS_[1],
					cIS_[2],
					x.getConstRawPtr(),
					y.getRawPtr() );
		}
		else {

			//		y.initField();
			int dir = fType;

			if( spaceC_->getProcGrid()->participating() ) {
				switch( fType ) {
					case EField::U:
						x.exchange(1);
						x.exchange(2);
						x.exchange(0);
						break;
					case EField::V:
						x.exchange(2);
						x.exchange(0);
						x.exchange(1);
						break;
					case EField::W:
						x.exchange(0);
						x.exchange(1);
						x.exchange(2);
					case EField::S:
						break;
				}
			}

			if( nGather_[0]*nGather_[1]*nGather_[2]>1 ) {
				MG_InterpolateScatter(
						spaceC_->nLoc(),
						spaceC_->bl(),
						spaceC_->bu(),
						iimax_.getRawPtr(),
						nGather_.getRawPtr(),
						spaceC_->getProcGrid()->participating(),
						rankc2_,
						MPI_Comm_c2f(comm2_),
						dispI_,          
						offsI_,          
						x.getConstRawPtr() );
			}

			MG_interpolateV(
					spaceC_->dim(),
					dir+1,
					spaceC_->nLoc(),
					spaceC_->bl(),
					spaceC_->bu(),
					spaceC_->sIndB(fType),
					spaceC_->eIndB(fType),
					spaceF_->nLoc(),
					spaceF_->bl(),
					spaceF_->bu(),
					spaceF_->sIndB(fType),
					spaceF_->eIndB(fType),
					iimax_.getRawPtr(),
					dd_.getRawPtr(),
					cIV_[dir],
					cIS_[0],
					cIS_[1],
					cIS_[2],
					x.getConstRawPtr(),
					y.getRawPtr() );

		}

		y.changed();

	}

  void print(  std::ostream& out=std::cout ) const {

		out << "=== Interpolation OP ===\n";
		out << "nGather:\t" << nGather_ << "\n";
		out << "dd:\t" << dd_ << "\n";
		out << "iimax:\t" << iimax_ << "\n";
		out << "rankc2:\t" << rankc2_ << "\n";
		out << "comm2:\t" << comm2_ << "\n";
    out << "\n";
    for( int j=0; j<3; ++j ) {
      out << "\n Scalar dir: " << j << ":\n";
      out << "i:\tcI(1,i)\tcI(2,i)\n";
      for( int i=0; i<( spaceC_->eInd(EField::S)[j]-spaceC_->sInd(EField::S)[j]+1 ); ++i) {
        out <<  i + spaceC_->sInd(EField::S)[j] << "\t";
        for( int k=0; k<2; ++k ) {
          out << cIS_[j][i*2+k] << "\t";
        }
        out << "\n";
      }

    }

    out << "\n";
    for( int j=0; j<3; ++j ) {
      out << "\n Vector dir: " << j << ":\n";
      out << "i:\tcV(1,i)\tcV(2,i)\n";
      for( int i=0; i<( spaceF_->nLoc(j)-0+1 ); ++i) {
        out << i << "\t";
        for( int k=0; k<2; ++k ) {
          out << cIV_[j][i*2+k] << "\t";
        }
        out << "\n";
      }
    }
    out << "\n";

  }

}; // end of class InterpolationOp




} // end of namespace Pimpact

#ifdef COMPILE_ETI
extern template class Pimpact::InterpolationOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::InterpolationOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::InterpolationOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::InterpolationOp< Pimpact::Space<double,int,4,4> >;
#endif

#endif // end of #ifndef PIMPACT_INTERPOLATIONOP_HPP
