#pragma once
#ifndef PIMPACT_DIVGRADO2INV_HPP
#define PIMPACT_DIVGRADO2INV_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialQRDenseSolver.hpp"

#include "Pimpact_DivGradO2Op.hpp"




namespace Pimpact{




/// \brief inverse for second Order DivGradOp.
///
/// \relates DivGradO2Op
/// \ingroup BaseOperator
template<class OperatorT>
class DivGradO2Inv {

public:

  using SpaceT = typename OperatorT::SpaceT;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

	using MatrixT = Teuchos::SerialDenseMatrix<Ordinal,Scalar>;
	//using SolverT = Teuchos::SerialQRDenseSolver<Ordinal,Scalar>;
	using SolverT = Teuchos::SerialDenseSolver<Ordinal,Scalar>;

	//Ordinal N_;
	bool levelYes_;
  const Teuchos::RCP<const OperatorT> op_;

	Teuchos::RCP< MatrixT > A_;
	Teuchos::RCP< SolverT > Asov_;

	Teuchos::RCP< MatrixT > X_;
	Teuchos::RCP< MatrixT > B_;

	void init() {
		Ordinal N_ = 1;
		for( int i=0; i<3; ++i ) {
			N_ *= ( space()->eIndB( EField::S, i ) - space()->sIndB( EField::S, i ) + 1 );
		}

		A_ = Teuchos::rcp( new MatrixT(N_,N_,true) );
		Asov_ = Teuchos::rcp( new SolverT() );

		X_ = Teuchos::rcp( new MatrixT(N_,1,false) );
		B_ = Teuchos::rcp( new MatrixT(N_,1,false) );

		Teuchos::Tuple<Ordinal,3> cw;
		for( int i=0; i<3; ++i ) {
			cw[i] = ( space()->eIndB( EField::S, i ) - space()->sIndB( EField::S, i ) + 1 );
		}

		// inner stencil
		for( Ordinal k=space()->sIndB(EField::S,Z); k<=space()->eIndB(EField::S,Z); ++k ) {
			for( Ordinal j=space()->sIndB(EField::S,Y); j<=space()->eIndB(EField::S,Y); ++j ) {
				for( Ordinal i=space()->sIndB(EField::S,X); i<=space()->eIndB(EField::S,X); ++i ) {

					Ordinal I = (i-space()->sIndB(EField::S,X)) +             
                      (j-space()->sIndB(EField::S,Y))*cw[0] +
                      (k-space()->sIndB(EField::S,Z))*cw[0]*cw[1];
					for( int o=-1; o<=1; ++o ) {
						if( (i+o)>=space()->sIndB(EField::S,X) && (i+o)<=space()->eIndB(EField::S,X) ) {

							Ordinal Io = (i+o-space()->sIndB(EField::S,0)) +             
											 		 (j-space()->sIndB(EField::S,1))*cw[0] +
											 		 (k-space()->sIndB(EField::S,2))*cw[0]*cw[1];
							(*A_)(I,Io) += op_->getC( X, i, o) ;
						}
						if( (j+o)>=space()->sIndB(EField::S,1) && (j+o)<=space()->eIndB(EField::S,1) ) {

							Ordinal Io = (i-space()->sIndB(EField::S,0)) +             
											 		 (j+o-space()->sIndB(EField::S,1))*cw[0] +
											 		 (k-space()->sIndB(EField::S,2))*cw[0]*cw[1];
							(*A_)(I,Io) += op_->getC( Y, j, o) ;
						}
						if( (k+o)>=space()->sIndB(EField::S,2) && (k+o)<=space()->eIndB(EField::S,2) ) {

							Ordinal Io = (i-space()->sIndB(EField::S,0)) +             
											 		 (j-space()->sIndB(EField::S,1))*cw[0] +
											 		 (k+o-space()->sIndB(EField::S,2))*cw[0]*cw[1];
							(*A_)(I,Io) += op_->getC( Z, k, o) ;
						}

					}
				}
			}
		}

		// boundary conditions
		if( space()->getBCLocal()->getBCL(X)>0 ) {
			//if( BCL(1) > 0 ) then
			//i = 1
			//!do k = S3R, N3R
			//do k = S33R, N33R
			//!pgi$ unroll = n:8
			//do j = S22R, N22R
			//Lap(i,j,k) = cdg1(0,i)*phi(i,j,k) + cdg1(1,i)*phi(i+1,j,k)
			//end do
			//end do
			//end if
			Ordinal i = space()->sIndB(EField::S,X);
			for( Ordinal k=space()->sIndB(EField::S,Z); k<=space()->eIndB(EField::S,Z); ++k ) {
				for( Ordinal j=space()->sIndB(EField::S,Y); j<=space()->eIndB(EField::S,Y); ++j ) {

					Ordinal I = (i-space()->sIndB(EField::S,0)) +             
						(j-space()->sIndB(EField::S,1))*cw[0] +
						(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

					// set row zero
					for( Ordinal l=0; l<A_->numCols(); ++ l ) 
						(*A_)(I,l) = 0.;

					for( int o=-1; o<=1; ++o ) {
						if( (i+o)>=space()->sIndB(EField::S,0) && (i+o)<=space()->eIndB(EField::S,0) ) {
							Ordinal Io = (i+o-space()->sIndB(EField::S,0)) +             
								(j-space()->sIndB(EField::S,1))*cw[0] +
								(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];
							(*A_)(I,Io) = op_->getC( X, i, o) ;
						}
					}
				}
			}
		}

    //if( BCU(1) > 0 )then
      //i = N(1)
      //do k = S33R, N33R
        //!pgi$ unroll = n:8
        //do j = S22R, N22R
          //Lap(i,j,k) = cdg1(-1,i)*phi(i-1,j,k) + cdg1(0,i)*phi(i,j,k)
        //end do
      //end do
    //end if

		if( space()->getBCLocal()->getBCU(X)>0 ) {
			Ordinal i = space()->eIndB(EField::S,X);
			for( Ordinal k=space()->sIndB(EField::S,Z); k<=space()->eIndB(EField::S,Z); ++k ) {
				for( Ordinal j=space()->sIndB(EField::S,Y); j<=space()->eIndB(EField::S,Y); ++j ) {

					Ordinal I = (i-space()->sIndB(EField::S,0)) +             
						(j-space()->sIndB(EField::S,1))*cw[0] +
						(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

					// set row zero
					for( Ordinal l=0; l<A_->numCols(); ++ l ) 
						(*A_)(I,l) = 0.;

					// set stencil
					for( int o=-1; o<=1; ++o ) {
						if( (i+o)>=space()->sIndB(EField::S,0) && (i+o)<=space()->eIndB(EField::S,0) ) {
							Ordinal Io = (i+o-space()->sIndB(EField::S,0)) +             
								(j-space()->sIndB(EField::S,1))*cw[0] +
								(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];
							(*A_)(I,Io) = op_->getC( X, i, o) ;
						}
					}
				}
			}
		}

    //if( BCL(2) > 0 ) then
      //j = 1
      //do k = S33R, N33R
        //!pgi$ unroll = n:8
        //do i = S11R, N11R
          //Lap(i,j,k) = cdg2(0,j)*phi(i,j,k) + cdg2(1,j)*phi(i,j+1,k)
        //end do
      //end do
    //end if
		if( space()->getBCLocal()->getBCL(Y)>0 ) {
			Ordinal j = space()->sIndB(EField::S,Y);
			for( Ordinal k=space()->sIndB(EField::S,Z); k<=space()->eIndB(EField::S,Z); ++k ) {
				for( Ordinal i=space()->sIndB(EField::S,X); i<=space()->eIndB(EField::S,X); ++i ) {

					Ordinal I = (i-space()->sIndB(EField::S,X)) +             
						(j-space()->sIndB(EField::S,Y))*cw[0] +
						(k-space()->sIndB(EField::S,Z))*cw[0]*cw[1];
					
					// set row zero
					for( Ordinal l=0; l<A_->numCols(); ++ l ) 
						(*A_)(I,l) = 0.;

					for( int o=-1; o<=1; ++o ) {
						if( (j+o)>=space()->sIndB(EField::S,1) && (j+o)<=space()->eIndB(EField::S,1) ) {

							Ordinal Io = (i-space()->sIndB(EField::S,0)) +             
								(j+o-space()->sIndB(EField::S,1))*cw[0] +
								(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];
							(*A_)(I,Io) = op_->getC( Y, j, o) ;
						}

					}
				}
			}
		}
    //if( BCU(2) > 0 ) then
      //j = N(2)
      //do k = S33R, N33R
        //!pgi$ unroll = n:8
        //do i = S11R, N11R
          //Lap(i,j,k) = cdg2(-1,j)*phi(i,j-1,k) + cdg2(0,j)*phi(i,j,k)
        //end do
      //end do
    //end if
		if( space()->getBCLocal()->getBCU(Y)>0 ) {
			Ordinal j = space()->eIndB(EField::S,Y);
			for( Ordinal k=space()->sIndB(EField::S,Z); k<=space()->eIndB(EField::S,Z); ++k ) {
				for( Ordinal i=space()->sIndB(EField::S,X); i<=space()->eIndB(EField::S,X); ++i ) {

					Ordinal I = (i-space()->sIndB(EField::S,X)) +             
						(j-space()->sIndB(EField::S,Y))*cw[0] +
						(k-space()->sIndB(EField::S,Z))*cw[0]*cw[1];
					
					// set row zero
					for( Ordinal l=0; l<A_->numCols(); ++ l ) 
						(*A_)(I,l) = 0.;

					for( int o=-1; o<=1; ++o ) {
						if( (j+o)>=space()->sIndB(EField::S,1) && (j+o)<=space()->eIndB(EField::S,1) ) {

							Ordinal Io = (i-space()->sIndB(EField::S,0)) +             
								(j+o-space()->sIndB(EField::S,1))*cw[0] +
								(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];
							(*A_)(I,Io) = op_->getC( Y, j, o) ;
						}

					}
				}
			}
		}

    //if( BCL(3) > 0 ) then
      //k = 1
      //do j = S22R, N22R
        //!pgi$ unroll = n:8
        //do i = S11R, N11R
          //Lap(i,j,k) = cdg3(0,k)*phi(i,j,k) + cdg3(1,k)*phi(i,j,k+1)
        //end do
      //end do
    //end if
		if( space()->getBCLocal()->getBCL(Z)>0 ) {
			Ordinal k = space()->sIndB(EField::S,Z);
			for( Ordinal j=space()->sIndB(EField::S,Y); j<=space()->eIndB(EField::S,Y); ++j ) {
				for( Ordinal i=space()->sIndB(EField::S,X); i<=space()->eIndB(EField::S,X); ++i ) {

					Ordinal I = (i-space()->sIndB(EField::S,X)) +             
						(j-space()->sIndB(EField::S,Y))*cw[0] +
						(k-space()->sIndB(EField::S,Z))*cw[0]*cw[1];

					// set row zero
					for( Ordinal l=0; l<A_->numCols(); ++ l ) 
						(*A_)(I,l) = 0.;
					
					for( int o=-1; o<=1; ++o ) {
						if( (k+o)>=space()->sIndB(EField::S,2) && (k+o)<=space()->eIndB(EField::S,2) ) {

							Ordinal Io = (i-space()->sIndB(EField::S,0)) +             
								(j-space()->sIndB(EField::S,1))*cw[0] +
								(k+o-space()->sIndB(EField::S,2))*cw[0]*cw[1];
							(*A_)(I,Io) += op_->getC( Z, k, o) ;
						}
					}

				}
			}
		}

    //if( BCU(3) > 0 ) then
      //k = N(3)
      //do j = S22R, N22R
        //!pgi$ unroll = n:8
        //do i = S11R, N11R
          //Lap(i,j,k) = cdg3(-1,k)*phi(i,j,k-1) + cdg3(0,k)*phi(i,j,k)
        //end do
      //end do
    //end if
		if( space()->getBCLocal()->getBCU(Z)>0 ) {
			Ordinal k = space()->eIndB(EField::S,Z);
			for( Ordinal j=space()->sIndB(EField::S,Y); j<=space()->eIndB(EField::S,Y); ++j ) {
				for( Ordinal i=space()->sIndB(EField::S,X); i<=space()->eIndB(EField::S,X); ++i ) {

					Ordinal I = (i-space()->sIndB(EField::S,X)) +             
						(j-space()->sIndB(EField::S,Y))*cw[0] +
						(k-space()->sIndB(EField::S,Z))*cw[0]*cw[1];

					// set row zero
					for( Ordinal l=0; l<A_->numCols(); ++ l ) 
						(*A_)(I,l) = 0.;
					
					for( int o=-1; o<=1; ++o ) {
						if( (k+o)>=space()->sIndB(EField::S,2) && (k+o)<=space()->eIndB(EField::S,2) ) {

							Ordinal Io = (i-space()->sIndB(EField::S,0)) +             
								(j-space()->sIndB(EField::S,1))*cw[0] +
								(k+o-space()->sIndB(EField::S,2))*cw[0]*cw[1];
							(*A_)(I,Io) += op_->getC( Z, k, o) ;
						}
					}

				}
			}
		}

    //if( bcl(1) > 0 .and. bcl(2) > 0 ) phi(1   ,1   ,1:n(3)) = 0. ! test!!! verifizieren ...
		if( space()->getBCLocal()->getBCL(X)>0 && space()->getBCLocal()->getBCL(Y)>0 ) {
			Ordinal i = space()->sIndB(EField::S,X);
			Ordinal j = space()->sIndB(EField::S,Y);
			for( Ordinal k=space()->sIndB(EField::S,Z); k<=space()->eIndB(EField::S,Z); ++k ) {

				Ordinal I = (i-space()->sIndB(EField::S,0)) +             
					(j-space()->sIndB(EField::S,1))*cw[0] +
					(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

				// set row zero
				for( Ordinal l=0; l<A_->numCols(); ++ l ) (*A_)(I,l) = 0.;

				(*A_)(I,I) = 1.;
			}
		}
    //if( bcl(1) > 0 .and. bcu(2) > 0 ) phi(1   ,n(2),1:n(3)) = 0.
		if( space()->getBCLocal()->getBCL(X)>0 && space()->getBCLocal()->getBCU(Y)>0 ) {
			Ordinal i = space()->sIndB(EField::S,X);
			Ordinal j = space()->eIndB(EField::S,Y);
			for( Ordinal k=space()->sIndB(EField::S,Z); k<=space()->eIndB(EField::S,Z); ++k ) {

				Ordinal I = (i-space()->sIndB(EField::S,0)) +             
					(j-space()->sIndB(EField::S,1))*cw[0] +
					(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

				// set row zero
				for( Ordinal l=0; l<A_->numCols(); ++ l ) (*A_)(I,l) = 0.;

				(*A_)(I,I) = 1.;
			}
		}
    //if( bcu(1) > 0 .and. bcl(2) > 0 ) phi(n(1),1   ,1:n(3)) = 0.
		if( space()->getBCLocal()->getBCU(X)>0 && space()->getBCLocal()->getBCL(Y)>0 ) {
			Ordinal i = space()->eIndB(EField::S,X);
			Ordinal j = space()->sIndB(EField::S,Y);
			for( Ordinal k=space()->sIndB(EField::S,Z); k<=space()->eIndB(EField::S,Z); ++k ) {

				Ordinal I = (i-space()->sIndB(EField::S,0)) +             
					(j-space()->sIndB(EField::S,1))*cw[0] +
					(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

				// set row zero
				for( Ordinal l=0; l<A_->numCols(); ++ l ) (*A_)(I,l) = 0.;

				(*A_)(I,I) = 1.;
			}
		}
    //if( bcu(1) > 0 .and. bcu(2) > 0 ) phi(n(1),n(2),1:n(3)) = 0.
		if( space()->getBCLocal()->getBCU(X)>0 && space()->getBCLocal()->getBCU(Y)>0 ) {
			Ordinal i = space()->eIndB(EField::S,X);
			Ordinal j = space()->eIndB(EField::S,Y);
			for( Ordinal k=space()->sIndB(EField::S,Z); k<=space()->eIndB(EField::S,Z); ++k ) {

				Ordinal I = (i-space()->sIndB(EField::S,0)) +             
					(j-space()->sIndB(EField::S,1))*cw[0] +
					(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

				// set row zero
				for( Ordinal l=0; l<A_->numCols(); ++ l ) (*A_)(I,l) = 0.;

				(*A_)(I,I) = 1.;
			}
		}

    //if( bcl(1) > 0 .and. bcl(3) > 0 ) phi(1   ,1:n(2),1   ) = 0.
		if( space()->getBCLocal()->getBCL(X)>0 && space()->getBCLocal()->getBCL(Z)>0 ) {
			Ordinal i = space()->sIndB(EField::S,X);
			Ordinal k = space()->sIndB(EField::S,Z);
			for( Ordinal j=space()->sIndB(EField::S,Y); j<=space()->eIndB(EField::S,Y); ++j ) {

				Ordinal I = (i-space()->sIndB(EField::S,0)) +             
					(j-space()->sIndB(EField::S,1))*cw[0] +
					(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

				// set row zero
				for( Ordinal l=0; l<A_->numCols(); ++ l ) (*A_)(I,l) = 0.;

				(*A_)(I,I) = 1.;
			}
		}
    //if( bcl(1) > 0 .and. bcu(3) > 0 ) phi(1   ,1:n(2),n(3)) = 0.
		if( space()->getBCLocal()->getBCL(X)>0 && space()->getBCLocal()->getBCU(Z)>0 ) {
			Ordinal i = space()->sIndB(EField::S,X);
			Ordinal k = space()->eIndB(EField::S,Z);
			for( Ordinal j=space()->sIndB(EField::S,Y); j<=space()->eIndB(EField::S,Y); ++j ) {

				Ordinal I = (i-space()->sIndB(EField::S,0)) +             
					(j-space()->sIndB(EField::S,1))*cw[0] +
					(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

				// set row zero
				for( Ordinal l=0; l<A_->numCols(); ++ l ) (*A_)(I,l) = 0.;

				(*A_)(I,I) = 1.;
			}
		}
    //if( bcu(1) > 0 .and. bcl(3) > 0 ) phi(n(1),1:n(2),1   ) = 0.
		if( space()->getBCLocal()->getBCU(X)>0 && space()->getBCLocal()->getBCL(Z)>0 ) {
			Ordinal i = space()->eIndB(EField::S,X);
			Ordinal k = space()->sIndB(EField::S,Z);
			for( Ordinal j=space()->sIndB(EField::S,Y); j<=space()->eIndB(EField::S,Y); ++j ) {

				Ordinal I = (i-space()->sIndB(EField::S,0)) +             
					(j-space()->sIndB(EField::S,1))*cw[0] +
					(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

				// set row zero
				for( Ordinal l=0; l<A_->numCols(); ++ l ) (*A_)(I,l) = 0.;

				(*A_)(I,I) = 1.;
			}
		}
    //if( bcu(1) > 0 .and. bcu(3) > 0 ) phi(n(1),1:n(2),n(3)) = 0.
		if( space()->getBCLocal()->getBCU(X)>0 && space()->getBCLocal()->getBCU(Z)>0 ) {
			Ordinal i = space()->eIndB(EField::S,X);
			Ordinal k = space()->eIndB(EField::S,Z);
			for( Ordinal j=space()->sIndB(EField::S,Y); j<=space()->eIndB(EField::S,Y); ++j ) {

				Ordinal I = (i-space()->sIndB(EField::S,0)) +             
					(j-space()->sIndB(EField::S,1))*cw[0] +
					(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

				// set row zero
				for( Ordinal l=0; l<A_->numCols(); ++ l ) (*A_)(I,l) = 0.;

				(*A_)(I,I) = 1.;
			}
		}

    //if( bcl(2) > 0 .and. bcl(3) > 0 ) phi(1:n(1),1   ,1   ) = 0.
		if( space()->getBCLocal()->getBCL(Y)>0 && space()->getBCLocal()->getBCL(Z)>0 ) {
			Ordinal j = space()->sIndB(EField::S,Y);
			Ordinal k = space()->sIndB(EField::S,Z);
			for( Ordinal i=space()->sIndB(EField::S,X); i<=space()->eIndB(EField::S,X); ++i ) {

				Ordinal I = (i-space()->sIndB(EField::S,0)) +             
					(j-space()->sIndB(EField::S,1))*cw[0] +
					(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

				// set row zero
				for( Ordinal l=0; l<A_->numCols(); ++ l ) (*A_)(I,l) = 0.;

				(*A_)(I,I) = 1.;
			}
		}
    //if( bcl(2) > 0 .and. bcu(3) > 0 ) phi(1:n(1),1   ,n(3)) = 0.
		if( space()->getBCLocal()->getBCL(Y)>0 && space()->getBCLocal()->getBCU(Z)>0 ) {
			Ordinal j = space()->sIndB(EField::S,Y);
			Ordinal k = space()->eIndB(EField::S,Z);
			for( Ordinal i=space()->sIndB(EField::S,X); i<=space()->eIndB(EField::S,X); ++i ) {

				Ordinal I = (i-space()->sIndB(EField::S,0)) +             
					(j-space()->sIndB(EField::S,1))*cw[0] +
					(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

				// set row zero
				for( Ordinal l=0; l<A_->numCols(); ++ l ) (*A_)(I,l) = 0.;

				(*A_)(I,I) = 1.;
			}
		}
    //if( bcu(2) > 0 .and. bcl(3) > 0 ) phi(1:n(1),n(2),1   ) = 0.
		if( space()->getBCLocal()->getBCU(Y)>0 && space()->getBCLocal()->getBCL(Z)>0 ) {
			Ordinal j = space()->eIndB(EField::S,Y);
			Ordinal k = space()->sIndB(EField::S,Z);
			for( Ordinal i=space()->sIndB(EField::S,X); i<=space()->eIndB(EField::S,X); ++i ) {

				Ordinal I = (i-space()->sIndB(EField::S,0)) +             
					(j-space()->sIndB(EField::S,1))*cw[0] +
					(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

				// set row zero
				for( Ordinal l=0; l<A_->numCols(); ++ l ) (*A_)(I,l) = 0.;

				(*A_)(I,I) = 1.;
			}
		}
    //if( bcu(2) > 0 .and. bcu(3) > 0 ) phi(1:n(1),n(2),n(3)) = 0.
		if( space()->getBCLocal()->getBCU(Y)>0 && space()->getBCLocal()->getBCU(Z)>0 ) {
			Ordinal j = space()->eIndB(EField::S,Y);
			Ordinal k = space()->eIndB(EField::S,Z);
			for( Ordinal i=space()->sIndB(EField::S,X); i<=space()->eIndB(EField::S,X); ++i ) {

				Ordinal I = (i-space()->sIndB(EField::S,0)) +             
					(j-space()->sIndB(EField::S,1))*cw[0] +
					(k-space()->sIndB(EField::S,2))*cw[0]*cw[1];

				// set row zero
				for( Ordinal l=0; l<A_->numCols(); ++ l ) (*A_)(I,l) = 0.;

				(*A_)(I,I) = 1.;
			}
		}


		//std::cout << *A_;

		// set solver
		Asov_->factorWithEquilibration( true );
		Asov_->setMatrix( A_ );
		Asov_->factor();
	}

public:

	DivGradO2Inv( const Teuchos::RCP<const SpaceT>& space ):
		op_( Teuchos::rcp( new OperatorT(space) ) ) {
			init();
		}

	/// \brief constructor
	///
	/// \param[in] op pointer to operator that is smoothed
	/// \param[in] pl  Parameter list of options for the multi grid solver.
	///   These are the options accepted by the solver manager:
	DivGradO2Inv( const Teuchos::RCP<const OperatorT>& op,
      const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList() ):
    op_(op) {
			init();
		}


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}( x - N y_k ) \f]
	void apply( const DomainFieldT& x, RangeFieldT& y, Belos::ETrans
			trans=Belos::NOTRANS ) const {

		Teuchos::Tuple<Ordinal,3> cw;
		for( int i=0; i<3; ++i ) {
			cw[i] = space()->eIndB( EField::S, i ) - space()->sIndB( EField::S, i ) + 1;
		}

		x.setCornersZero();
		x.exchange();

		for( Ordinal k=space()->sIndB(EField::S,Z); k<=space()->eIndB(EField::S,Z); ++k ) {
			for( Ordinal j=space()->sIndB(EField::S,Y); j<=space()->eIndB(EField::S,Y); ++j ) {
				for( Ordinal i=space()->sIndB(EField::S,X); i<=space()->eIndB(EField::S,X); ++i ) {
					Ordinal I = (i-space()->sIndB(EField::S,X)) +             
						(j-space()->sIndB(EField::S,Y))*cw[0] +
						(k-space()->sIndB(EField::S,Z))*cw[0]*cw[1];
					(*B_)(I,0) = x.at(i,j,k);
				}
			}
		}

		//std::cout << *B_;
		Asov_->setVectors( X_, B_ );
		Asov_->solve();
		//X_->multiply(  Teuchos::NO_TRANS,  Teuchos::NO_TRANS, 1., *A_, *B_, 0. );
		//*X_ = *B_;

		for( Ordinal k=space()->sIndB(EField::S,Z); k<=space()->eIndB(EField::S,Z); ++k ) {
			for( Ordinal j=space()->sIndB(EField::S,Y); j<=space()->eIndB(EField::S,Y); ++j ) {
				for( Ordinal i=space()->sIndB(EField::S,X); i<=space()->eIndB(EField::S,X); ++i ) {
					Ordinal I = (i-space()->sIndB(EField::S,X)) +             
						(j-space()->sIndB(EField::S,Y))*cw[0] +
						(k-space()->sIndB(EField::S,Z))*cw[0]*cw[1];
					y.at(i,j,k) = (*X_)(I,0);
				}
			}
		}
		//std::cout << *X_ << "\n";
		//std::cout << "Teuchos::SerialSpdDenseSolver::solve() returned : " << info << std::endl;
		y.setCornersZero();
		y.changed();

		if( levelYes_ )
			y.level();

	}

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    op_->print( out );
		out << "\n" << *A_ << "\n";
  }

	const std::string getLabel() const { return( "DivGradO2Inv" ); };

}; // end of class DivGradO2Inv



} // end of namespace Pimpact



#ifdef COMPILE_ETI
extern template class Pimpact::DivGradO2Inv< Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::DivGradO2Inv< Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::DivGradO2Inv< Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::DivGradO2Inv< Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,4> > >;
#endif


#endif // end of #ifndef PIMPACT_DIVGRADO2INV_HPP
