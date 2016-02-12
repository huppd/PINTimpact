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
	using SolverT = Teuchos::SerialQRDenseSolver<Ordinal,Scalar>;

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

		for( Ordinal k=space()->sIndB(EField::S,2); k<=space()->eIndB(EField::S,2); ++k ) {
			for( Ordinal j=space()->sIndB(EField::S,1); j<=space()->eIndB(EField::S,1); ++j ) {
				for( Ordinal i=space()->sIndB(EField::S,0); i<=space()->eIndB(EField::S,0); ++i ) {
					Ordinal I = (i-space()->sIndB(EField::S,0)) +             
                      (j-space()->sIndB(EField::S,1))*cw[0] +
                      (k-space()->sIndB(EField::S,2))*cw[0]*cw[1];
					for( int o=-1; o<=1; ++o ) {
						if( (i+o)>space()->sIndB(EField::S,0) && (i+o)<space()->eIndB(EField::S,0) ) {

							Ordinal Io = (i+o-space()->sIndB(EField::S,0)) +             
											 		 (j-space()->sIndB(EField::S,1))*cw[0] +
											 		 (k-space()->sIndB(EField::S,2))*cw[0]*cw[1];
							(*A_)(I,Io) += op_->getC( X, i, o) ;
						}
						if( (j+o)>space()->sIndB(EField::S,1) && (j+o)<space()->eIndB(EField::S,1) ) {

							Ordinal Io = (i-space()->sIndB(EField::S,0)) +             
											 		 (j+o-space()->sIndB(EField::S,1))*cw[0] +
											 		 (k-space()->sIndB(EField::S,2))*cw[0]*cw[1];
							(*A_)(I,Io) += op_->getC( Y, j, o) ;
						}
						if( (k+o)>space()->sIndB(EField::S,2) && (k+o)<space()->eIndB(EField::S,2) ) {

							Ordinal Io = (i-space()->sIndB(EField::S,0)) +             
											 		 (j-space()->sIndB(EField::S,1))*cw[0] +
											 		 (k+o-space()->sIndB(EField::S,2))*cw[0]*cw[1];
							(*A_)(I,Io) += op_->getC( Z, k, o) ;
						}
					}
				}
			}
		}

		//Asov_->setMatrix( A_ );
		//Asov_->factor();
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
			cw[i] = ( space()->eIndB( EField::S, i ) - space()->sIndB( EField::S, i ) + 1 );
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

		std::cout << *B_;
		//Asov_->setVectors( X_, B_ );
		//Asov_->solve();
		//X_->multiply(  Teuchos::NO_TRANS,  Teuchos::NO_TRANS, 1., *A_, *B_, 0. );
		//*X_ = *B_;

		for( Ordinal k=space()->sIndB(EField::S,Z); k<=space()->eIndB(EField::S,Z); ++k ) {
			for( Ordinal j=space()->sIndB(EField::S,Y); j<=space()->eIndB(EField::S,Y); ++j ) {
				for( Ordinal i=space()->sIndB(EField::S,X); i<=space()->eIndB(EField::S,X); ++i ) {
					Ordinal I = (i-space()->sIndB(EField::S,X)) +             
						(j-space()->sIndB(EField::S,Y))*cw[0] +
						(k-space()->sIndB(EField::S,Z))*cw[0]*cw[1];
					y.at(i,j,k) = (*B_)(I,0);
				}
			}
		}
		//std::cout << *X_ << "\n";
		//std::cout << "Teuchos::SerialSpdDenseSolver::solve() returned : " << info << std::endl;
		//y.setCornersZero();
		//y.changed();

		//if( levelYes_ )
			//y.level();

	}

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    op_->print( out );
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
