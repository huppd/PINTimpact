#pragma once
#ifndef PIMPACT_DIVGRADO2INV_HPP
#define PIMPACT_DIVGRADO2INV_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialQRDenseSolver.hpp"

#include "Pimpact_DivGradO2Op.hpp"
#include "Pimpact_TeuchosTransfer.hpp"




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

	using VectorT = Teuchos::SerialDenseVector<Ordinal,Scalar>;
	using MatrixT = Teuchos::SerialDenseMatrix<Ordinal,Scalar>;

	//using SolverT = Teuchos::SerialQRDenseSolver<Ordinal,Scalar>;
	using SolverT = Teuchos::SerialDenseSolver<Ordinal,Scalar>;

	//Ordinal N_;
	bool levelYes_;
  const Teuchos::RCP<const OperatorT> op_;

	const Teuchos::RCP<const TeuchosTransfer<SpaceT> > trans_;

	Teuchos::RCP< MatrixT > A_;
	Teuchos::RCP< SolverT > Asov_;

	Teuchos::RCP< VectorT > X_;
	Teuchos::RCP< VectorT > B_;

	void init() {
		Ordinal N_ = trans_->getN();

		A_ = Teuchos::rcp( new MatrixT( N_, N_, true ) );
		Asov_ = Teuchos::rcp( new SolverT() );

		X_ = Teuchos::rcp( new VectorT(N_,false) );
		B_ = Teuchos::rcp( new VectorT(N_,false) );

		trans_->apply( op_, A_ );

		// set solver
		Asov_->factorWithEquilibration( true );
		Asov_->setMatrix( A_ );
		Asov_->factor();
	}

public:

	DivGradO2Inv( const Teuchos::RCP<const SpaceT>& space ):
		op_( Teuchos::rcp( new OperatorT(space) ) ),
		trans_( Teuchos::rcp( new TeuchosTransfer<SpaceT>(op_->space(), op_->space()->sIndB(EField::S), op_->space()->eIndB(EField::S) ) ) ) {
			init();
		}

	/// \brief constructor
	///
	/// \param[in] op pointer to operator that is smoothed
	/// \param[in] pl  Parameter list of options for the multi grid solver.
	///   These are the options accepted by the solver manager:
	DivGradO2Inv( const Teuchos::RCP<const OperatorT>& op,
      const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList() ):
    op_(op),
		trans_( Teuchos::rcp( new TeuchosTransfer<SpaceT>(op_->space(), op_->space()->sIndB(EField::S), op_->space()->eIndB(EField::S) ) ) ) {
			init();
		}


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}( x - N y_k ) \f]
	void apply( const DomainFieldT& x, RangeFieldT& y, Belos::ETrans
			trans=Belos::NOTRANS ) const {

		x.setCornersZero();

		//x.exchange();
		trans_->apply( x, B_ );

		y.exchange();
		trans_->updateRHS( op_, y, B_ );

		Asov_->setVectors( X_, B_ );
		Asov_->solve();
		//X_->multiply(  Teuchos::NO_TRANS,  Teuchos::NO_TRANS, 1., *A_, *B_, 0. );
		//*X_ = *B_;

		trans_->apply( X_, y);

		//std::cout << "Teuchos::SerialSpdDenseSolver::solve() returned : " << info << std::endl;
		
		y.setCornersZero(); // ???
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
