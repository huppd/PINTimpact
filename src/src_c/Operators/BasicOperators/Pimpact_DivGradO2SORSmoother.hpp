#pragma once
#ifndef PIMPACT_DIVGRADO2SORSMOOTHER_HPP
#define PIMPACT_DIVGRADO2SORSMOOTHER_HPP


#include "Pimpact_DivGradO2Op.hpp"




namespace Pimpact{


extern "C"
void OP_DivGradO2SORSmoother(
    const int* const N,
    const int* const BL,
    const int* const BU,
    const int* const BCL,
    const int* const BCU,
    const double* const cdg1,
    const double* const cdg2,
    const double* const cdg3,
    const double& omega,
    const int& rbgsMode,
    const double* const b,
    double* const x );



/// \brief \f$\omega\f$-SORacobian smoother for second Order DivGradOp.
///
///
/// \relates DivGradO2Op
/// \ingroup BaseOperator
/// \todo instead of hardcode 2nd Order it would be pretty to use new space with StencilWidth<3,2>
template<class OperatorT>
class DivGradO2SORSmoother {

public:

  typedef typename OperatorT::SpaceT SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;

protected:

  Scalar omega_;
	int rbgsMode_;
  int nIter_;
	bool levelYes_;

  const Teuchos::RCP<const OperatorT> op_;

public:

  DivGradO2SORSmoother(
      const Teuchos::RCP<const SpaceT>& space):
    omega_( 1. ),
		rbgsMode_( 2 ),
    nIter_( 2 ),
    levelYes_( false ),
    op_( Teuchos::rcp( new OperatorT(space) ) ) {}

	/// \brief constructor
	///
	/// \param[in] op pointer to operator that is smoothed
	/// \param[in] pl  Parameter list of options for the multi grid solver.
	///   These are the options accepted by the solver manager:
	///   - "omega" - a \c Scalar damping factor. Default: for 2D 0.8 for 3D 6./7.  /
	///   - "numIters" - a \c int number of smoothing steps . Default: 4  /
	///   - "level" - a \c bool number of smoothing steps . Default: false  /
  DivGradO2SORSmoother(
      const Teuchos::RCP<const OperatorT>& op,
      const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList() ):
    omega_( pl->get<Scalar>("omega", 1. ) ),
    rbgsMode_( pl->get<int>( "RBGS mode", 2 ) ),
    nIter_( pl->get<int>( "numIters", 2 ) ),
    levelYes_( pl->get<bool>( "level", false ) ),
    op_(op) {}


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}( x - N y_k ) \f]
	void apply(const DomainFieldT& x, RangeFieldT& y,
			Belos::ETrans trans=Belos::NOTRANS ) const {

		for( int i=0; i<nIter_; ++i) {

      y.exchange();

			OP_DivGradO2SORSmoother(
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					space()->getBCLocal()->getBCL(),
					space()->getBCLocal()->getBCU(),
					op_->getC(X),
					op_->getC(Y),
					op_->getC(Z),
					omega_,
					rbgsMode_,
					x.getConstRawPtr(),
					y.getRawPtr() );

			//y.setCornersZero();
			y.changed();
		}
		if( levelYes_ )
			y.level();

	}

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    out << "\t omega: " << omega_ << "\n";
    out << "\t numIter: " << nIter_ << "\n";
    op_->print( out );
  }

	const std::string getLabel() const { return( "DivGradO2SORSmoother" ); };

}; // end of class DivGradO2SORSmoother



} // end of namespace Pimpact



#ifdef COMPILE_ETI
extern template class Pimpact::DivGradO2SORSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::DivGradO2SORSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::DivGradO2SORSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::DivGradO2SORSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,4> > >;
#endif


#endif // end of #ifndef PIMPACT_DIVGRADO2SORSMOOTHER_HPP
