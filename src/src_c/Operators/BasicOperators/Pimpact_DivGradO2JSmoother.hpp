#pragma once
#ifndef PIMPACT_DIVGRADO2JSMOOTHER_HPP
#define PIMPACT_DIVGRADO2JSMOOTHER_HPP


#include "Pimpact_DivGradO2Op.hpp"




namespace Pimpact{


extern "C" {

void OP_DivGradO2JSmoother(
    const int& dimens,
    const int* const N,
    const int* const BL,
    const int* const BU,
    const int* const BCL,
    const int* const BCU,
    const double* const cdg1,
    const double* const cdg2,
    const double* const cdg3,
    const double& omega,
    const double* const b,
    const double* const x,
          double* const temp );

void SF_handle_corner(
    const int* const N,
    const int* const BL,
    const int* const BU,
    const int* const BCL,
    const int* const BCU,
    double* const phi );

}



/// \brief \f$\omega\f$-Jacobian smoother for second Order DivGradOp.
///
///
/// \ingroup BaseOperator
/// \todo instead of hardcode 2nd Order it would be pretty to use new space with StencilWidth<3,2>
/// \todo add ParameterList
/// \todo handle corner
template<class ST>
class DivGradO2JSmoother {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;

protected:

  Scalar omega_;
  int nIter_;

  const Teuchos::RCP<const SpaceT> space_;
  Teuchos::RCP<DomainFieldT> temp_;
  const Teuchos::RCP<const DivGradO2Op<SpaceT> > op_;

public:

  DivGradO2JSmoother(
      const Teuchos::RCP<const DivGradO2Op<SpaceT> >& op,
      Teuchos::RCP<Teuchos::ParameterList> pl ):
    omega_( pl->get("ommega",0.8) ),
    nIter_( pl->get("numIters",10) ),
    space_(op->space_),
    temp_( createScalarField<SpaceT>( space_ ) ),
    op_(op) {}


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}( x - A y_k ) \f]
  void apply(const DomainFieldT& x, RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS ) const {

    for( int i=0; i<nIter_; ++ i) {
      y.exchange();

      OP_DivGradO2JSmoother(
          space_->dim(),
          space_->nLoc(),
          space_->bl(),
          space_->bu(),
          space_->getDomain()->getBCLocal()->getBCL(),
          space_->getDomain()->getBCLocal()->getBCU(),
          op_->c_[0],
          op_->c_[1],
          op_->c_[2],
          omega_,
          x.s_,
          y.s_,
          temp_->s_);

      SF_handle_corner(
          space_->nLoc(),
          space_->bl(),
          space_->bu(),
          space_->getDomain()->getBCLocal()->getBCL(),
          space_->getDomain()->getBCLocal()->getBCU(),
          temp_->s_);

      // attention: could lead to problems when ScalarField is used as part of a higherlevel class
      std::swap( y.s_, temp_->s_ );
      y.changed();
    }

  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

  void print( std::ostream& out=std::cout ) const {
    out << "\n --- Jacobian smoother ---\n";
    out << "\t omega: " << omega_ << "\n";
    out << "\t numIter: " << nIter_ << "\n";
    op_->print( out );
  }

}; // end of class DivGradO2JSmoother



/// \relates DivGradO2JSmoother
template<class SpaceT>
Teuchos::RCP< DivGradO2JSmoother<SpaceT> >
createDivGradO2JSmoother(
    const Teuchos::RCP<const DivGradO2Op<SpaceT> >& op,
    Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ) {

  return(
      Teuchos::rcp( new DivGradO2JSmoother<SpaceT>( op, pl ) ) );

}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DIVGRADO2JSMOOTHER_HPP
