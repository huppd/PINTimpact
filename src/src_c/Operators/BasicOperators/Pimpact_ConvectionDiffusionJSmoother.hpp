#pragma once
#ifndef PIMPACT_CONVECTIONDIFFUSIONJSMOOTHER_HPP
#define PIMPACT_CONVECTIONDIFFUSIONJSMOOTHER_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_ScalarField.hpp"

#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_ConvectionSOp.hpp"




namespace Pimpact {


extern "C" {

void OP_convectionDiffusionJSmoother(
    const int& dimens,
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const nL,
    const int* const nU,
    const int* const SS,
    const int* const NN,
//    const int* const dir,
    const double* const c1D,
    const double* const c2D,
    const double* const c3D,
    const double* const c1U,
    const double* const c2U,
    const double* const c3U,
    const double* const c11,
    const double* const c22,
    const double* const c33,
    const double* const phiU,
    const double* const phiV,
    const double* const phiW,
    const double* const b,
    const double* const phi,
          double* const phio,
    const double& mulI,
    const double& mulL,
    const double& mul,
    const double& om );

}


/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup BaseOperator
template<class OperatorT>
class ConvectionDiffusionJSmoother {

public:

  typedef typename OperatorT::SpaceT SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef Teuchos::Tuple< Teuchos::RCP< ScalarField<SpaceT> >, 3 > FluxFieldT;
  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;

protected:

  Scalar omega_;
  int nIter_;

  Teuchos::RCP<const SpaceT> space_;

  const Teuchos::RCP<const OperatorT> op_;

  Teuchos::RCP<DomainFieldT> temp_;

public:

  ConvectionDiffusionJSmoother(
      const Teuchos::RCP<const OperatorT>& op,
      Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ):
        omega_( pl->get("omega", 1. ) ),
        nIter_( pl->get("numIters", 4 ) ),
        space_(op->space_),
        op_(op),
        temp_( create<DomainFieldT>(op_->space_) ) {}


  void assignField( const RangeFieldT& mv ) {};


  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, Scalar mul=1. ) const {

    int m = (int)z.fType_;

    TEUCHOS_TEST_FOR_EXCEPTION(
        z.getType() != y.getType(),
        std::logic_error,
        "Pimpact::ConvectionSOP can only be applied to same fieldType !!!\n");


    for( int i =0; i<space_->dim(); ++i ) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          x[i]->fType_ != y.fType_,
          std::logic_error,
          "Pimpact::ConvectionSOP can only be applied to same fieldType !!!\n");
    }

    for( int vel_dir=0; vel_dir<space_->dim(); ++vel_dir )
      x[vel_dir]->exchange();

    for( int i=0; i<nIter_; ++i ) {

      y.exchange();
      z.exchange();

      OP_convectionDiffusionJSmoother(
          space_->dim(),
          space_->nLoc(),
          space_->bl(),
          space_->bu(),
          space_->nl(),
          space_->nu(),
          space_->sInd(m),
          space_->eInd(m),
          op_->convSOp_->getCD(X,z.getType()),
          op_->convSOp_->getCD(Y,z.getType()),
          op_->convSOp_->getCD(Z,z.getType()),
          op_->convSOp_->getCU(X,z.getType()),
          op_->convSOp_->getCU(Y,z.getType()),
          op_->convSOp_->getCU(Z,z.getType()),
          op_->helmOp_->getC(X,z.getType()),
          op_->helmOp_->getC(Y,z.getType()),
          op_->helmOp_->getC(Z,z.getType()),
          x[0]->s_,
          x[1]->s_,
          x[2]->s_,
          y.s_,
          z.s_,
          temp_->s_,
          0.,
          1./space_->getDomain()->getDomainSize()->getRe(),
          mul,
          omega_ );

      // attention: could lead to problems when ScalarField is used as part of a higherlevel class (s is shared)
      std::swap( z.s_, temp_->s_ );
      z.changed();

    }

  }

  void print( std::ostream& out=std::cout ) const {
    op_->print();
  }


  bool hasApplyTranspose() const { return( false ); }

  Teuchos::RCP<const SpaceT>  space() const { return( space_ ); }

}; // end of class ConvectionDiffusionJSmoother





} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONDIFFUSIONJSMOOTHER_HPP
