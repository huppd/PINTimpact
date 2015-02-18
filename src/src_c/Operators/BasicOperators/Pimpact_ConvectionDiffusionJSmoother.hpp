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
    const double& mulC,
    const double& mulL,
    const double& om );

}


/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup NonliearOperator
/// \todo merge with SORSmoother or make interface
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

  const Teuchos::RCP<const OperatorT> op_;

  Teuchos::RCP<DomainFieldT> temp_;

public:

	/// \brief constructor
	///
  /// These options include the following:
	/// - "omega" - damping parameter
  /// - "numIters" - an \c int specifying the maximum number of iterations the 
  ConvectionDiffusionJSmoother(
      const Teuchos::RCP<const OperatorT>& op,
      Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ):
        omega_( pl->get("omega", 0.8 ) ),
        nIter_( pl->get("numIters", 10 ) ),
        op_(op),
        temp_( create<DomainFieldT>(op_->space()) ) {}


  void assignField( const RangeFieldT& mv ) {};



  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, Scalar mul, Scalar mulI, Scalar mulC, Scalar mulL ) const { std::cout << "not implmented\n"; }

  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, Scalar mul=0. ) const {

    //int m = (int)z.getType();

    TEUCHOS_TEST_FOR_EXCEPT( z.getType() != y.getType() );


    for( int i =0; i<space()->dim(); ++i )
      TEUCHOS_TEST_FOR_EXCEPT( x[i]->getType() != y.getType() );

    for( int vel_dir=0; vel_dir<space()->dim(); ++vel_dir )
      x[vel_dir]->exchange();

    for( int i=0; i<nIter_; ++i ) {

      y.exchange();
      z.exchange();

      OP_convectionDiffusionJSmoother(
          space()->dim(),
          space()->nLoc(),
          space()->bl(),
          space()->bu(),
          space()->nl(),
          space()->nu(),
          space()->sInd(z.getType()),
          space()->eInd(z.getType()),
          op_->convSOp_->getCD(X,z.getType()),
          op_->convSOp_->getCD(Y,z.getType()),
          op_->convSOp_->getCD(Z,z.getType()),
          op_->convSOp_->getCU(X,z.getType()),
          op_->convSOp_->getCU(Y,z.getType()),
          op_->convSOp_->getCU(Z,z.getType()),
        	op_->helmOp_->getC( X, z.getType() ),
        	op_->helmOp_->getC( Y, z.getType() ),
        	op_->helmOp_->getC( Z, z.getType() ),
          x[0]->getConstRawPtr(),
          x[1]->getConstRawPtr(),
          x[2]->getRawPtr(),
          y.getConstRawPtr(),
          z.getConstRawPtr(),
          temp_->getRawPtr(),
          0.,
          1.,
          1./space()->getDomain()->getDomainSize()->getRe(),
          omega_ );

      // attention: could lead to problems when ScalarField is used as part of a higherlevel class (s is shared)
      std::swap( z.s_, temp_->s_ );
      z.changed();

    }
    if( 0!=(nIter_%2) )
      std::swap( z.s_, temp_->s_ );
  }

	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    op_->print();
  }


  bool hasApplyTranspose() const { return( false ); }


}; // end of class ConvectionDiffusionJSmoother





} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONDIFFUSIONJSMOOTHER_HPP
