#pragma once
#ifndef PIMPACT_CONVECTIONDIFFUSIONSOP_HPP
#define PIMPACT_CONVECTIONDIFFUSIONSOP_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_ScalarField.hpp"

#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_ConvectionSOp.hpp"




namespace Pimpact {


extern "C" {

void OP_convectionDiffusion(
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
    const double* const phi,
    const double* nlu,
    const double& mulI,
    const double& mulL,
    const double& mul );

}


/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup BaseOperator
template<class ST>
class ConvectionDiffusionSOp {

  template<class OT>
  friend class ConvectionDiffusionSORSmoother;
  template<class OT>
  friend class ConvectionDiffusionJSmoother;

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef Teuchos::Tuple< Teuchos::RCP< ScalarField<SpaceT> >, 3 > FluxFieldT;
  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;

protected:

  Teuchos::RCP<const SpaceT> space_;

  Teuchos::RCP<const ConvectionSOp<SpaceT> > convSOp_;
  Teuchos::RCP<const HelmholtzOp<SpaceT> > helmOp_;

public:

  ConvectionDiffusionSOp( const Teuchos::RCP<const SpaceT>& space  ):
    space_(space),
    convSOp_( create<ConvectionSOp>(space_) ),
    helmOp_( create<HelmholtzOp>(space_) ) {};

  void assignField( const RangeFieldT& mv ) {};


  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, Scalar mul=0. ) const {

//    int m = (int)z.getType();

    TEUCHOS_TEST_FOR_EXCEPT( z.getType() != y.getType() );

    for( int i =0; i<space_->dim(); ++i )
      TEUCHOS_TEST_FOR_EXCEPT( x[i]->getType() != y.getType() );

    for( int vel_dir=0; vel_dir<space_->dim(); ++vel_dir )
      x[vel_dir]->exchange();

    y.exchange();

    // why not use default parameter 1? because one has to init z equal 0.
    if( std::abs(mul) < 1.e-16 ) {
      z.init( 0. );
      mul = 1.;
    }

    OP_convectionDiffusion(
        space_->dim(),
        space_->nLoc(),
        space_->bl(),
        space_->bu(),
        space_->nl(),
        space_->nu(),
        space_->sInd(z.getType()),
        space_->eInd(z.getType()),
        convSOp_->getCD(X,z.getType()),
        convSOp_->getCD(Y,z.getType()),
        convSOp_->getCD(Z,z.getType()),
        convSOp_->getCU(X,z.getType()),
        convSOp_->getCU(Y,z.getType()),
        convSOp_->getCU(Z,z.getType()),
//        convSOp_->getCD(X,U),
//        convSOp_->getCD(Y,V),
//        convSOp_->getCD(Z,W),
//        convSOp_->getCU(X,U),
//        convSOp_->getCU(Y,V),
//        convSOp_->getCU(Z,W),
        helmOp_->getC(X,z.getType()),
        helmOp_->getC(Y,z.getType()),
        helmOp_->getC(Z,z.getType()),
        x[0]->getConstRawPtr(),
        x[1]->getConstRawPtr(),
        x[2]->getConstRawPtr(),
        y.getConstRawPtr(),
        z.getRawPtr(),
//        mul,
        1.,
        1./space_->getDomain()->getDomainSize()->getRe(),
        mul );

    z.changed();

  }

  void print( std::ostream& out=std::cout ) const {
    convSOp_->print(out);
    helmOp_->print(out);
   }


  bool hasApplyTranspose() const { return( false ); }

  Teuchos::RCP<const SpaceT>  space() const { return( space_ ); }

}; // end of class ConvectionDiffusionSOp




} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONDIFFUSIONSOP_HPP
