#pragma once
#ifndef PIMPACT_CONVECTIONDIFFUSIONSORSMOOTHER_HPP
#define PIMPACT_CONVECTIONDIFFUSIONSORSMOOTHER_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_ScalarField.hpp"

#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_ConvectionSOp.hpp"




namespace Pimpact {


extern "C" {

void OP_convectionDiffusionSOR(
    const int& dimens,
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const nL,
    const int* const nU,
    const int* const SS,
    const int* const NN,
    const short int* const dir,
    const short int* const loopOrder,
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
    double* const phi,
    const double& mulI,
    const double& mulL,
    const double& mul,
    const double& om );

}


/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup NonliearOperator
template<class OperatorT>
class ConvectionDiffusionSORSmoother {

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

  int ordering_;

  Teuchos::Tuple<short int,3> dirs_;

  Teuchos::Tuple<short int,3> loopOrder_;

  Teuchos::RCP<const SpaceT> space_;

  const Teuchos::RCP<const OperatorT> op_;

public:

  ConvectionDiffusionSORSmoother(
      const Teuchos::RCP<const OperatorT>& op,
      Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ):
        omega_( pl->get("omega", 1. ) ),
        nIter_( pl->get("numIters", 1 ) ),
        ordering_( pl->get("Ordering",1 ) ),
        loopOrder_( Teuchos::tuple<short int>(1,2,3) ),
        space_(op->space_),
        op_(op) {

    if( 4==SpaceT::dimNC )
      if( 0==op->space_->rankST() )
        std::cout << "Warning!!! ConvectionDiffusionSORSmoother strange behavior for dimNC=4, problems at outflow\n";

    if( 0==ordering_ ) {
      dirs_[0] = pl->get<short int>( "dir X", 1 );
      dirs_[1] = pl->get<short int>( "dir Y", 1 );
      dirs_[2] = pl->get<short int>( "dir Z", 1 );
    }

  }


  /// \todo finding direction here.
  void assignField( const RangeFieldT& mv ) {};


  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, Scalar mul=1. ) const {

//    x[0]->write(666);
    // testing field consistency
    TEUCHOS_TEST_FOR_EXCEPT(
        z.getType() != y.getType() );

    for( int i=0; i<space_->dim(); ++i )
      TEUCHOS_TEST_FOR_EXCEPT( x[i]->getType() != y.getType() );

    // exchange wind and "rhs"
    for( int vel_dir=0; vel_dir<space_->dim(); ++vel_dir )
      x[vel_dir]->exchange();

    y.exchange();

    for( int i=0; i<nIter_; ++i ) {

      if( ordering_==0 )
        apply(x,y,z,dirs_,loopOrder_,mul);
      else
        applyNPoint(x,y,z,mul);

    }

  }

protected:

  void applyNPoint( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z,
      Scalar mul=1. ) const {

      if( 3==space_->dim() )
        for( dirs_[2]=-1; dirs_[2]<2; dirs_[2]+=2 )
          for( dirs_[1]=-1; dirs_[1]<2; dirs_[1]+=2 )
            for( dirs_[0]=-1; dirs_[0]<2; dirs_[0]+=2 )
              apply( x, y, z, dirs_, loopOrder_, mul );
      else {
        dirs_[2] = 1 ;
        for( dirs_[1]=-1; dirs_[1]<2; dirs_[1]+=2 )
          for( dirs_[0]=-1; dirs_[0]<2; dirs_[0]+=2 )
            apply( x, y, z, dirs_, loopOrder_, mul );
      }

  }


  /// \brief little helper
  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z,
      const Teuchos::Tuple<short int,3>& dirs,
      const Teuchos::Tuple<short int,3>& loopOrder,
      Scalar mul=1. ) const {

    z.exchange();
    OP_convectionDiffusionSOR(
        space_->dim(),
        space_->nLoc(),
        space_->bl(),
        space_->bu(),
        space_->nl(),
        space_->nu(),
        space_->sInd(z.getType()),
        space_->eInd(z.getType()),
        dirs.getRawPtr(),
        loopOrder.getRawPtr(),
        op_->convSOp_->getCD( X, z.getType() ),
        op_->convSOp_->getCD( Y, z.getType() ),
        op_->convSOp_->getCD( Z, z.getType() ),
        op_->convSOp_->getCU( X, z.getType() ),
        op_->convSOp_->getCU( Y, z.getType() ),
        op_->convSOp_->getCU( Z, z.getType() ),
//        op_->convSOp_->getCD( Y, (z.getType()==U)?V:U ),
//        op_->convSOp_->getCD( Y, z.getType() ),
//        op_->convSOp_->getCD( Z, z.getType() ),
//        op_->convSOp_->getCU( Y, (z.getType()==U)?V:U ),
//        op_->convSOp_->getCU( Y, z.getType() ),
//        op_->convSOp_->getCU( Z, z.getType() ),
//        op_->convSOp_->getCD( X, U ),
//        op_->convSOp_->getCD( X, U ),
//        op_->convSOp_->getCD( Z, U ),
//        op_->convSOp_->getCU( X, U ),
//        op_->convSOp_->getCU( X, U ),
//        op_->convSOp_->getCU( Z, U ),
        op_->helmOp_->getC( X,z.getType() ),
        op_->helmOp_->getC( Y,z.getType() ),
        op_->helmOp_->getC( Z,z.getType() ),
        x[0]->getRawPtr(),
        x[1]->getRawPtr(),
        x[2]->getRawPtr(),
        y.getConstRawPtr(),
        z.getRawPtr(),
        0.,
        1./space_->getDomain()->getDomainSize()->getRe(),
        mul,
        omega_ );

    z.changed();

  }

public:

  void print( std::ostream& out=std::cout ) const {

    out << "--- ConvectionDiffusionSORSmoother ---\n";
    op_->print();

  }


  bool hasApplyTranspose() const { return( false ); }

  Teuchos::RCP<const SpaceT>  space() const { return( space_ ); }

}; // end of class ConvectionDiffusionSORSmoother





} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONDIFFUSIONSORSMOOTHER_HPP
