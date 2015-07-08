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
    const double& mul,
    const double& mulI,
    const double& mulC,
    const double& mulL );

}


/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup NonliearOperator
template<class ST>
class ConvectionDiffusionSOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef Teuchos::Tuple< Teuchos::RCP< ScalarField<SpaceT> >, 3 > FluxFieldT;
  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;

protected:


  Teuchos::RCP<const ConvectionSOp<SpaceT> > convSOp_;
  Teuchos::RCP<const HelmholtzOp<SpaceT> > helmOp_;

	Scalar mul_;
	Scalar mulI_;
	Scalar mulC_;
	Scalar mulL_;

public:

  ConvectionDiffusionSOp( const Teuchos::RCP<const SpaceT>& space  ):
    convSOp_( create<ConvectionSOp>(space) ),
    helmOp_( create<HelmholtzOp>(space) ),
 		mul_(0.),
		mulI_(0.),
		mulC_(1.),
		mulL_(1./space()->getDomain()->getDomainSize()->getRe())	{};

  void assignField( const RangeFieldT& mv ) {};


	/// \f[ z =   (x\cdot\nabla) y - \frac{1}{Re} \Delta z \f]
  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z ) const {
		
		apply( x, y, z, mul_, mulI_, mulC_, mulL_ );

  }

	/// \f[ z = mul z + (x\cdot\nabla) y - \frac{1}{Re} \Delta z \f]
  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, Scalar mul ) const {
		
		apply( x, y, z, mul, mulI_, mulC_, mulL_ );

  }

	/// \f[ z = mul z + mulI y + mulC(x\cdot\nabla)y - mulL \Delta y \f]
  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z,
			Scalar mul, Scalar mulI, Scalar mulC, Scalar mulL ) const {

    TEUCHOS_TEST_FOR_EXCEPT( z.getType() != y.getType() );

    for( int i =0; i<space()->dim(); ++i )
      TEUCHOS_TEST_FOR_EXCEPT( x[i]->getType() != y.getType() );

    for( int vel_dir=0; vel_dir<space()->dim(); ++vel_dir )
      x[vel_dir]->exchange();

    y.exchange();

    OP_convectionDiffusion(
        space()->dim(),
        space()->nLoc(),
        space()->bl(),
        space()->bu(),
        space()->nl(),
        space()->nu(),
        space()->sInd(z.getType()),
        space()->eInd(z.getType()),
        convSOp_->getCD(X,z.getType()),
        convSOp_->getCD(Y,z.getType()),
        convSOp_->getCD(Z,z.getType()),
        convSOp_->getCU(X,z.getType()),
        convSOp_->getCU(Y,z.getType()),
        convSOp_->getCU(Z,z.getType()),
        helmOp_->getC(X,z.getType()),
        helmOp_->getC(Y,z.getType()),
        helmOp_->getC(Z,z.getType()),
        x[0]->getConstRawPtr(),
        x[1]->getConstRawPtr(),
        x[2]->getConstRawPtr(),
        y.getConstRawPtr(),
        z.getRawPtr(),
				mul,
        mulI,
				mulC,
        mulL );

    z.changed();

  }

	Teuchos::RCP<const SpaceT> space() const { return(helmOp_->space()); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		mul_ = para->get<Scalar>( "mul", 0. );
		mulI_ = para->get<Scalar>( "mulI", 0. );
		mulC_ = para->get<Scalar>( "mulC", 1. );
		mulL_ = para->get<Scalar>( "mulL", 1./space()->getDomain()->getDomainSize()->getRe() );
//		std::cout << "mul: " << mul_ << "\n";
//		std::cout << "mulI: " << mulI_ << "\n";
//		std::cout << "mulC: " << mulC_ << "\n";
//		std::cout << "mulL: " << mulL_ << "\n";
	}

  Teuchos::RCP<const ConvectionSOp<SpaceT> > getConvSOp() const { return( convSOp_ ); }
  Teuchos::RCP<const HelmholtzOp<SpaceT> > getHelmOp() const { return( helmOp_ ); }

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    convSOp_->print(out);
    helmOp_->print(out);
   }

	const Scalar& getMul() const { return( mul_ ); }
	const Scalar& getMulI() const { return( mulI_); }
	const Scalar& getMulC() const { return( mulC_); }
	const Scalar& getMulL() const { return( mulL_); }

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "ConvectionDiffusion" ); };


}; // end of class ConvectionDiffusionSOp




} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_CONVECTIONDIFFUSIONSOP_HPP
