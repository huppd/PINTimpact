#pragma once
#ifndef PIMPACT_FORCINGOP_HPP
#define PIMPACT_FORCINGOP_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_Types.hpp"




namespace Pimpact{


/// \brief forcing operator.
/// \ingroup Operator
template<class Field>
class ForcingOp {

  using Scalar = typename Field::Scalar;
  using Ordinal = typename Field::Ordinal;

public:

  using DomainFieldT = Field;
  using RangeFieldT = Field;

  using SpaceT = typename DomainFieldT::SpaceT;

protected:

  Teuchos::RCP<DomainFieldT> forcing_;
  Scalar mul_;

public:

  ForcingOp( const Teuchos::RCP<DomainFieldT>& forcing=Teuchos::null, Scalar mul=1 ):
    forcing_(forcing),mul_(mul) {};

  void setForcing( const Teuchos::RCP<DomainFieldT> forcing ) {
    forcing_ = forcing;
  }

  void setMultiplicator( Scalar mul ){
    mul_ = mul;
  }

  /// \brief \f[ y = force*x \f]
  void apply( const DomainFieldT& x, RangeFieldT& y ) const {
		if( std::abs(mul_-1.) < 1.e-16 ) {
			y.assign( x );
		}
		else
			y.add( mul_, x, 0., y );
    y.scale( *forcing_ );
  }

	Teuchos::RCP<const SpaceT> space() const { return(forcing_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void assignField( const DomainFieldT& field ) {};

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "Forcing" ); };

  void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
  }

}; // end of ForcingOp



/// \relates ForcingOp
template<class F>
Teuchos::RCP<ForcingOp<F> > createForcingOp(
    const Teuchos::RCP<F>& forcing, typename F::Scalar mul=1. ) {
  return(
      Teuchos::rcp( new ForcingOp<F>( forcing, mul ) )
  );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_FORCINGOP_HPP
