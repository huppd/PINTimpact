#pragma once
#ifndef PIMPACT_FORCINGOP_HPP
#define PIMPACT_FORCINGOP_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_Types.hpp"

//#include "Pimpact_VectorField.hpp"



namespace Pimpact{


/// \brief forcing operator.
/// \ingroup Operator
template<class Field>
class ForcingOp {

  typedef typename Field::Scalar Scalar;
  typedef typename Field::Ordinal Ordinal;

public:

  typedef Field  DomainFieldT;
  typedef Field  RangeFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

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
    y.add( mul_, x, 0., y );
    y.scale( *forcing_ );
  }

	Teuchos::RCP<const SpaceT> space() const { return(forcing_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void assignField( const DomainFieldT& field ) {};

  bool hasApplyTranspose() const { return( false ); }

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
