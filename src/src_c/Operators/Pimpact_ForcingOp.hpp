#pragma once
#ifndef PIMPACT_FORCINGOP_HPP
#define PIMPACT_FORCINGOP_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"



namespace Pimpact{


/// \brief forcing operator.
/// \ingroup BaseOperator
template<class Scalar,class Ordinal>
class ForcingOp {

public:

  typedef VectorField<Scalar,Ordinal>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal>  RangeFieldT;

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

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    y.add( mul_, x, 0., y );
    y.scale( *forcing_ );
  }

  void assignField( const DomainFieldT& field ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of ForcingOp



/// \relates ForcingOp
template<class Scalar,class Ordinal>
Teuchos::RCP<ForcingOp<Scalar,Ordinal> > createForcingOp(
    const Teuchos::RCP< VectorField<Scalar,Ordinal> >& forcing, Scalar mul=1. ) {
  return(
      Teuchos::rcp( new ForcingOp<Scalar,Ordinal>( forcing, mul ) )
  );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_FORCINGOP_HPP
