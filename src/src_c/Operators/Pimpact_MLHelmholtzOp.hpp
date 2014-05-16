#pragma once
#define PIMPACT_MLHELMHOLTDOP_HPP
#ifndef PIMPACT_MLHELMHOLTZOP_HPP

#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiField.hpp"

#include "Pimpact_HelmholtzOp.hpp"

#include "Pimpact_LinearProblem.hpp"



namespace Pimpact{


/// \ingroup BaseOperator
template< class Scalar, class Ordinal >
class MLHelmholtzOp {

public:

  typedef VectorField<Scalar,Ordinal>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal>  RangeFieldT;

private:

  Teuchos::RCP< LinearProblem<BVF> > lap_prob_;

public:

  MLHelmholtzOp( const Teuchos::RCP< LinearProblem<BVF> >& lap_prob=Teuchos::null ):
    lap_prob_(lap_prob)
         {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    auto X = createMultiField<DomainFieldT>( Teuchos::rcpFromRef( const_cast<DomainFieldT&>(x) ) );
    auto Y = createMultiField<RangeFieldT>( Teuchos::rcpFromRef( y ) );
    lap_prob_->solve( Y, X );
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class MLHelmholtzOp



/// \relates MLHelmholtzOp
template< class Scalar, class Ordinal>
Teuchos::RCP< MLHelmholtzOp<Scalar,Ordinal> > createMLHelmholtzOp(
    const Teuchos::RCP<LinearProblem<MultiField<VectorField<Scalar,Ordinal> > > > lap_prob ) {

  return( Teuchos::rcp( new MLHelmholtzOp<Scalar,Ordinal>( lap_prob ) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_MLHELMHOLTDOP_HPP
