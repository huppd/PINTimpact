#pragma once
#ifndef PIMPACT_HELMHOLTZOP_HPP
#define PIMPACT_HELMHOLTDOP_HPP

#include "Pimpact_Types.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiField.hpp"

#include "Pimpact_HelmholtzOp.hpp"
#include "Pimpact_OperatorMV.hpp"

#include "Pimpact_LinearProblem.hpp"




namespace Pimpact{


/// \ingroup BaseOperator
template< class Scalar, class Ordinal >
class Linv {
public:
  typedef VectorField<Scalar,Ordinal>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal>  RangeFieldT;
//  typedef MultiField<ModeField<VectorField<Scalar,Ordinal> > > BVF;
  typedef MultiField<VectorField<Scalar,Ordinal> > BVF;
  typedef OperatorMV<Helmholtz<Scalar,Ordinal> > LapType;
  typedef NonModeOp OpType;
private:
  Teuchos::RCP< LinearProblem<Scalar, BVF, LapType > > lap_prob_;

public:
  Linv( const Teuchos::RCP< LinearProblem<Scalar, BVF, LapType > >& lap_prob ):
    lap_prob_(lap_prob)
         {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    auto X = createMultiField<DomainFieldT>( Teuchos::rcpFromRef( const_cast<DomainFieldT&>(x) ) );
    auto Y = createMultiField<RangeFieldT>( Teuchos::rcpFromRef( y ) );
    lap_prob_->solve( Y, X );
  }

  bool hasApplyTranspose() const { return( false ); }

}; // end of class Linv


template< class Scalar, class Ordinal>
Teuchos::RCP<OperatorMV< Linv<Scalar,Ordinal> > > createLinv(
    const Teuchos::RCP<LinearProblem<Scalar,MultiField<VectorField<Scalar,Ordinal> >,OperatorMV<Helmholtz<Scalar,Ordinal> > > > lap_prob ) {

  return( Teuchos::rcp( new OperatorMV<Linv<Scalar,Ordinal> >( Teuchos::rcp( new Linv<Scalar,Ordinal>( lap_prob ) ) ) ) );
}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_HELMHOLTZOP_HPP
