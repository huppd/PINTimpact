#pragma once
#ifndef PIMPACT_OPERATORFACTORY_HPP
#define PIMPACT_OPERATORFACTORY_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_Space.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_Operator.hpp"

#include "Pimpact_MultiHarmonicMLPrec.hpp"



namespace Pimpact {


/// \relates MultiOpWrap
/// \relates ModeOpWrap
template<class Op>
Teuchos::RCP< MultiOpWrap<ModeOpWrap<Op> > > createMultiModeOpWrap( const Teuchos::RCP<Op>& op=Teuchos::null ) {
  if( Teuchos::is_null(op) )
    return( createMultiOpWrap( createModeOpWrap( Teuchos::rcp( new Op() ) ) ) );
  else
    return( createMultiOpWrap( createModeOpWrap( op ) ) );
}



/// \relates MultiOpWrap
/// \relates OperatorBase
template<class MF, class Op>
Teuchos::RCP< OperatorBase<MF> > createMultiOperatorBase( const Teuchos::RCP<Op>& op=Teuchos::null ) {
  if( op.is_null() )
    return(
        Teuchos::rcp_dynamic_cast< OperatorBase<MF> >(
            Teuchos::rcp(
                new OperatorPimpl< MF, MultiOpWrap<Op> >(
                    createMultiOpWrap<Op>( Teuchos::rcp( new Op() ) ) ) ) )
    );
  else
    return(
        Teuchos::rcp_dynamic_cast< OperatorBase<MF> >(
            Teuchos::rcp( new OperatorPimpl< MF, MultiOpWrap<Op> >( createMultiOpWrap<Op>(op) ) ) )
    );
}


/// \relates OperatorBase
/// \relates MultiOpWrap
/// \relates ModeOpWrap
template<class MF, class Op>
Teuchos::RCP<OperatorBase<MF> > createMultiModeOperatorBase( const Teuchos::RCP<Op>& op=Teuchos::null ) {
  if( Teuchos::is_null(op) )
    return(
        Teuchos::rcp_dynamic_cast< OperatorBase<MF> >(
            Teuchos::rcp(
                new OperatorPimpl< MF, MultiOpWrap<ModeOpWrap<Op> > >(
                    createMultiOpWrap(
                        createModeOpWrap(
                            Teuchos::rcp( new Op() ) ) )
                        ) ) ) );
  else
    return(
        Teuchos::rcp_dynamic_cast< OperatorBase<MF> >(
            Teuchos::rcp(
                new OperatorPimpl< MF, MultiOpWrap<ModeOpWrap<Op> > >(
                    createMultiOpWrap(
                        createModeOpWrap(op) )
                        ) ) ) );
}


/// \relates OperatorBase
/// \relates MultiOpWrap
/// \relates InverseOperator
template<class MF>
Teuchos::RCP< OperatorBase<MF> > createInverseOperatorBase( const Teuchos::RCP< LinearProblem<MF> >& linProb ) {
  return( createOperatorBase<MF,InverseOperator<MF>>( createInverseOperator<MF>(linProb) ) );
}


/// \relates OperatorBase
/// \relates MultiHarmonicMLPrec
/// \relates EddyPrec
/// \relates MLHelmholtzOp
template<class S,class O>
Teuchos::RCP< OperatorBase<MultiField<MultiHarmonicField<VectorField<S,O> > > > >
createMultiHarmonicMLEddy(
    const Teuchos::RCP<const Space<O> >& space,
    int nGrids = 20,
    int nf = 1,
    S alpha2=1.,
    S iRe = 1. ) {

  typedef VectorField<S,O> VF;
  typedef MultiField<VF> MF;
  typedef ModeField<VF> MoF;
  typedef MultiField<MoF> MMF;
  typedef MultiHarmonicField<VF> HF;
  typedef MultiField<HF> MHF;

  typedef typename MultiHarmonicMLPrec<VF>::Op0 Op0;
  typedef typename MultiHarmonicMLPrec<VF>::Ops Ops;

  Teuchos::Array< Teuchos::RCP<Ops> > array( nf );

  for( int i=0; i<nf; ++i ) {
    auto tempOp = Pimpact::createMultiModeOperatorBase<MMF>(
            Pimpact::createMLHelmholtzOp<S,O>( space, nGrids, alpha2*(i+1), iRe  ) );
    array[i] =
        Pimpact::createMultiOperatorBase<MMF>(
            Pimpact::createEddyPrec<S,O>( Teuchos::null, tempOp ) );
//    lprec = Pimpact::createMultiOperatorBase<MVF>( op2 );
//        Pimpact::createMultiOperatorBase<MultiField<ModeField<VectorField<S,O> > > >(
//            Pimpact::createModeOpWrap(
//                Pimpact::createMLHelmholtzOp<S,O>( space, nGrids, alpha2*(i+1), iRe  ) ) );
  }

  auto prec =
      Pimpact::createMultiHarmonicMLPrec<VF>(
          createMultiOperatorBase<MF>( createMLHelmholtzOp<S,O>( space, nGrids, 0., iRe ) ),
          array ) ;
  return( createMultiOperatorBase<MHF>( prec) );

//  return Teuchos::null;

}


} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_OPERATORFACTORY_HPP
