#pragma once
#ifndef PIMPACT_MGSMOOTHERS_HPP
#define PIMPACT_MGSMOOTHERS_HPP


#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_Types.hpp"




namespace Pimpact {


/// \ingroup MG
/// \todo add constructro from FOperator
/// \todo generalize MGContainer, specify to Fields and Operators
template<class MGOperatorsT, template<class> class ST>
class MGSmoothers {

public:

  typedef typename MGOperatorsT::MGSpacesT MGSpacesT;
  typedef typename MGOperatorsT::COperatorT COperatorT;

  typedef ST<COperatorT> SmootherT;

protected:

//  template< template<class> class STT, class MGOpsT >
//  friend
//  Teuchos::RCP<const MGSmoothers<MGOpsT,STT> >
//  createMGSmoothers( const Teuchos::RCP<const MGOpsT>& mgOps, Teuchos::RCP<Teuchos::ParameterList> pl );

  Teuchos::RCP<const MGSpacesT> mgSpaces_;

  std::vector< Teuchos::RCP<SmootherT> >  smoothers_;

public:

  MGSmoothers(
      const Teuchos::RCP<const MGOperatorsT>& mgOperators,
      Teuchos::RCP<Teuchos::ParameterList> pl ):
    mgSpaces_( mgOperators->getMGSpaces() ),
    smoothers_( mgSpaces_->getNGrids() ) {

    for( int i=0; i<mgSpaces_->getNGrids(); ++i )
			if( mgSpaces_->participating(i) )
				smoothers_[i] = Teuchos::rcp( new SmootherT( mgOperators->get(i), pl ) );
//      smoothers_[i] = create<SmootherT>( mgOperators->get(i), pl );

	// not working on brutus
    //smoothers_.shrink_to_fit();

  }

//public:

  /// \brief gets ith smoother, similar to python i=-1 is gets you the coarses space
  Teuchos::RCP<SmootherT> get( int i ) const {
    if( i<0 )
      return( smoothers_[mgSpaces_->getNGrids()+i] );
  else
      return( smoothers_[i] );
  }

	//  void print(  std::ostream& out=std::cout ) const { }

}; // end of class MGSmoothers



/// \relates MGSmoothers
template< template<class> class SmootherT, class MGOperatorsT >
Teuchos::RCP<const MGSmoothers< MGOperatorsT, SmootherT> >
createMGSmoothers(
    const Teuchos::RCP<const MGOperatorsT>& mgOperators,
    Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ) {

  return(
      Teuchos::rcp(
          new MGSmoothers<MGOperatorsT,SmootherT>(
              mgOperators,
              pl
          )
      )
  );

}



} // end of namespace Pimpact



#include "Pimpact_MGSpaces.hpp"
#include "Pimpact_MGOperators.hpp"
#include "Pimpact_DivGradO2JSmoother.hpp"
#include "Pimpact_NonlinearSmoother.hpp"
#include "Pimpact_ConvectionDiffusionSORSmoother.hpp"
#include "Pimpact_ConvectionDiffusionJSmoother.hpp"
template<class T> using ConvDiffSORT = Pimpact::NonlinearSmoother<T,Pimpact::ConvectionDiffusionSORSmoother >;
template<class T> using ConvDiffJT = Pimpact::NonlinearSmoother<T,Pimpact::ConvectionDiffusionJSmoother >;

#ifdef COMPILE_ETI
extern template class Pimpact::MGSmoothers< Pimpact::MGOperators< Pimpact::MGSpaces< Pimpact::Space<double,int,3,4>, Pimpact::Space<double,int,3,2> >, Pimpact::DivGradOp, Pimpact::DivGradO2Op >, Pimpact::DivGradO2JSmoother >;
extern template class Pimpact::MGSmoothers< Pimpact::MGOperators< Pimpact::MGSpaces< Pimpact::Space<double,int,4,4>, Pimpact::Space<double,int,4,2> >, Pimpact::DivGradOp, Pimpact::DivGradO2Op >, Pimpact::DivGradO2JSmoother >;
#endif


#endif // end of #ifndef PIMPACT_MGSMOOTHERS_HPP
