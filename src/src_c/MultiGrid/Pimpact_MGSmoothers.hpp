#pragma once
#ifndef PIMPACT_MGSMOOTHERS_HPP
#define PIMPACT_MGSMOOTHERS_HPP


#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_Utils.hpp"




namespace Pimpact {


/// \ingroup MG
/// \note todo add constructro from FOperator
/// \note todo generalize MGContainer, specify to Fields and Operators
template<class MGOperatorsT, template<class> class ST>
class MGSmoothers {

public:

  using MGSpacesT = typename MGOperatorsT::MGSpacesT;
  using COperatorT = typename MGOperatorsT::COperatorT;

  using SmootherT = ST<COperatorT>;

protected:

  Teuchos::RCP<const MGSpacesT> mgSpaces_;

  std::vector< Teuchos::RCP<SmootherT> >  smoothers_;

public:

  MGSmoothers(
    const Teuchos::RCP<const MGOperatorsT>& mgOperators,
    Teuchos::RCP<Teuchos::ParameterList> pl ):
    mgSpaces_( mgOperators->getMGSpaces() ),
    smoothers_( mgSpaces_->getNGrids()-1 ) {

    for( int i=0; i<mgSpaces_->getNGrids()-1; ++i )
      if( mgSpaces_->participating(i) )
        smoothers_[i] = Teuchos::rcp( new SmootherT( mgOperators->get(i), pl ) );

    // not working on brutus
    //smoothers_.shrink_to_fit();
  }

//public:

  /// \brief gets ith smoother, similar to python i=-1 is gets you the coarsest space
  const Teuchos::RCP<SmootherT>& get( int i ) const {
    assert( -1!=i );
    if( i<0 )
      return smoothers_[mgSpaces_->getNGrids()+i];
    else
      return smoothers_[i];
  }

  //  void print(  std::ostream& out=std::cout ) const { }

}; // end of class MGSmoothers



/// \relates MGSmoothers
template< template<class> class SmootherT, class MGOperatorsT >
Teuchos::RCP<const MGSmoothers< MGOperatorsT, SmootherT> >
createMGSmoothers(
  const Teuchos::RCP<const MGOperatorsT>& mgOperators,
  Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ) {

  return Teuchos::rcp( new MGSmoothers<MGOperatorsT,SmootherT>( mgOperators, pl ) );
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


#endif // end of #ifndef PIMPACT_MGSMOOTHERS_HPP
