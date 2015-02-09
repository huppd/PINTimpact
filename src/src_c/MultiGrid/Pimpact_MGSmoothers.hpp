#pragma once
#ifndef PIMPACT_MGSMOOTHERS_HPP
#define PIMPACT_MGSMOOTHERS_HPP

#include <vector>

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
//  template<class MGOperatorsT>
  MGSmoothers(
      const Teuchos::RCP<const MGOperatorsT>& mgOperators,
      Teuchos::RCP<Teuchos::ParameterList> pl ):
    mgSpaces_( mgOperators->getMGSpaces() ),
    smoothers_( mgSpaces_->getNGrids() ) {

    for( int i=0; i<mgSpaces_->getNGrids(); ++i )
      smoothers_[i] = Teuchos::rcp( new SmootherT( mgOperators->get(i), pl ) );
//      smoothers_[i] = create<SmootherT>( mgOperators->get(i), pl );

    smoothers_.shrink_to_fit();

  }

//public:

  /// \brief gets ith smoother, similar to python i=-1 is gets you the coarses space
  Teuchos::RCP<SmootherT> get( int i ) const {
    if( i<0 )
      return( smoothers_[mgSpaces_->getNGrids()+i] );
  else
      return( smoothers_[i] );
  }

  //  void print(  std::ostream& out=std::cout ) const {
  //
  //  }

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

#endif // end of #ifndef PIMPACT_MGSMOOTHERS_HPP
