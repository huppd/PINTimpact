#pragma once
#ifndef PIMPACT_SPACEFACTORY_HPP
#define PIMPACT_SPACEFACTORY_HPP


//#include "Teuchos_Tuple.hpp"
//#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"
//#include "Teuchos_ParameterList.hpp"
//#include "Teuchos_XMLParameterListCoreHelpers.hpp"




namespace Pimpact {



/// \brief creates non const template from const
template< template<class> class SpaceObjectT, class SpaceT>
Teuchos::RCP< SpaceObjectT<SpaceT> >
create(
		const Teuchos::RCP<const SpaceT>& space ) {
  return(
			Teuchos::rcp( new SpaceObjectT<SpaceT>(space) ) );
}


/// \brief creates non const template from non const
template< template<class> class SpaceObjectT, class SpaceT>
Teuchos::RCP< SpaceObjectT<SpaceT> >
create( const Teuchos::RCP< SpaceT>& space ) {
  return(
			Teuchos::rcp( new SpaceObjectT<SpaceT>(space) ) );
}


/// \brief creates non const from non const
template< class SpaceObjectT>
Teuchos::RCP<SpaceObjectT>
create( const Teuchos::RCP<const typename SpaceObjectT::SpaceT>& space ) {
  return(
			Teuchos::rcp( new SpaceObjectT( space ) ) );
}


/// \brief creates non const from non const
template< class SpaceObjectT, class SpaceT>
Teuchos::RCP< SpaceObjectT >
create( const Teuchos::RCP< SpaceT>& space ) {
	return( Teuchos::rcp( new SpaceObjectT( space ) ) );
}


template< template<class> class SpaceObjectT, class SpaceT>
Teuchos::RCP<const SpaceObjectT<SpaceT> >
createConst( const Teuchos::RCP<const SpaceT>& space ) {
	return( Teuchos::rcp( new SpaceObjectT<SpaceT>(space) ) );
}


template< template<class> class SpaceObjectT, class SpaceT>
Teuchos::RCP<const SpaceObjectT<SpaceT> >
createConst( const Teuchos::RCP< SpaceT>& space ) {
	return( Teuchos::rcp( new SpaceObjectT<SpaceT>(space) ) );
}


template< class SpaceObjectT, class SpaceT>
Teuchos::RCP<const SpaceObjectT >
createConst( const Teuchos::RCP<const SpaceT>& space ) {
	return( Teuchos::rcp( new SpaceObjectT( space ) ) );
}


template< class SpaceObjectT, class SpaceT>
Teuchos::RCP<const SpaceObjectT >
createConst( const Teuchos::RCP< SpaceT>& space ) {
	return( Teuchos::rcp( new SpaceObjectT( space ) ) );
}


/// \relates TransferOp
template<class OpT>
Teuchos::RCP<const OpT >
create(
    const Teuchos::RCP<const typename OpT::FSpaceT>& fSpace,
    const Teuchos::RCP<const typename OpT::CSpaceT>& cSpace ) {

	return( Teuchos::rcp( new OpT( fSpace, cSpace ) ) );
}




} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SPACEFACTORY_HPP
