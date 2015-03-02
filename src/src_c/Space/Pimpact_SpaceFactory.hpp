#pragma once
#ifndef PIMPACT_SPACEFACTORY_HPP
#define PIMPACT_SPACEFACTORY_HPP


//#include "Teuchos_Tuple.hpp"
//#include "Teuchos_ArrayRCP.hpp"
//#include "Teuchos_RCP.hpp"
//#include "Teuchos_ParameterList.hpp"
//#include "Teuchos_XMLParameterListCoreHelpers.hpp"
//
//#include "Pimpact_StencilWidths.hpp"
//#include "Pimpact_GridSizeGlobal.hpp"
//#include "Pimpact_GridSizeLocal.hpp"
//#include "Pimpact_IndexSpace.hpp"
//
//#include "Pimpact_ProcGridSize.hpp"
//#include "Pimpact_ProcGrid.hpp"
//
//#include "Pimpact_Domain.hpp"
//
//#include "Pimpact_GridCoordinatesGlobal.hpp"
//#include "Pimpact_GridCoordinatesLocal.hpp"
//
//#include "Pimpact_InterpolateV2SOp.hpp"
//
//#include "pimpact.hpp"
//
//#include <iostream>






namespace Pimpact {



template< template<class> class SpaceObjectT, class SpaceT>
Teuchos::RCP< SpaceObjectT<SpaceT> >
create( const Teuchos::RCP<const SpaceT>& space ) {
  return(
      Teuchos::rcp( new SpaceObjectT<SpaceT>(space) )
  );
}


template< template<class> class SpaceObjectT, class SpaceT>
Teuchos::RCP< SpaceObjectT<SpaceT> >
create( const Teuchos::RCP< SpaceT>& space ) {
  return(
      Teuchos::rcp( new SpaceObjectT<SpaceT>(space) )
  );
}


template< class SpaceObjectT, class SpaceT>
Teuchos::RCP< SpaceObjectT >
create( const Teuchos::RCP<const SpaceT>& space ) {
  return(
      Teuchos::rcp( new SpaceObjectT( space ) )
  );
}


template< class SpaceObjectT, class SpaceT>
Teuchos::RCP< SpaceObjectT >
create( const Teuchos::RCP< SpaceT>& space ) {
  return(
      Teuchos::rcp( new SpaceObjectT( space ) )
  );
}


template< template<class> class SpaceObjectT, class SpaceT>
Teuchos::RCP<const SpaceObjectT<SpaceT> >
createConst( const Teuchos::RCP<const SpaceT>& space ) {
  return(
      Teuchos::rcp( new SpaceObjectT<SpaceT>(space) )
  );
}


template< template<class> class SpaceObjectT, class SpaceT>
Teuchos::RCP<const SpaceObjectT<SpaceT> >
createConst( const Teuchos::RCP< SpaceT>& space ) {
  return(
      Teuchos::rcp( new SpaceObjectT<SpaceT>(space) )
  );
}


template< class SpaceObjectT, class SpaceT>
Teuchos::RCP<const SpaceObjectT >
createConst( const Teuchos::RCP<const SpaceT>& space ) {
  return(
      Teuchos::rcp( new SpaceObjectT( space ) )
  );
}


template< class SpaceObjectT, class SpaceT>
Teuchos::RCP<const SpaceObjectT >
createConst( const Teuchos::RCP< SpaceT>& space ) {
  return(
      Teuchos::rcp( new SpaceObjectT( space ) )
  );
}


/// \relates TransferOp
template<class OpT>
Teuchos::RCP<const OpT >
create(
    const Teuchos::RCP<const typename OpT::FSpaceT>& fSpace,
    const Teuchos::RCP<const typename OpT::CSpaceT>& cSpace ) {

  return(
      Teuchos::rcp(
          new OpT( fSpace, cSpace )
      )
  );

}




} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_SPACEFACTORY_HPP
