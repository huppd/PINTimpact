#pragma once
#ifndef PIMPACT_SPACE_HPP
#define PIMPACT_SPACE_HPP


#include "Teuchos_Tuple.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_GridSizeGlobal.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_IndexSpace.hpp"

#include "Pimpact_ProcGridSize.hpp"
#include "Pimpact_ProcGrid.hpp"

#include "Pimpact_Domain.hpp"

#include <iostream>


// \defgroup Space Space
///
/// overloaded class managing indexing, grid ...


namespace Pimpact {


/// \brief manages all Space component, one big composition
/// \ingroup Space
template< class Scalar=double, class Ordinal=int, int dimension=3 >
class Space {

public:

  typedef Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<Ordinal> > >  IndexSpaces;

  Space(
      const Teuchos::RCP<const FieldSpace<Ordinal,dimension> >& fieldSpace,
      const Teuchos::RCP<const IndexSpace<Ordinal> >& scalarIS,
      const IndexSpaces& innerIS,
      const IndexSpaces& fullIS,
      const Teuchos::RCP< GridSizeGlobal<Ordinal,dimension> >& gridSizeGlobal,
      const Teuchos::RCP< GridSizeLocal<Ordinal,dimension> >& gridSizeLocal,
      const Teuchos::RCP< ProcGridSize<Ordinal,dimension> >& procGridSize,
      const Teuchos::RCP< ProcGrid<Ordinal,dimension> >& procGrid,
      const Teuchos::RCP< Domain<Scalar> >& domain ):
        fieldSpace_(fieldSpace),
        scalarIS_(scalarIS),
        innerIS_(innerIS),
        fullIS_(fullIS),
        gridSizeGlobal_(gridSizeGlobal),
        gridSizeLocal_(gridSizeLocal),
        procGridSize_(procGridSize),
        procGrid_(procGrid),
        domain_(domain)
  {}

protected:

  Teuchos::RCP<const FieldSpace<Ordinal,dimension> > fieldSpace_;

  Teuchos::RCP<const IndexSpace<Ordinal> > scalarIS_;

  IndexSpaces innerIS_;
  IndexSpaces fullIS_;

  Teuchos::RCP< GridSizeGlobal<Ordinal,dimension> > gridSizeGlobal_;
  Teuchos::RCP< GridSizeLocal<Ordinal,dimension> > gridSizeLocal_;

  Teuchos::RCP< ProcGridSize<Ordinal,dimension> > procGridSize_;

  Teuchos::RCP< ProcGrid<Ordinal,dimension> > procGrid_;

  Teuchos::RCP< Domain<Scalar> > domain_;

public:


  const MPI_Fint& commf() const { return( procGrid_->commSpacef_ ); }
  const MPI_Comm& comm()  const { return( procGrid_->commSpace_  ); }

  const MPI_Comm& commST()  const { return( procGrid_->commSpaceTime_  ); }

  int rankST() const { return( procGrid_->rankST_ ); }

  const int&      dim()   const { return( fieldSpace_->dim_   ); }

  const Ordinal* nGlo()  const { return( gridSizeGlobal_->getPtr() ); }
  const Ordinal* nLoc()  const { return( gridSizeLocal_->getPtr() ); }
  const Ordinal* bl  ()  const { return( fieldSpace_->bl_.getRawPtr()   ); }
  const Ordinal* bu  ()  const { return( fieldSpace_->bu_.getRawPtr()   ); }


  const Ordinal* sInd() const { return( scalarIS_->sInd_.getRawPtr()  ); }
  const Ordinal* eInd() const { return( scalarIS_->eInd_.getRawPtr()  ); }

  const Ordinal* sInd(  int fieldType ) const { return( innerIS_[fieldType]->sInd_.getRawPtr() ); }
  const Ordinal* eInd(  int fieldType ) const { return( innerIS_[fieldType]->eInd_.getRawPtr() ); }

  const Ordinal* sIndB( int fieldType ) const { return( fullIS_[fieldType]->sInd_.getRawPtr()  ); }
  const Ordinal* eIndB( int fieldType ) const { return( fullIS_[fieldType]->eInd_.getRawPtr()  ); }

  const Ordinal* procCoordinate() const { return( procGrid_->iB_.getRawPtr()  ); }

  const Ordinal* shift() const { return( procGrid_->shift_.getRawPtr()  ); }

  const Ordinal& getNProc(int i) const { return( procGridSize_->get(i) ); }

  Teuchos::RCP<const ProcGrid<Ordinal,dimension> > getProcGrid() const { return( procGrid_ ); }

  Teuchos::RCP<const Domain<Scalar> > getDomain() const { return( domain_ ); }

  void print(  std::ostream& out=std::cout ) const {

    out << "\t---Space: ---\n";

    if( !fieldSpace_.is_null() )
      fieldSpace_->print( out );
    else
      out << "fieldSpace_ is null\n";
    //    MPI_Barrier( commST() );

    if( !gridSizeGlobal_.is_null() ) {
      out <<"---GridSizeGlobal: ---\n";
      gridSizeGlobal_->print( out );
    }
    else
      out << "gridSizeGlobal_ is null\n";
    //    MPI_Barrier( commST() );

    if( !gridSizeLocal_.is_null() ) {
      out <<"---GridSizeLocal: ---\n";
      gridSizeLocal_->print( out );
    }
    else
      out << "gridSizeLocal_ is null\n";
    //    MPI_Barrier( commST() );

    if( !scalarIS_.is_null() )
      scalarIS_->print(out);
    else
      out << "scalarIS_ is null\n";
    //    MPI_Barrier( commST() );

    for( int i=0; i<2; ++i ) {
      innerIS_[i]->print( out );
      fullIS_[i]->print( out );
    }
    //    MPI_Barrier( comm() );

    if( !procGridSize_.is_null() )
      procGridSize_->print( out );
    else
      out << "procGridSize_ is null\n";
    //    MPI_Barrier( commST() );

    if( !procGrid_.is_null() )
      procGrid_->print( out );
    else
      out << "procGrid_ is null\n";
    //    MPI_Barrier( commST() );

  }

  static Teuchos::RCP<const Teuchos::ParameterList>  getValidParameters()  {
    typedef Scalar S;
    typedef Ordinal O;

    static Teuchos::RCP<const Teuchos::ParameterList> validPL;
    // Set all the valid parameters and their default values.
    if(is_null(validPL)) {
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList("Space");
      pl->set("Re", 1., "Re");
      pl->set("alpha2", 1.,
          "\alpha^2");
      // domain type
      pl->set( "domain", 2,
          "Domain type: 0:all dirichlet, 1:dirichlet 2d channel, 2: periodic 2d channel" );

      // domain size
      S l1 = 1.;
      pl->set( "lx", l1, "length in x-direction" );

      S l2 = 1.;
      pl->set("ly", l2, "length in y-direction" );

      S l3 = 1.;
      pl->set("lz", l3, "length in z-direction" );

      int dim = 2;
      pl->set("dim", dim, "dimension of problem" );

      // grid size
      O n1 = 33;
      pl->set("nx", n1, "amount of grid points in x-direction: a*2**q+1" );

      O n2 = 33;
      pl->set("ny", n2, "amount of grid points in y-direction: a*2**q+1" );

      O n3 = 2.;
      pl->set("nz", n3, "amount of grid points in z-direction: a*2**q+1" );

      O nf = 4.;
      pl->set("nf", nf, "amount of grid points in f-direction" );

      O nfs = 1.;
      pl->set("nfs", nfs, "start amount of grid points in f-direction" );

      O nfe = 1.;
      pl->set("nfe", nfe, "end amount of grid points in f-direction" );

      // processor grid size
      O np1 = 2;
      pl->set("npx", np1, "amount of processors in x-direction" );

      O np2 = 2;
      pl->set("npy", np2, "amount of processors in y-direction" );

      O np3 = 1.;
      pl->set("npz", np3, "amount of processors in z-direction" );
      validPL = pl;

      O npf = 1.;
      pl->set("npf", npf, "amount of processors in f-direction" );
      validPL = pl;
    }
    return( validPL );
  }

}; // end of class Space


/// \relates Space
/// \todo wünschenswert initialization from parameterlist
template<class S=double, class O=int, int d=3>
Teuchos::RCP< const Space<S,O,d> > createSpace(
    const Teuchos::RCP<const FieldSpace<O,d> >& fieldSpace,
    const Teuchos::RCP<const IndexSpace<O> >& scalarIS,
    const Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<O> > >& innerIS,
    const Teuchos::ArrayRCP< Teuchos::RCP<const IndexSpace<O> > >& fullIS,
    const Teuchos::RCP< GridSizeGlobal<O,d> >& gridSizeGlobal,
    const Teuchos::RCP< GridSizeLocal<O,d> >& gridSizeLocal,
    const Teuchos::RCP< ProcGridSize<O,d> >& procGridSize,
    const Teuchos::RCP< ProcGrid<O,d> >& procGrid,
    const Teuchos::RCP< Domain<S> >& domain=Teuchos::null ) {

  return(
      Teuchos::rcp(
          new Space<S,O,d>(
              fieldSpace,
              scalarIS,
              innerIS,
              fullIS,
              gridSizeGlobal,
              gridSizeLocal,
              procGridSize,
              procGrid,
              domain ) ) );
}

/// \relates Space
/// \todo wünschenswert initialization from parameterlist
template<class S=double, class O=int, int d=3>
Teuchos::RCP< const Space<S,O,d> > createSpace(
    const Teuchos::RCP<Teuchos::ParameterList> pl ) {

//  int rank =
      Pimpact::init_impact_pre();

  pl->validateParametersAndSetDefaults( *Space<S,O,d>::getValidParameters() );

  Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );

  auto domainSize = Pimpact::createDomainSize<S>( pl->get("Re",1.), pl->get("alpha2",1.), pl->get("lx",1.), pl->get("ly",1.), pl->get("lz",1.) );
  auto boundaryConditionsGlobal = Pimpact::createBoudaryConditionsGlobal( Pimpact::EDomainType( pl->get("domain",2) ) );
  auto procGridSize = Pimpact::createProcGridSize<O,d>( pl->get("npx",2), pl->get("npy",2), pl->get("npz",1), pl->get("npf",1) );
  auto gridSizeGlobal = Pimpact::createGridSizeGlobal<O,d>( pl->get("nx",33), pl->get("ny",33), pl->get("nz",33) );
  auto gridSizeLocal = Pimpact::createGridSizeLocal<O,d>( gridSizeGlobal, procGridSize );

  Pimpact::init_impact_mid();

  auto procGrid = Pimpact::createProcGrid<O,d>( gridSizeLocal, boundaryConditionsGlobal, procGridSize );

  auto boundaryConditionsLocal = Pimpact::createBoudaryConditionsLocal( boundaryConditionsGlobal, procGridSize, procGrid );

  auto fieldSpace = Pimpact::createFieldSpace<O,d>();

  auto scalarIndexSpace = Pimpact::createScalarIndexSpace<O,d>( fieldSpace, gridSizeLocal, boundaryConditionsLocal );

  auto innerIndexSpace = Pimpact::createInnerFieldIndexSpaces<O,d>( fieldSpace, gridSizeLocal, boundaryConditionsLocal );
  auto  fullIndexSpace = Pimpact::createFullFieldIndexSpaces<O,d>(  fieldSpace, gridSizeLocal, boundaryConditionsLocal );

  Pimpact::init_impact_postpost();

//  auto fullIndexSpace = Pimpact::createFullFieldIndexSpaces<O>();

  auto domain = Pimpact::createDomain<S>( domainSize, boundaryConditionsGlobal, boundaryConditionsLocal );

  return( Pimpact::createSpace<S,O>(
      fieldSpace,
      scalarIndexSpace,
      innerIndexSpace,
      fullIndexSpace,
      gridSizeGlobal,
      gridSizeLocal,
      procGridSize,
      procGrid,
      domain) );


}


/// \relates Space
template< class S=double, class O=int, int d=3>
Teuchos::RCP< const Space<S,O,d> > createSpace() {

  return(
      Teuchos::rcp(
          new Space<S,O,d>(
              createFieldSpace<O>(),
              createScalarIndexSpace<O>(),
              createInnerFieldIndexSpaces<O>(),
              createFullFieldIndexSpaces<O>(),
              createGridSizeGlobal<O>(),
              createGridSizeLocal<O,3>(),
              createProcGridSize<O>(),
              createProcGrid<O,d>(),
              createDomain<S>() ) ) );
}

} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SPACE_HPP
