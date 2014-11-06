#pragma once
#ifndef PIMPACT_SPACE_HPP
#define PIMPACT_SPACE_HPP


#include "Teuchos_Tuple.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

#include "Pimpact_StencilWidths.hpp"
#include "Pimpact_GridSizeGlobal.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_IndexSpace.hpp"

#include "Pimpact_ProcGridSize.hpp"
#include "Pimpact_ProcGrid.hpp"

#include "Pimpact_Domain.hpp"

#include "Pimpact_GridCoordinatesGlobal.hpp"
#include "Pimpact_GridCoordinatesLocal.hpp"

#include "Pimpact_InterpolateV2SOp.hpp"

#include "pimpact.hpp"

#include <iostream>


/// \defgroup Space Space
///
/// overloaded class managing indexing, grid ...


namespace Pimpact {


extern "C" {
  void openH5F();
  void closeH5F();
}
/// \brief Space in the sense of a VectorSpace, it is the connection between Field and Operators
///
/// \ingroup Space
template< class Scalar=double, class Ordinal=int, int dimension=3 >
class Space {

public:

  Space(
      const Teuchos::RCP<const StencilWidths<dimension> >& stencilWidths,
      const Teuchos::RCP<const IndexSpace<Ordinal> >& indexSpace,
      const Teuchos::RCP<const GridSizeGlobal<Ordinal,dimension> >& gridSizeGlobal,
      const Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> >& gridSizeLocal,
      const Teuchos::RCP<const ProcGridSize<Ordinal,dimension> >& procGridSize,
      const Teuchos::RCP<const ProcGrid<Ordinal,dimension> >& procGrid,
      const Teuchos::RCP<const GridCoordinatesGlobal<Scalar,Ordinal,dimension> >& coordGlobal,
      const Teuchos::RCP<const GridCoordinatesLocal<Scalar,Ordinal,dimension> >& coordLocal,
      const Teuchos::RCP<const Domain<Scalar> >& domain,
      const Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension> >& interV2S ):
        stencilWidths_(stencilWidths),
        indexSpace_(indexSpace),
        gridSizeGlobal_(gridSizeGlobal),
        gridSizeLocal_(gridSizeLocal),
        procGridSize_(procGridSize),
        procGrid_(procGrid),
        coordGlobal_(coordGlobal),
        coordLocal_(coordLocal),
        domain_(domain),
        interV2S_(interV2S)
  {
    openH5F();
  }

//  ~Space(){ closeH5F(); }

protected:

  Teuchos::RCP<const StencilWidths<dimension> > stencilWidths_;

  Teuchos::RCP<const IndexSpace<Ordinal> > indexSpace_;

  Teuchos::RCP<const GridSizeGlobal<Ordinal,dimension> > gridSizeGlobal_;

  Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> > gridSizeLocal_;

  Teuchos::RCP<const ProcGridSize<Ordinal,dimension> > procGridSize_;

  Teuchos::RCP<const ProcGrid<Ordinal,dimension> > procGrid_;

  Teuchos::RCP<const GridCoordinatesGlobal<Scalar,Ordinal,dimension> > coordGlobal_;

  Teuchos::RCP<const GridCoordinatesLocal<Scalar,Ordinal,dimension> > coordLocal_;

  Teuchos::RCP<const Domain<Scalar> > domain_;

  Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension> > interV2S_;

public:

  Teuchos::RCP<const StencilWidths<dimension> > getStencilWidths() const { return( stencilWidths_ ); }

  Teuchos::RCP<const IndexSpace<Ordinal> > getIndexSpace() const { return( indexSpace_ ); }

  Teuchos::RCP<const GridSizeGlobal<Ordinal,dimension> > getGridSizeGlobal() const { return( gridSizeGlobal_ );  }

  Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> > getGridSizeLocal() const { return( gridSizeLocal_ );  }

  Teuchos::RCP<const ProcGridSize<Ordinal,dimension> > getProcGridSize() const { return( procGridSize_ ); }

  Teuchos::RCP<const ProcGrid<Ordinal,dimension> > getProcGrid() const { return( procGrid_ ); }

  Teuchos::RCP<const GridCoordinatesGlobal<Scalar,Ordinal,dimension> > getCoordinatesGlobal() const {
    return( coordGlobal_ );
  }

  Teuchos::RCP<const GridCoordinatesLocal<Scalar,Ordinal,dimension> > getCoordinatesLocal() const {
    return( coordLocal_ );
  }

  Teuchos::RCP<const Domain<Scalar> > getDomain() const { return( domain_ ); }

  Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension> > getInterpolateV2S() const { return( interV2S_ ); }


  const MPI_Fint& commf() const { return( procGrid_->commSpacef_ ); }
  const MPI_Comm& comm()  const { return( procGrid_->commSpace_  ); }

  const MPI_Comm& commST()  const { return( procGrid_->commSpaceTime_  ); }

  int rankST() const { return( procGrid_->getRank() ); }
//  int rank  () const { return( procGrid_->rank_   ); }

  const int&      dim()   const { return( domain_->getDomainSize()->getDim() ); }

  const Ordinal* nGlo()        const { return( gridSizeGlobal_->get()  ); }
  const Ordinal& nGlo( int i ) const { return( gridSizeGlobal_->get(i) ); }

  const Ordinal* nLoc()        const { return( gridSizeLocal_->get()  ); }
  const Ordinal& nLoc( int i ) const { return( gridSizeLocal_->get(i) ); }

  const Ordinal* bl()         const { return( stencilWidths_->getBL()   ); }
  const Ordinal& bl( int i )  const { return( stencilWidths_->getBL(i)  ); }

  const Ordinal* bu()         const { return( stencilWidths_->getBU()   ); }
  const Ordinal& bu( int i )  const { return( stencilWidths_->getBU(i)  ); }

  const Ordinal* dl()         const { return( stencilWidths_->getDL()   ); }
  const Ordinal& dl( int i )  const { return( stencilWidths_->getDL(i)  ); }

  const Ordinal* du()         const { return( stencilWidths_->getDU()   ); }
  const Ordinal& du( int i )  const { return( stencilWidths_->getDU(i)  ); }

  const Ordinal* gl()         const { return( stencilWidths_->getGL()   ); }
  const Ordinal& gl( int i )  const { return( stencilWidths_->getGL(i)  ); }

  const Ordinal* gu()         const { return( stencilWidths_->getGU()   ); }
  const Ordinal& gu( int i )  const { return( stencilWidths_->getGU(i)  ); }

  const Ordinal* nl()         const { return( stencilWidths_->getNL()   ); }
  const Ordinal& nl( int i )  const { return( stencilWidths_->getNL(i)  ); }

  const Ordinal* nu()         const { return( stencilWidths_->getNU()   ); }
  const Ordinal& nu( int i )  const { return( stencilWidths_->getNU(i)  ); }


  const Ordinal* sInd( int fieldType ) const {
    return( indexSpace_->sInd( fieldType ) );
  }
  const Ordinal* eInd(  int fieldType ) const {
    return( indexSpace_->eInd( fieldType ) );
  }

  const Ordinal* sIndB( int fieldType ) const {
    return( indexSpace_->sIndB( fieldType ) );
  }
  const Ordinal* eIndB( int fieldType ) const {
    return( indexSpace_->eIndB( fieldType ) );
  }

  const Ordinal& sInd( int fieldType, int dir ) const {
    return( indexSpace_->sInd( fieldType, dir ) );
  }
  const Ordinal& eInd(  int fieldType, int dir ) const {
    return( indexSpace_->eInd( fieldType, dir ) );
  }

  const Ordinal& sIndB( int fieldType, int dir ) const {
    return( indexSpace_->sIndB( fieldType, dir ) );
  }
  const Ordinal& eIndB( int fieldType, int dir ) const {
    return( indexSpace_->eIndB( fieldType, dir ) );
  }

  const Ordinal* procCoordinate() const { return( procGrid_->getIB()  ); }

  const Ordinal* getShift()      const { return( procGrid_->shift_.getRawPtr()  ); }
  const Ordinal& getShift(int i) const { return( procGrid_->shift_[i]  ); }

//  const Ordinal* getNProc()      const { return( procGridSize_->get() ); }
  const Ordinal& getNProc(int i) const { return( procGridSize_->get(i) ); }

  void print(  std::ostream& out=std::cout ) const {

      out << "\t---Space: ---\n";

      stencilWidths_->print( out );

      out <<"---GridSizeGlobal: ---\n";
      gridSizeGlobal_->print( out );

      out <<"---GridSizeLocal: ---\n";
      gridSizeLocal_->print( out );

      indexSpace_->print(out);

      procGridSize_->print( out );

      procGrid_->print( out );

      coordGlobal_->print(out);

      coordLocal_->print(out);

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
      int dim = 2;
      pl->set("dim", dim, "dimension of problem" );

      S l1 = 1.;
      pl->set( "lx", l1, "length in x-direction" );

      S l2 = 1.;
      pl->set("ly", l2, "length in y-direction" );

      S l3 = 1.;
      pl->set("lz", l3, "length in z-direction" );


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
/// \deprecated \param setImpact should be uneccessary in the future
template<class S=double, class O=int, int d=3>
Teuchos::RCP<const Space<S,O,d> > createSpace(
    Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::null,
    bool setImpact=false ) {

  if( setImpact ) Pimpact::init_impact_pre();
  if( pl.is_null() ) pl = Teuchos::parameterList();

  pl->validateParametersAndSetDefaults( *Space<S,O,d>::getValidParameters() );

  Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );

  auto domainSize = Pimpact::createDomainSize<S>(
      pl->get("dim",2),
      pl->get("Re",1.),
      pl->get("alpha2",1.),
      pl->get("lx",2.),
      pl->get("ly",2.),
      pl->get("lz",1.) );
  if( setImpact ) domainSize->set_Impact();

  auto boundaryConditionsGlobal = Pimpact::createBoudaryConditionsGlobal( Pimpact::EDomainType( pl->get("domain",2) ) );
  if( setImpact ) boundaryConditionsGlobal->set_Impact();

  auto procGridSize = Pimpact::createProcGridSize<O,d>( pl->get("npx",2), pl->get("npy",2), pl->get("npz",1), pl->get("npf",1) );
  if( setImpact ) procGridSize->set_Impact();

  auto gridSizeGlobal = Pimpact::createGridSizeGlobal<O,d>( pl->get("nx",33), pl->get("ny",33), pl->get("nz",2), pl->get("nf",32) );
  if( setImpact ) gridSizeGlobal->set_Impact();

  auto gridSizeLocal = Pimpact::createGridSizeLocal<O,d>( gridSizeGlobal, procGridSize );
  if( setImpact ) gridSizeLocal->set_Impact();

  if( setImpact ) Pimpact::init_impact_mid();

  auto procGrid = Pimpact::createProcGrid<O,d>( gridSizeLocal, boundaryConditionsGlobal, procGridSize );
  if( setImpact ) procGrid->set_Impact();

  auto boundaryConditionsLocal = Pimpact::createBoudaryConditionsLocal( boundaryConditionsGlobal, procGridSize, procGrid );
  if( setImpact ) boundaryConditionsLocal->set_Impact();

  auto fieldSpace = Pimpact::createStencilWidths<d>();

  auto indexSpace = Pimpact::createIndexSpace<O,d>( fieldSpace, gridSizeLocal, boundaryConditionsLocal, setImpact );

  if( setImpact ) Pimpact::init_impact_postpost();

  auto domain =
      Pimpact::createDomain<S>( domainSize, boundaryConditionsGlobal, boundaryConditionsLocal );

  auto  coordGlobal =
      Pimpact::createGridCoordinatesGlobal<S,O,d>( gridSizeGlobal, domainSize );

  auto  coordLocal =
      Pimpact::createGridCoordinatesLocal<S,O,d>(
          fieldSpace,
          domainSize,
          gridSizeGlobal,
          gridSizeLocal,
          boundaryConditionsGlobal,
          boundaryConditionsLocal,
          procGrid,
          coordGlobal );

  auto interV2S =
      Pimpact::createInterpolateV2S<S,O,d>(
          procGrid,
          gridSizeLocal,
          fieldSpace,
          domain,
          coordLocal );

  return(
       Teuchos::rcp(
           new Space<S,O,d>(
               fieldSpace,
               indexSpace,
               gridSizeGlobal,
               gridSizeLocal,
               procGridSize,
               procGrid,
               coordGlobal,
               coordLocal,
               domain,
               interV2S ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SPACE_HPP
