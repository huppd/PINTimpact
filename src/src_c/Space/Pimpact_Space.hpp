#pragma once
#ifndef PIMPACT_SPACE_HPP
#define PIMPACT_SPACE_HPP


#include <iostream>

#include "Teuchos_Tuple.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"

#include "pimpact.hpp"

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


#include "Pimpact_SpaceFactory.hpp"



/// \defgroup SpaceObject Space Objects
///
/// overloaded class managing indexing, grid ...



namespace Pimpact {

extern "C" {
void openH5F();
void closeH5F();
}



/// \brief Space in the sense of a vector space, it is the connection between Field and Operators
///
/// \ingroup SpaceObject
template<class S=double, class O=int, int d=3, int dNC=4>
class Space {

public:

  typedef S Scalar;
  typedef O Ordinal;

  static const int dimension = d;

  static const int dimNC = dNC;


	Space( Teuchos::RCP<Teuchos::ParameterList> pl ) {

    pl->validateParametersAndSetDefaults( *getValidParameters() );

    Teuchos::writeParameterListToXmlFile( *pl, "parameterSpace.xml" );

    stencilWidths_ = Pimpact::createStencilWidths<dimension,dimNC>();

    auto domainSize = Pimpact::createDomainSize<S>(
        pl->get<int>("dim"),
        pl->get<S>("Re"),
        pl->get<S>("alpha2"),
        pl->get<S>("lx"),
        pl->get<S>("ly"),
        pl->get<S>("lz") );

    // are all template paramter needed here?
		int domain = pl->get<int>("domain");
		domain = ( 2==pl->get<int>("dim") && 0==domain )?1:domain;
    auto boundaryConditionsGlobal =
        Pimpact::createBoudaryConditionsGlobal<d>( Pimpact::EDomainType( domain ) );

    procGridSize_ =
        Pimpact::createProcGridSize<O,d>( pl->get<O>("npx"), pl->get<O>("npy",2), pl->get<O>("npz",1), pl->get<O>("npf") );

    gridSizeGlobal_ =
        Pimpact::createGridSizeGlobal<O,d>(
            pl->get<O>("nx"),
            pl->get<O>("ny"),
            ( 2==pl->get<int>("dim") )?2:pl->get<O>("nz"),
            pl->get<O>("nf") );

    gridSizeLocal_ =
        Pimpact::createGridSizeLocal<O,d,dNC>( gridSizeGlobal_, procGridSize_, stencilWidths_ );


    procGrid_ =
        Pimpact::createProcGrid<O,d>(
            gridSizeLocal_,
            boundaryConditionsGlobal,
            procGridSize_ );

    auto boundaryConditionsLocal =
        Pimpact::createBoudaryConditionsLocal(
            boundaryConditionsGlobal,
            procGridSize_,
            procGrid_ );

    indexSpace_ =
        Pimpact::createIndexSpace<O,d>(
            stencilWidths_,
            gridSizeLocal_,
            boundaryConditionsLocal );


    domain_ =
        Pimpact::createDomain<S,d>(
            domainSize,
            boundaryConditionsGlobal,
            boundaryConditionsLocal );

    coordGlobal_ =
        Pimpact::createGridCoordinatesGlobal<S,O,d>( gridSizeGlobal_, domainSize );

    coordLocal_ =
        Pimpact::createGridCoordinatesLocal<S,O,d>(
            stencilWidths_,
            domainSize,
            gridSizeGlobal_,
            gridSizeLocal_,
            boundaryConditionsGlobal,
            boundaryConditionsLocal,
            procGrid_,
            coordGlobal_ );

    interV2S_ =
        Pimpact::createInterpolateV2S<Scalar,Ordinal,dimension,dimNC>(
            procGrid_,
            gridSizeLocal_,
            stencilWidths_,
            domain_,
            coordLocal_ );

    openH5F();

  }


  Space(
      const Teuchos::RCP<const StencilWidths<dimension,dimNC> >& stencilWidths,
      const Teuchos::RCP<const IndexSpace<Ordinal,dimension> >& indexSpace,
      const Teuchos::RCP<const GridSizeGlobal<Ordinal,dimension> >& gridSizeGlobal,
      const Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> >& gridSizeLocal,
      const Teuchos::RCP<const ProcGridSize<Ordinal,dimension> >& procGridSize,
      const Teuchos::RCP<const ProcGrid<Ordinal,dimension> >& procGrid,
      const Teuchos::RCP<const GridCoordinatesGlobal<Scalar,Ordinal,dimension> >& coordGlobal,
      const Teuchos::RCP<const GridCoordinatesLocal<Scalar,Ordinal,dimension> >& coordLocal,
      const Teuchos::RCP<const Domain<Scalar,dimension> >& domain,
      const Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension,dimNC> >& interV2S ):
        stencilWidths_(stencilWidths),
        indexSpace_(indexSpace),
        gridSizeGlobal_(gridSizeGlobal),
        gridSizeLocal_(gridSizeLocal),
        procGridSize_(procGridSize),
        procGrid_(procGrid),
        coordGlobal_(coordGlobal),
        coordLocal_(coordLocal),
        domain_(domain),
        interV2S_(interV2S) {
    //    openH5F();
  }

  //  ~Space(){ closeH5F(); }

protected:

  Teuchos::RCP<const StencilWidths<dimension,dimNC> > stencilWidths_;

  Teuchos::RCP<const IndexSpace<Ordinal,dimension> > indexSpace_;

  Teuchos::RCP<const GridSizeGlobal<Ordinal,dimension> > gridSizeGlobal_;

  Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> > gridSizeLocal_;

  Teuchos::RCP<const ProcGridSize<Ordinal,dimension> > procGridSize_;

  Teuchos::RCP<const ProcGrid<Ordinal,dimension> > procGrid_;

  Teuchos::RCP<const GridCoordinatesGlobal<Scalar,Ordinal,dimension> > coordGlobal_;

  Teuchos::RCP<const GridCoordinatesLocal<Scalar,Ordinal,dimension> > coordLocal_;

  Teuchos::RCP<const Domain<Scalar,dimension> > domain_;

  Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension,dimNC> > interV2S_;

public:

  Teuchos::RCP<const StencilWidths<dimension,dimNC> > getStencilWidths() const { return( stencilWidths_ ); }

  Teuchos::RCP<const IndexSpace<Ordinal,dimension> > getIndexSpace() const { return( indexSpace_ ); }

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

  Teuchos::RCP<const Domain<Scalar,dimension> > getDomain() const { return( domain_ ); }

  Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension,dimNC> > getInterpolateV2S() const { return( interV2S_ ); }


  const MPI_Comm& comm()  const { return( procGrid_->getCommS()  ); }

//  const MPI_Comm& commST()  const { return( procGrid_->getCommWorld()  ); }

  int rankST() const { return( procGrid_->getRank() ); }
  int rankS () const { return( procGrid_->getRankS() ); }

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

  const Ordinal* getShift()      const { return( procGrid_->getShift()  ); }
  const Ordinal& getShift(int i) const { return( procGrid_->getShift(i)  ); }

	const Ordinal* getNProc()      const { return( procGridSize_->get() ); }
  const Ordinal& getNProc(int i) const { return( procGridSize_->get(i) ); }

  void print(  std::ostream& out=std::cout ) const {

    out << "\t---Space: ---\n";

    stencilWidths_->print( out );

    out <<"\t---GridSizeGlobal: ---\n";
    gridSizeGlobal_->print( out );

    out <<"\t---GridSizeLocal: ---\n";
    gridSizeLocal_->print( out );

    indexSpace_->print(out);

    procGridSize_->print( out );

    procGrid_->print( out );

		getDomain()->getBCLocal()->print( out );

//    coordGlobal_->print(out);
//
//		coordLocal_->print(out);

  }

  static Teuchos::RCP<const Teuchos::ParameterList>  getValidParameters()  {

    static Teuchos::RCP<const Teuchos::ParameterList> validPL;

    // Set all the valid parameters and their default values.
    if(is_null(validPL)) {

      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList("Space");
      pl->set<S>("Re", 1., "Reynolds number");
      pl->set<S>("alpha2", 1.,
          "Womersley square alpha^2");
      // domain type
      pl->set<int>( "domain", 2,
          "Domain type: 0:all dirichlet, 1:dirichlet 2d channel, 2: periodic 2d channel" );

      // domain size
      pl->set<int>("dim", 3, "dimension of problem" );

      pl->set<S>( "lx", 1., "length in x-direction" );
      pl->set<S>( "ly", 1., "length in y-direction" );
      pl->set<S>( "lz", 1., "length in z-direction" );


      // grid size
      pl->set<O>("nx", 33, "amount of grid points in x-direction: a*2**q+1" );
      pl->set<O>("ny", 33, "amount of grid points in y-direction: a*2**q+1" );
      pl->set<O>("nz", 33, "amount of grid points in z-direction: a*2**q+1" );
      pl->set<O>("nf", 4, "amount of grid points in f-direction" );

      pl->set<O>("nfs", 1, "start amount of grid points in f-direction" );
      pl->set<O>("nfe", 1, "end amount of grid points in f-direction" );

      // processor grid size
      pl->set<O>("npx", 2, "amount of processors in x-direction" );
      pl->set<O>("npy", 2, "amount of processors in y-direction" );
      pl->set<O>("npz", 1, "amount of processors in z-direction" );
      pl->set<O>("npf", 1, "amount of processors in f-direction" );

      validPL = pl;
    }
    return( validPL );
  }


//  /// \todo make something like this work
//  static Teuchos::RCP<Teuchos::ParameterList>  getParameters( Teuchos::CommandLineProcessor& clp )  {
//
//    // processor grid size
//    O npx = 2;
//    clp.setOption( "npx", &npx, "amount of processors in x-direction" );
//
//    O npy = 2;
//    clp.setOption( "npy", &npy, "amount of processors in y-direction" );
//
//    O npz = 1.;
//    clp.setOption( "npz", &npz, "amount of processors in z-direction" );
//
//
//      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList("Space");
////      pl->set<S>("Re", 1., "Reynolds number");
////      pl->set<S>("alpha2", 1.,
////          "Womersley square \alpha^2");
////      // domain type
////      pl->set<int>( "domain", 2,
////          "Domain type: 0:all dirichlet, 1:dirichlet 2d channel, 2: periodic 2d channel" );
////
////      // domain size
////      pl->set<int>("dim", 2, "dimension of problem" );
////
////      pl->set<S>( "lx", 1., "length in x-direction" );
////      pl->set<S>( "ly", 1., "length in y-direction" );
////      pl->set<S>( "lz", 1., "length in z-direction" );
//
//
////      // grid size
////      pl->set<O>("nx", 33, "amount of grid points in x-direction: a*2**q+1" );
////      pl->set<O>("ny", 33, "amount of grid points in y-direction: a*2**q+1" );
////      pl->set<O>("nz", 2, "amount of grid points in z-direction: a*2**q+1" );
////      pl->set<O>("nf", 4, "amount of grid points in f-direction" );
////
////      pl->set<O>("nfs", 1, "start amount of grid points in f-direction" );
////      pl->set<O>("nfe", 1, "end amount of grid points in f-direction" );
//
//      // processor grid size
//      pl->set<O>("npx", npx, "amount of processors in x-direction" );
//      pl->set<O>("npy", npy, "amount of processors in y-direction" );
//      pl->set<O>("npz", npz, "amount of processors in z-direction" );
////      pl->set<O>("npf", n, "amount of processors in f-direction" );
//
//    return( pl );
//  }

}; // end of class Space



/// \relates Space
/// \deprecated \c setImpact should be uneccessary in the future
template<class S=double, class O=int, int d=3, int dimNC=4>
Teuchos::RCP<const Space<S,O,d,dimNC> >
createSpace( Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ) {

	return( Teuchos::rcp( new Space<S,O,d,dimNC>( pl ) ) );

}

Teuchos::RCP<std::ostream> createOstream( const std::string& fname, int rank);

#ifdef COMPILE_ETI
extern template class Space<double,int,3,2>;
extern template class Space<double,int,3,4>;
extern template class Space<double,int,4,2>;
extern template class Space<double,int,4,4>;
#endif



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_SPACE_HPP
