#pragma once
#ifndef PIMPACT_SPACE_HPP
#define PIMPACT_SPACE_HPP


#include <iostream>

#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"


#include "Pimpact_BoundaryConditionsGlobal.hpp"
#include "Pimpact_BoundaryConditionsLocal.hpp"
#include "Pimpact_DomainSize.hpp"
#include "Pimpact_GridCoordinatesGlobal.hpp"
#include "Pimpact_GridCoordinatesLocal.hpp"
#include "Pimpact_GridSizeGlobal.hpp"
#include "Pimpact_GridSizeLocal.hpp"
#include "Pimpact_IndexSpace.hpp"
#include "Pimpact_InterpolateV2SOp.hpp"
#include "Pimpact_ProcGrid.hpp"
#include "Pimpact_StencilWidths.hpp"
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

		stencilWidths_ =
			Pimpact::createStencilWidths<dimension,dimNC>(
					pl->get<bool>("spectral in time") );


		int dim = pl->get<int>("dim");
		TEUCHOS_TEST_FOR_EXCEPT( 2!=dim && 3!=dim );
		domainSize_ = Pimpact::createDomainSize<S>(
				dim,
				pl->get<S>("Re"),
				pl->get<S>("alpha2"),
				pl->get<S>("lx"),
				pl->get<S>("ly"),
				pl->get<S>("lz") );

		// are all template paramter needed here?
		int domain = pl->get<int>("domain");
		domain = ( 2==pl->get<int>("dim") && 0==domain )?1:domain;

		boundaryConditionsGlobal_ =
			Pimpact::createBoudaryConditionsGlobal<d>(
					Pimpact::EDomainType( domain ) );

		gridSizeGlobal_ =
			Pimpact::createGridSizeGlobal<O,d>(
					pl->get<O>("nx"),
					pl->get<O>("ny"),
					( 2==pl->get<int>("dim") )?2:pl->get<O>("nz"),
					pl->get<O>("nf") );

		Teuchos::Tuple<O,d> procGridSize;
		procGridSize[0] = pl->get<O>("npx");
		procGridSize[1] = pl->get<O>("npy");
		procGridSize[2] = pl->get<O>("npz");
		if( d>3 )
			procGridSize[3] = pl->get<O>("npf");

		procGrid_ =
			Pimpact::createProcGrid<O,d>(
					procGridSize,
					boundaryConditionsGlobal_ );

		boundaryConditionsLocal_ =
			Pimpact::createBoudaryConditionsLocal(
					boundaryConditionsGlobal_,
					procGrid_ );

		gridSizeLocal_ =
			Pimpact::createGridSizeLocal<O,d,dNC>(
					gridSizeGlobal_,
					procGrid_,
					stencilWidths_ );

		indexSpace_ =
			Pimpact::createIndexSpace<O,d>(
					stencilWidths_,
					gridSizeLocal_,
					boundaryConditionsLocal_,
				 	procGrid_ );

		coordGlobal_ =
			Pimpact::createGridCoordinatesGlobal<S,O,d>(
					gridSizeGlobal_,
					domainSize_,
					Teuchos::tuple( None, None, None) );

		coordLocal_ =
			Pimpact::createGridCoordinatesLocal<S,O,d>(
					stencilWidths_,
					domainSize_,
					gridSizeGlobal_,
					gridSizeLocal_,
					boundaryConditionsGlobal_,
					boundaryConditionsLocal_,
					procGrid_,
					coordGlobal_ );

		interV2S_ =
			Pimpact::createInterpolateV2S<Scalar,Ordinal,dimension,dimNC>(
					indexSpace_,
					gridSizeLocal_,
					stencilWidths_,
					domainSize_,
					boundaryConditionsLocal_,
					coordLocal_ );

		openH5F();

	}


	Space(
			const Teuchos::RCP<const StencilWidths<dimension,dimNC> >& stencilWidths,
			const Teuchos::RCP<const IndexSpace<Ordinal,dimension> >& indexSpace,
			const Teuchos::RCP<const GridSizeGlobal<Ordinal,dimension> >& gridSizeGlobal,
			const Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> >& gridSizeLocal,
			const Teuchos::RCP<const ProcGrid<Ordinal,dimension> >& procGrid,
			const Teuchos::RCP<const GridCoordinatesGlobal<Scalar,Ordinal,dimension> >& coordGlobal,
			const Teuchos::RCP<const GridCoordinatesLocal<Scalar,Ordinal,dimension> >& coordLocal,
			const Teuchos::RCP<const DomainSize<Scalar> > domainSize,
			const Teuchos::RCP<const BoundaryConditionsGlobal<dimension> > boundaryConditionsGlobal,
			const Teuchos::RCP<const BoundaryConditionsLocal> boundaryConditionsLocal,
			const Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension,dimNC> >& interV2S ):
		stencilWidths_(stencilWidths),
		indexSpace_(indexSpace),
		gridSizeGlobal_(gridSizeGlobal),
		gridSizeLocal_(gridSizeLocal),
		procGrid_(procGrid),
		coordGlobal_(coordGlobal),
		coordLocal_(coordLocal),
		domainSize_(domainSize),
		boundaryConditionsGlobal_(boundaryConditionsGlobal),
		boundaryConditionsLocal_(boundaryConditionsLocal),
		interV2S_(interV2S) {}

	//  ~Space(){ closeH5F(); }

protected:

	Teuchos::RCP<const StencilWidths<dimension,dimNC> >
		stencilWidths_;

	Teuchos::RCP<const IndexSpace<Ordinal,dimension> >
		indexSpace_;

	Teuchos::RCP<const GridSizeGlobal<Ordinal,dimension> >
		gridSizeGlobal_;

	Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> >
		gridSizeLocal_;

	Teuchos::RCP<const ProcGrid<Ordinal,dimension> >
		procGrid_;

	Teuchos::RCP<const GridCoordinatesGlobal<Scalar,Ordinal,dimension> >
		coordGlobal_;

	Teuchos::RCP<const GridCoordinatesLocal<Scalar,Ordinal,dimension> >
		coordLocal_;

	Teuchos::RCP<const DomainSize<Scalar> >
		domainSize_;

	Teuchos::RCP<const BoundaryConditionsGlobal<dimension> >
		boundaryConditionsGlobal_;

	Teuchos::RCP<const BoundaryConditionsLocal>
		boundaryConditionsLocal_;

	Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension,dimNC> >
		interV2S_;

public:

	/// \name getter methods
	/// \{

	Teuchos::RCP<const StencilWidths<dimension,dimNC> >
		getStencilWidths() const { return( stencilWidths_ ); }

	Teuchos::RCP<const IndexSpace<Ordinal,dimension> >
		getIndexSpace() const { return( indexSpace_ ); }

	Teuchos::RCP<const GridSizeGlobal<Ordinal,dimension> >
		getGridSizeGlobal() const { return( gridSizeGlobal_ );  }

	Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> >
		getGridSizeLocal() const { return( gridSizeLocal_ );  }

	Teuchos::RCP<const ProcGrid<Ordinal,dimension> >
		getProcGrid() const { return( procGrid_ ); }

	Teuchos::RCP<const GridCoordinatesGlobal<Scalar,Ordinal,dimension> >
		getCoordinatesGlobal() const { return( coordGlobal_ ); }

	Teuchos::RCP<const GridCoordinatesLocal<Scalar,Ordinal,dimension> >
		getCoordinatesLocal() const { return( coordLocal_ ); }

	Teuchos::RCP<const DomainSize<Scalar> >     
		getDomainSize() const { return( domainSize_ ); }

	Teuchos::RCP<const BoundaryConditionsGlobal<dimension> >
		getBCGlobal()   const { return( boundaryConditionsGlobal_ ); }

	Teuchos::RCP<const BoundaryConditionsLocal>
		getBCLocal()    const { return( boundaryConditionsLocal_ ); }

	Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension,dimNC> >
		getInterpolateV2S() const { return( interV2S_ ); }

	/// \}

	/// \name getter methods IMPACT style
	/// \{

	const MPI_Comm& comm()  const { return( procGrid_->getCommS()  ); }

	const MPI_Comm& commST()  const { return( procGrid_->getCommWorld()  ); }

	int rankST() const { return( procGrid_->getRank() ); }
	int rankS () const { return( procGrid_->getRankS() ); }

	const int&      dim()   const { return( getDomainSize()->getDim() ); }

	const Ordinal* nGlo()        const { return( gridSizeGlobal_->getRawPtr()  ); }
	const Ordinal& nGlo( int i ) const { return( gridSizeGlobal_->get(i) ); }

	const Ordinal* nLoc()        const { return( gridSizeLocal_->getRawPtr()  ); }
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

	const Ordinal* ib() const { return( procGrid_->getIB().getRawPtr() ); }

	const Ordinal* getShift()      const { return( indexSpace_->getShift()  ); }
	const Ordinal& getShift(int i) const { return( indexSpace_->getShift(i)  ); }

	const Ordinal* np()      const { return( procGrid_->getNP().getRawPtr() ); }
	const Ordinal& np(int i) const { return( procGrid_->getNP(i) ); }

	/// \}

	void print(  std::ostream& out=std::cout ) const {

		out << "\t---Space: ---\n";

		stencilWidths_->print( out );

		out <<"\t---GridSizeGlobal: ---\n";
		gridSizeGlobal_->print( out );

		out <<"\t---GridSizeLocal: ---\n";
		gridSizeLocal_->print( out );

		indexSpace_->print(out);

		procGrid_->print( out );

		getBCLocal()->print( out );

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
			pl->set<bool>("spectral in time", false, "enables spectral time discretization" );

			pl->set<O>("nx", 33, "amount of grid points in x-direction: a*2**q+1" );
			pl->set<O>("ny", 33, "amount of grid points in y-direction: a*2**q+1" );
			pl->set<O>("nz", 33, "amount of grid points in z-direction: a*2**q+1" );
			pl->set<O>("nf", 4, "amount of grid points in f-direction" );

			// grid stretching
			pl->set<int>("grid stretching in x", 0, "");
			pl->set<int>("grid stretching in y", 0, "");
			pl->set<int>("grid stretching in z", 0, "");

			// processor grid size
			pl->set<O>("npx", 2, "amount of processors in x-direction" );
			pl->set<O>("npy", 2, "amount of processors in y-direction" );
			pl->set<O>("npz", 1, "amount of processors in z-direction" );
			pl->set<O>("npf", 1, "amount of processors in f-direction" );

			validPL = pl;
		}
		return( validPL );
	}


}; // end of class Space



/// \relates Space
template<class S=double, class O=int, int d=3, int dimNC=4>
Teuchos::RCP<const Space<S,O,d,dimNC> >
createSpace( Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ) {

	return( Teuchos::rcp( new Space<S,O,d,dimNC>( pl ) ) );

}


/// \relates Space
/// \relates TransferOp
/// "same" space with changing dimNC
template< class OSpaceT, class ISpaceT> 
Teuchos::RCP< const OSpaceT > createSpace(
		const Teuchos::RCP<const ISpaceT>& space ) {

	auto stencilWidths =
		createStencilWidths< OSpaceT::dimension, OSpaceT::dimNC >(
				space->getStencilWidths()->spectralT() );

	auto domainSize = space->getDomainSize();

	auto boundaryConditionsGlobal = space->getBCGlobal();

	auto boundaryConditionsLocal = space->getBCLocal();

	auto gridSizeGlobal = space->getGridSizeGlobal();

	auto gridSizeLocal = space->getGridSizeLocal();

	auto procGrid = space->getProcGrid();

	auto indexSpace =
		Pimpact::createIndexSpace(
				stencilWidths,
				gridSizeLocal,
				boundaryConditionsLocal,
				procGrid );

	auto  coordGlobal = space->getCoordinatesGlobal();

	auto  coordLocal =
		Pimpact::createGridCoordinatesLocal(
				stencilWidths,
				domainSize,
				gridSizeGlobal,
				gridSizeLocal,
				boundaryConditionsGlobal,
				boundaryConditionsLocal,
				procGrid,
				coordGlobal );

	auto interV2S =
		Pimpact::createInterpolateV2S(
				indexSpace,
				gridSizeLocal,
				stencilWidths,
				domainSize,
				boundaryConditionsLocal,
				coordLocal );

	return(
			Teuchos::rcp(
				new OSpaceT(
					stencilWidths,
					indexSpace,
					gridSizeGlobal,
					gridSizeLocal,
					procGrid,
					coordGlobal,
					coordLocal,
					domainSize,
					boundaryConditionsGlobal,
					boundaryConditionsLocal,
					interV2S ) ) );

}



template<class SpaceT>
static Teuchos::RCP< const SpaceT > createSpace(
		const Teuchos::RCP<const SpaceT>& space,
		const GridSizeGlobal<typename SpaceT::Ordinal,SpaceT::dimension>& newGridSizeGlobal ) {

	using Scalar = typename SpaceT::Scalar;
	using Ordinal = typename SpaceT::Ordinal;
	const int dimension = SpaceT::dimension;
	const int dimNC = SpaceT::dimNC;

	auto stencilWidths = space->getStencilWidths();

	auto domainSize = space->getDomainSize();

	auto boundaryConditionsGlobal = space->getBCGlobal();
	auto boundaryConditionsLocal = space->getBCLocal();


	auto gridSizeGlobal = createGridSizeGlobal<Ordinal,dimension>(
			newGridSizeGlobal );

	auto procGrid = space->getProcGrid();

	auto gridSizeLocal =
		Pimpact::createGridSizeLocal<Ordinal,dimension,dimNC>(
				gridSizeGlobal,
				procGrid,
				stencilWidths );

	auto indexSpace =
		Pimpact::createIndexSpace<Ordinal,dimension>(
				stencilWidths,
				gridSizeLocal,
				boundaryConditionsLocal,
				procGrid	);

	auto coordGlobal =
		Pimpact::createGridCoordinatesGlobal<Scalar,Ordinal,dimension>(
				gridSizeGlobal,
				domainSize,
				Teuchos::tuple( None, None, None) );

	auto coordLocal =
		Pimpact::createGridCoordinatesLocal<Scalar,Ordinal,dimension>(
				stencilWidths,
				domainSize,
				gridSizeGlobal,
				gridSizeLocal,
				boundaryConditionsGlobal,
				boundaryConditionsLocal,
				procGrid,
				coordGlobal );

	auto interV2S =
		Pimpact::createInterpolateV2S<Scalar,Ordinal,dimension>(
				indexSpace,
				gridSizeLocal,
				stencilWidths,
				domainSize,
				boundaryConditionsLocal,
				coordLocal );

	return(
			Teuchos::rcp(
				new SpaceT(
					stencilWidths,
					indexSpace,
					gridSizeGlobal,
					gridSizeLocal,
					procGrid,
					coordGlobal,
					coordLocal,
					domainSize,
					boundaryConditionsGlobal,
					boundaryConditionsLocal,
					interV2S ) ) );

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
