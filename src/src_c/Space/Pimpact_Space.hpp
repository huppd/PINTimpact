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
#include "Pimpact_CoordinatesGlobal.hpp"
#include "Pimpact_CoordinatesLocal.hpp"
#include "Pimpact_DomainSize.hpp"
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
template< class S, class O, int d, int dNC >
class Space {

public:

	using Scalar  = S;
	using Ordinal = O;

	static const int dimension = d;
	static const int dimNC = dNC;


	Space( Teuchos::RCP<Teuchos::ParameterList> pl ) {

		pl->validateParametersAndSetDefaults( *getValidParameters(), 0 );

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
			Pimpact::createGridSizeGlobal<O>(
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
			Pimpact::createIndexSpace<O,d,dNC>(
					stencilWidths_,
					gridSizeLocal_,
					boundaryConditionsLocal_,
					procGrid_ );

		coordGlobal_ =
			Pimpact::createCoordinatesGlobal<S,O,d>(
					gridSizeGlobal_,
					domainSize_,
					Teuchos::tuple< Teuchos::RCP<Teuchos::ParameterList> >(
						Teuchos::rcpFromRef( pl->sublist("Stretching in X") ),
						Teuchos::rcpFromRef( pl->sublist("Stretching in Y") ),
						Teuchos::rcpFromRef( pl->sublist("Stretching in Z") ) ) );

		coordLocal_ =
			Pimpact::createCoordinatesLocal<S,O,d>(
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
			const Teuchos::RCP<const GridSizeGlobal<Ordinal> >& gridSizeGlobal,
			const Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> >& gridSizeLocal,
			const Teuchos::RCP<const ProcGrid<Ordinal,dimension> >& procGrid,
			const Teuchos::RCP<const CoordinatesGlobal<Scalar,Ordinal,dimension> >& coordGlobal,
			const Teuchos::RCP<const CoordinatesLocal<Scalar,Ordinal,dimension> >& coordLocal,
			const Teuchos::RCP<const DomainSize<Scalar> > domainSize,
			const Teuchos::RCP<const BoundaryConditionsGlobal<dimension> > boundaryConditionsGlobal,
			const Teuchos::RCP<const BoundaryConditionsLocal<dimension> > boundaryConditionsLocal,
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

	Teuchos::RCP<const GridSizeGlobal<Ordinal> >
	gridSizeGlobal_;

	Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> >
	gridSizeLocal_;

	Teuchos::RCP<const ProcGrid<Ordinal,dimension> >
	procGrid_;

	Teuchos::RCP<const CoordinatesGlobal<Scalar,Ordinal,dimension> >
	coordGlobal_;

	Teuchos::RCP<const CoordinatesLocal<Scalar,Ordinal,dimension> >
	coordLocal_;

	Teuchos::RCP<const DomainSize<Scalar> >
	domainSize_;

	Teuchos::RCP<const BoundaryConditionsGlobal<dimension> >
	boundaryConditionsGlobal_;

	Teuchos::RCP<const BoundaryConditionsLocal<dimension> >
	boundaryConditionsLocal_;

	Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension,dimNC> >
	interV2S_;

public:

	/// \name getter methods
	/// \{

	const Teuchos::RCP<const StencilWidths<dimension,dimNC> >&
	getStencilWidths() const { return( stencilWidths_ ); }

	const Teuchos::RCP<const IndexSpace<Ordinal,dimension> >&
	getIndexSpace() const { return( indexSpace_ ); }

	const Teuchos::RCP<const GridSizeGlobal<Ordinal> >&
	getGridSizeGlobal() const { return( gridSizeGlobal_ );  }

	const Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> >&
	getGridSizeLocal() const { return( gridSizeLocal_ );  }

	const Teuchos::RCP<const ProcGrid<Ordinal,dimension> >&
	getProcGrid() const { return( procGrid_ ); }

	const Teuchos::RCP<const CoordinatesGlobal<Scalar,Ordinal,dimension> >&
	getCoordinatesGlobal() const { return( coordGlobal_ ); }

	const Teuchos::RCP<const CoordinatesLocal<Scalar,Ordinal,dimension> >&
	getCoordinatesLocal() const { return( coordLocal_ ); }

	const Teuchos::RCP<const DomainSize<Scalar> >&
	getDomainSize() const { return( domainSize_ ); }

	const Teuchos::RCP<const BoundaryConditionsGlobal<dimension> >&
	getBCGlobal()   const { return( boundaryConditionsGlobal_ ); }

	const Teuchos::RCP<const BoundaryConditionsLocal<dimension> >&
	getBCLocal()    const { return( boundaryConditionsLocal_ ); }

	const Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension,dimNC> >&
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
			pl->sublist("Stretching in X");
			pl->sublist("Stretching in Y");
			pl->sublist("Stretching in Z");
			pl->sublist("Stretching in X").set<std::string>( "Stretch Type", "none" );
			pl->sublist("Stretching in Y").set<std::string>( "Stretch Type", "none" );
			pl->sublist("Stretching in Z").set<std::string>( "Stretch Type", "none" );

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
/// \relates CoarsenStrategy
/// \relates CoarsenStrategyGlobal
/// "same" space with changing dimNC
template< class OSpaceT, class ISpaceT> 
Teuchos::RCP< const OSpaceT >
createSpace(
		const Teuchos::RCP<const ISpaceT>& space ) {

	using Scalar    = typename OSpaceT::Scalar;
	using Ordinal   = typename OSpaceT::Ordinal;
	const int dim   = OSpaceT::dimension;
	const int dimNC = OSpaceT::dimNC;

	Teuchos::RCP<const StencilWidths< dim, dimNC > > stencilWidths =
		createStencilWidths< dim, dimNC >(
				space->getStencilWidths()->spectralT() );

	Teuchos::RCP<const DomainSize<Scalar> > domainSize = space->getDomainSize();

	Teuchos::RCP<const BoundaryConditionsGlobal<dim> > boundaryConditionsGlobal = space->getBCGlobal();

	Teuchos::RCP<const BoundaryConditionsLocal<dim> > boundaryConditionsLocal = space->getBCLocal();

	Teuchos::RCP<const GridSizeGlobal<Ordinal> > gridSizeGlobal = space->getGridSizeGlobal();

	Teuchos::RCP<const GridSizeLocal<Ordinal,dim> > gridSizeLocal = space->getGridSizeLocal();

	Teuchos::RCP<const ProcGrid<Ordinal, dim> > procGrid = space->getProcGrid();

	Teuchos::RCP<const IndexSpace<Ordinal,dim> > indexSpace =
		Pimpact::createIndexSpace<Ordinal,dim,dimNC>(
				stencilWidths,
				gridSizeLocal,
				boundaryConditionsLocal,
				procGrid );

	Teuchos::RCP<const CoordinatesGlobal<Scalar,Ordinal,dim> > coordGlobal = space->getCoordinatesGlobal();

	Teuchos::RCP<const CoordinatesLocal<Scalar,Ordinal,dim> >  coordLocal =
		Pimpact::createCoordinatesLocal(
				stencilWidths,
				domainSize,
				gridSizeGlobal,
				gridSizeLocal,
				boundaryConditionsGlobal,
				boundaryConditionsLocal,
				procGrid,
				coordGlobal );

	Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dim,dimNC> > interV2S =
		Pimpact::createInterpolateV2S<Scalar,Ordinal,dim,dimNC>(
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

} // end of createSpace



/// \brief creates coarse space from fine and new GridSize
///
/// \tparam SpaceT type of space
/// \param space space object
/// \param newGridSizeGlobal new grid size
///
/// \return new coarse space
template<class SpaceT>
static Teuchos::RCP< const SpaceT > createSpace(
		const Teuchos::RCP<const SpaceT>& space,
		const GridSizeGlobal<typename SpaceT::Ordinal>& newGridSizeGlobal ) {

	using Scalar = typename SpaceT::Scalar;
	using Ordinal = typename SpaceT::Ordinal;
	const int dim = SpaceT::dimension;
	const int dimNC = SpaceT::dimNC;

	Teuchos::RCP<const StencilWidths<dim,dimNC> > stencilWidths = space->getStencilWidths();

	Teuchos::RCP<const DomainSize<Scalar> > domainSize = space->getDomainSize();

	Teuchos::RCP<const BoundaryConditionsGlobal<dim> > boundaryConditionsGlobal = space->getBCGlobal();
	Teuchos::RCP<const BoundaryConditionsLocal <dim> > boundaryConditionsLocal = space->getBCLocal();


	Teuchos::RCP<const GridSizeGlobal<Ordinal> > gridSizeGlobal =
		createGridSizeGlobal<Ordinal>(
				newGridSizeGlobal );

	Teuchos::RCP<const ProcGrid<Ordinal,dim> > procGrid = space->getProcGrid();

	Teuchos::RCP<const GridSizeLocal<Ordinal,dim> > gridSizeLocal =
		Pimpact::createGridSizeLocal<Ordinal,dim,dimNC>(
				gridSizeGlobal,
				procGrid,
				stencilWidths );

	Teuchos::RCP<const IndexSpace<Ordinal,dim> > indexSpace =
		Pimpact::createIndexSpace<Ordinal,dim,dimNC>(
				stencilWidths,
				gridSizeLocal,
				boundaryConditionsLocal,
				procGrid	);

	Teuchos::RCP<const CoordinatesGlobal<Scalar,Ordinal,dim> > coordGlobal =
		Pimpact::createCoordinatesGlobal<Scalar,Ordinal,dim>(
				gridSizeGlobal,
				space->getCoordinatesGlobal() );

	Teuchos::RCP<const CoordinatesLocal<Scalar,Ordinal,dim> > coordLocal =
		Pimpact::createCoordinatesLocal<Scalar,Ordinal,dim>(
				stencilWidths,
				domainSize,
				gridSizeGlobal,
				gridSizeLocal,
				boundaryConditionsGlobal,
				boundaryConditionsLocal,
				procGrid,
				coordGlobal );

	Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dim,dimNC> > interV2S =
		Pimpact::createInterpolateV2S<Scalar,Ordinal,dim,dimNC>(
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

} // end of createSpace



/// \todo move ProcGrid createion to createProcGrid(...)
/// \todo redo input GridSizeGlobal, ProcGridSize
/// \relates CoarsenStrategyGlobal
template<class SpaceT>
static Teuchos::RCP< const SpaceT > createSpace(
		const Teuchos::RCP<const SpaceT>& space,
		const GridSizeGlobal<typename SpaceT::Ordinal>& newGridSizeGlobal,
		const Teuchos::Tuple<typename SpaceT::Ordinal,SpaceT::dimension>& npNew,
		Teuchos::Tuple<typename SpaceT::Ordinal,SpaceT::dimension>& stride,
		const Teuchos::Tuple<typename SpaceT::Ordinal,SpaceT::dimension>& npWorld,
		const Teuchos::Tuple<typename SpaceT::Ordinal,SpaceT::dimension>& ibWorld ) {

	typedef typename SpaceT::Scalar  Scalar;
	typedef typename SpaceT::Ordinal Ordinal;

	const int dimension = SpaceT::dimension;
	const int dimNC = SpaceT::dimNC;

	typedef typename Teuchos::Tuple<Ordinal,dimension> TO;

	Teuchos::RCP<const StencilWidths<dimension,dimNC> >
		stencilWidths = space->getStencilWidths();

	Teuchos::RCP<const DomainSize<Scalar> > domainSize = space->getDomainSize();
	Teuchos::RCP<const BoundaryConditionsGlobal<dimension> > boundaryConditionsGlobal = space->getBCGlobal();


	// --- coarsen gridSizeGlobal ---
	auto gridSizeGlobal =
		createGridSizeGlobal<Ordinal>( newGridSizeGlobal );


	/// redo procGrid create new communicator, make procgrid from communicator
	TO np = space->getProcGrid()->getNP();
	TO nGather;
	TO ib;
	for( Ordinal dir=0; dir<dimension; ++dir ) {
		nGather[dir] = np[dir]/npNew[dir];
		stride[dir] *= nGather[dir];
		ib[dir]      = ( ibWorld[dir] - 1 )*npNew[dir]/npWorld[dir] + 1;
	}

	bool participating = false;
	MPI_Comm commWorld = space->getProcGrid()->getCommWorld();
	MPI_Comm commSub ;//= space->getProcGrid()->getCommS();

	int rankWorld = space->getProcGrid()->getRank();
	int rankSub = space->getProcGrid()->getRank(); // necessary?

	Teuchos::Tuple<int,dimension> rankL;
	Teuchos::Tuple<int,dimension> rankU;
	for( int dir=0; dir<dimension; ++dir )
		MPI_Cart_shift(
				commWorld,      // communicator with Cartesian structure
				dir,            // coordinate dimension of shift
				stride[dir],    // displacement
				&rankL[dir],    // rank of source process
				&rankU[dir] );	// rank of destination process

	Ordinal gather_yes = 1;
	for( int i=0; i<dimension; ++i ) {
		gather_yes *= nGather[i];
	}
	if( gather_yes>1 ) {
		int n = 1;
		for( int i=0; i<dimension; ++i ) 
			n *= npNew[i];
		int* newRanks = new int[n];

		TO rankCoord;

		for( int i=0; i<npNew[0]; ++i) {
			rankCoord[0] = (i*stride[0])%npWorld[0];
			for( int j=0; j<npNew[1]; ++j ) {
				rankCoord[1] = (j*stride[1])%npWorld[1];
				for( int k=0; k<npNew[2]; ++k ) {
					rankCoord[2] = (k*stride[2])%npWorld[2];
					if( 4==dimension ) {
						for( int l=0; l<npNew[3]; ++l ) {
							rankCoord[3] = (l*stride[3])%npWorld[3];
							MPI_Cart_rank(
									commWorld,									          // comm
									rankCoord.getRawPtr(),                // processor coordinate
									&newRanks[i
									+j*npNew[0]
									+k*npNew[0]*npNew[1]
									+l*npNew[0]*npNew[1]*npNew[2] ] );   	// according rank to coordinate
							if( rankWorld==newRanks[i + j*npNew[0] + k*npNew[0]*npNew[1] + l*npNew[0]*npNew[1]*npNew[2] ] )
								participating = true;
						}
					}
					else {
						MPI_Cart_rank(
								commWorld,									                    // comm
								rankCoord.getRawPtr(),                          // processor coordinate
								&newRanks[i+j*npNew[0]+k*npNew[0]*npNew[1] ] ); // according rank to coordinate
						if( rankWorld==newRanks[i+j*npNew[0]+k*npNew[0]*npNew[1] ])
							participating = true;
					}
				}
			}
		}

		MPI_Comm commTemp;
		MPI_Group baseGroup, newGroup;

		MPI_Comm_group( commWorld, &baseGroup );
		MPI_Group_incl( baseGroup, n, newRanks, &newGroup );
		MPI_Comm_create( commWorld, newGroup, &commTemp );
		MPI_Group_free( &baseGroup );
		MPI_Group_free( &newGroup );

		Teuchos::Tuple<int,dimension> periodic = boundaryConditionsGlobal->periodic();

		if( participating ) {

			MPI_Cart_create(
					commTemp,		          // communicator without Cartesian information
					dimension,            // number of dimensions
					npNew.getRawPtr(),    // number of processors in each dimension
					periodic.getRawPtr(),	// array for mpi to signal which dimension is periodic
					false,                // false means ranking is not reordered
					&commSub );           // new communicator with Cartesian information
			if( 4==dimension ) {
				MPI_Comm commTemp_;
				int temp[] = {1,1,1,0};
				MPI_Cart_sub( commSub, temp, &commTemp_ );
				MPI_Comm_free( &commSub );
				commSub = commTemp_;
			}
			MPI_Comm_free( &commTemp );
		}
		else
			commSub=MPI_COMM_NULL;

		delete[] newRanks;

	}

	if( commSub==MPI_COMM_NULL )
		rankSub = -1;
	else 
		MPI_Comm_rank( commSub, &rankSub ); // get rank

	Teuchos::RCP<const ProcGrid<Ordinal,dimension> > procGrid =
		Teuchos::rcp(
				new ProcGrid<Ordinal,dimension>(
					npNew,
					participating,
					commWorld,
					commSub,
					rankWorld,
					rankSub,
					ib,
					rankL,
					rankU ) );


	Teuchos::RCP<const GridSizeLocal<Ordinal,dimension> > gridSizeLocal =
		Pimpact::createGridSizeLocal<Ordinal,dimension,dimNC>(
				gridSizeGlobal,
				procGrid,
				stencilWidths );

	Teuchos::RCP<const BoundaryConditionsLocal<dimension> > boundaryConditionsLocal =
		createBoudaryConditionsLocal<Ordinal,dimension>( 
				boundaryConditionsGlobal,
				procGrid );

	Teuchos::RCP<const IndexSpace<Ordinal,dimension> > indexSpace =
		Pimpact::createIndexSpace<Ordinal,dimension,dimNC>(
				stencilWidths,
				gridSizeLocal,
				boundaryConditionsLocal,
				procGrid );

	Teuchos::RCP<const CoordinatesGlobal<Scalar,Ordinal,dimension> > coordGlobal =
		Pimpact::createCoordinatesGlobal<Scalar,Ordinal,dimension>(
				gridSizeGlobal,
				space->getCoordinatesGlobal() );
	//Teuchos::RCP<const CoordinatesGlobal<Scalar,Ordinal,dimension> > coordGlobal =
		//Pimpact::createCoordinatesGlobal<Scalar,Ordinal,dimension>(
				//gridSizeGlobal,
				//domainSize,
				//Teuchos::tuple( None, None, None) );

	Teuchos::RCP<const CoordinatesLocal<Scalar,Ordinal,dimension> > coordLocal =
		Pimpact::createCoordinatesLocal<Scalar,Ordinal,dimension>(
				stencilWidths,
				domainSize,
				gridSizeGlobal,
				gridSizeLocal,
				boundaryConditionsGlobal,
				boundaryConditionsLocal,
				procGrid,
				coordGlobal );

	Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,dimension,dimNC> > interV2S =
		Pimpact::createInterpolateV2S<Scalar,Ordinal,dimension,dimNC>(
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
					interV2S )
				)
			);

}






#ifdef COMPILE_ETI
extern template class Space<double,int,3,2>;
extern template class Space<double,int,3,4>;
extern template class Space<double,int,4,2>;
extern template class Space<double,int,4,4>;
#endif



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SPACE_HPP
