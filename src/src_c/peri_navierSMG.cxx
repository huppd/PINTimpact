#include <mpi.h>

#include <ostream>
#include <fstream>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include "Teuchos_Range1D.hpp"
#include <Teuchos_CommHelpers.hpp>
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

#include "BelosOutputManager.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosTypes.hpp"

#include "pimpact.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_DomainSize.hpp"
#include "Pimpact_Space.hpp"
#include "Pimpact_ProcGridSize.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_ModeNonlinearOp.hpp"
#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "Pimpact_LinSolverParameter.hpp"

#include "Pimpact_CoarsenStrategyGlobal.hpp"
#include "Pimpact_RefinementStrategy.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_TransferMultiHarmonicOp.hpp"
#include "Pimpact_TransferCompoundOp.hpp"

#include "BelosPimpactAdapter.hpp"

#include "NOX_Pimpact_Vector.hpp"
#include "NOX_Pimpact_Interface.hpp"
#include "NOX_Pimpact_Group.hpp"

#include "NOX_Pimpact_StatusTest.hpp"

#include "NOX.H"



/// \brief gets \c ParameterList from command Line
template<class S, class O>
Teuchos::RCP<Teuchos::ParameterList>
getSpaceParametersFromCL( int argi, char** argv  )  {

	Teuchos::CommandLineProcessor my_CLP;

  std::string xmlFilename = "bla.xml";

	// Space parameters
  // physical constants
	S re = 100.;

	S alpha2 = 24.;

  // domain type
  int domain = 0;

  // domain size
  int dim = 3;
  S lx = 8.;
  S ly = 2.;
  S lz = 4.;

  // grid size
	O nx = 129;
	O ny = 33;
	O nz = 65;
	O nf = 2;

	O nfs = 1;
	O nfe = 1;

  // processor grid size
  O npx = 1;
  O npy = 1;
  O npz = 1.;
  O npf = 1.;

	// flow and forcing type
  int baseflow = 1;
  int flow = 0;
  int forcing = 1;

	// solver parameters
  std::string nonLinSolName = "Newton";

  std::string lineSearchName = "Backtrack";

  std::string linSolName = "";

  int withprec=2;
  int withoutput=1;

  int refinement=1;

  int maxIter = 10;

  int initZero = 0;

  S tolBelos = 1.e-1;
  S tolInnerBelos = 1.e-2;
  S tolNOX   = 1.e-6;

  // multi grid parameters
  int maxGrids = 10;
  int numCycles = 2;

  my_CLP.setOption("filename", &xmlFilename, "file name of the input xml paramerterlist");
  // Space parameters
  my_CLP.setOption( "re", &re, "Reynolds number" );

  my_CLP.setOption( "alpha2", &alpha2, "introduced frequency" );

  my_CLP.setOption( "domain", &domain,
      "Domain type: 0:all dirichlet, 1:dirichlet 2d channel, 2: periodic 2d channel" );

  // domain size
  my_CLP.setOption( "lx", &lx, "length in x-direction" );
  my_CLP.setOption( "ly", &ly, "length in y-direction" );
  my_CLP.setOption( "lz", &lz, "length in z-direction" );

  my_CLP.setOption( "dim", &dim, "dimension of problem" );

  // grid size
  my_CLP.setOption( "nx", &nx, "amount of grid points in x-direction: a*2**q+1" );
  my_CLP.setOption( "ny", &ny, "amount of grid points in y-direction: a*2**q+1" );
  my_CLP.setOption( "nz", &nz, "amount of grid points in z-direction: a*2**q+1" );
  my_CLP.setOption( "nf", &nf, "amount of grid points in f-direction" );

  my_CLP.setOption( "nfs", &nfs, "start amount of grid points in f-direction" );
  my_CLP.setOption( "nfe", &nfe, "end amount of grid points in f-direction" );

  // processor grid size
  my_CLP.setOption( "npx", &npx, "amount of processors in x-direction" );
  my_CLP.setOption( "npy", &npy, "amount of processors in y-direction" );
  my_CLP.setOption( "npz", &npz, "amount of processors in z-direction" );
  my_CLP.setOption( "npf", &npf, "amount of processors in f-direction" );

	// flow and forcint
  my_CLP.setOption( "baseflow", &baseflow,
      "Flow type" );
  my_CLP.setOption( "flow", &flow,
      "Flow type" );
  my_CLP.setOption( "force", &forcing,
      "forcing, ja?" );

	// solver parameters
  my_CLP.setOption("nonLinSolver", &nonLinSolName, "Newton");
  my_CLP.setOption("linSearch", &lineSearchName, "Backtrack");
  my_CLP.setOption("linSolver", &linSolName, "bla");
  my_CLP.setOption("withprec", &withprec, "0: no preconditioner, 1: block triangular, 2: block triangular with inner MG multigrid, 3: block triangular with mg solver");
  my_CLP.setOption("withoutput", &withoutput, "0: no preconditioner, 1: block triangular, 2: block triangular with inner MG multigrid, 3: block triangular with mg solver");
  my_CLP.setOption("refinement", &refinement, "");
  my_CLP.setOption("maxIter", &maxIter, "");
  my_CLP.setOption("initZero", &initZero, "");
  my_CLP.setOption("tolBelos", &tolBelos, "");
  my_CLP.setOption("tolInnerBelos", &tolInnerBelos, "");
  my_CLP.setOption("tolNOX", &tolNOX, "");

	// multi grid parameters
  my_CLP.setOption("maxGrids", &maxGrids, "");
  my_CLP.setOption("numCycles", &numCycles, "");


  my_CLP.recogniseAllOptions(true);
  my_CLP.throwExceptions(true);

  my_CLP.parse(argi,argv);


	// Set all the valid parameters and their default values.
	Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList("Pimpact");

  pl->set("XML", xmlFilename );

	// Space parameters
	pl->sublist("Space").set<S>("Re", re, "Reynolds number");
	pl->sublist("Space").set<S>("alpha2", alpha2,
			"Womersley square \alpha^2");
	// domain type
	pl->sublist("Space").set<int>( "domain", domain,
			"Domain type: 0:all dirichlet, 1:dirichlet 2d channel, 2: periodic 2d channel" );

	// domain size
	pl->sublist("Space").set<int>("dim", dim, "dimension of problem" );

	pl->sublist("Space").set<S>( "lx", lx, "length in x-direction" );
	pl->sublist("Space").set<S>( "ly", ly, "length in y-direction" );
	pl->sublist("Space").set<S>( "lz", lz, "length in z-direction" );


	// grid size
	pl->sublist("Space").set<O>("nx", nx, "amount of grid points in x-direction: a*2**q+1" );
	pl->sublist("Space").set<O>("ny", ny, "amount of grid points in y-direction: a*2**q+1" );
	pl->sublist("Space").set<O>("nz", nz, "amount of grid points in z-direction: a*2**q+1" );
	pl->sublist("Space").set<O>("nf", nf, "amount of grid points in f-direction" );

	pl->sublist("Space").set<O>("nfs", nfs, "start amount of grid points in f-direction" );
	pl->sublist("Space").set<O>("nfe", nfe, "end amount of grid points in f-direction" );

	// processor grid size
	pl->sublist("Space").set<O>("npx", npx, "amount of processors in x-direction" );
	pl->sublist("Space").set<O>("npy", npy, "amount of processors in y-direction" );
	pl->sublist("Space").set<O>("npz", npz, "amount of processors in z-direction" );
	pl->sublist("Space").set<O>("npf", npf, "amount of processors in f-direction" );

	// flow and forcint
	pl->set<int>( "baseflow", baseflow, "Flow type: depending on main" );
	pl->set<int>( "flow", flow, "Flow type: depending on main" );
	pl->set<int>( "forcing", forcing, "forcing, ja?" );

	// solver parameters
  pl->sublist("Solver").set("nonLinSolver", nonLinSolName );
  pl->sublist("Solver").set("linSearch", lineSearchName );
  pl->sublist("Solver").set("linSolver", linSolName );
  pl->sublist("Solver").set("withprec", withprec );
  pl->sublist("Solver").set("withoutput", withoutput );
  pl->sublist("Solver").set("refinement", refinement );
  pl->sublist("Solver").set("maxIter", maxIter );
  pl->sublist("Solver").set("initZero", initZero );
  pl->sublist("Solver").set("tolBelos", tolBelos );
  pl->sublist("Solver").set("tolInnerBelos", tolInnerBelos );
  pl->sublist("Solver").set("tolNOX", tolNOX );

	// Multi Grid parameters
  pl->sublist("Multi Grid").set("maxGrids", maxGrids );
  pl->sublist("Multi Grid").set("numCycles", numCycles );

	return( pl );
}



typedef double S;
typedef int O;

typedef Pimpact::Space<S,O,3,4> SpaceT;

typedef Pimpact::Space<S,O,3,4> FSpaceT;
typedef Pimpact::Space<S,O,3,2> CSpaceT;

typedef Pimpact::CoarsenStrategyGlobal<FSpaceT,CSpaceT> CS;

typedef Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> > VF;
typedef Pimpact::MultiHarmonicField< Pimpact::ScalarField<SpaceT> > SF;

typedef Pimpact::MultiField<VF> MVF;
typedef Pimpact::MultiField<SF> MSF;

typedef Pimpact::CompoundField< VF, SF> CF;
typedef Pimpact::MultiField<CF> MF;

typedef Pimpact::OperatorBase<MF> BOp;

typedef NOX::Pimpact::Interface<MF> Inter;
typedef NOX::Pimpact::Vector<typename Inter::Field> NV;


template<class T> using FT =
	Pimpact::CompoundField<
		Pimpact::MultiHarmonicField< Pimpact::VectorField<T> >,
		Pimpact::MultiHarmonicField< Pimpact::ScalarField<T> > >;

template<class T1,class T2>
using TransF = Pimpact::TransferCompoundOp<
	Pimpact::TransferMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::TransferOp<T1,T2>  >  >,
	Pimpact::TransferMultiHarmonicOp<                             Pimpact::TransferOp<T1,T2>  >  >;


template<class T>     
using RestrF = Pimpact::TransferCompoundOp<
	Pimpact::TransferMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::RestrictionHWOp<T> > >,
	Pimpact::TransferMultiHarmonicOp<                             Pimpact::RestrictionHWOp<T> > >;


template<class T> 
using InterF = Pimpact::TransferCompoundOp<
	Pimpact::TransferMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp<T> > >,
	Pimpact::TransferMultiHarmonicOp<                             Pimpact::InterpolationOp<T> > >;


template<class T> using DTConvDiffOpT =
	Pimpact::CompoundOpWrap<
		Pimpact::MultiDtConvectionDiffusionOp<T>,
		Pimpact::MultiHarmonicOpWrap< Pimpact::GradOp<T> >,
		Pimpact::MultiHarmonicOpWrap< Pimpact::DivOp<T > > >;

template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;


//template<class T> using SmootherT =
//	Pimpact::InverseTriangularOp<
//
//		,Pimpact::MultiHarmonicOpWrap<Pimpact::GradOp<T> >,
//		,Pimpact::MultiHarmonicOpWrap<Pimpact::DivGradO2JSmoother<Pimpact::DivGradO2Op<T>	> > >;






int main(int argi, char** argv ) {

	/////////////////////////////////////////// set up parameters ///////////////////////////
  // intialize MPI
  MPI_Init( &argi, &argv );

	Teuchos::CommandLineProcessor my_CLP;
  std::string xmlFilename = "parameterIn.xml";

//	Teuchos::RCP<Teuchos::ParameterList> plup = getSpaceParametersFromCL<S,O>( argi, argv  ) ;

	auto plXML = plup->get<std::string>("XML");

	auto pl = Teuchos::getParametersFromXmlFile( xmlFielname );
	pl->print();


	//	solver stuff
	std::string nonLinSolName = pl->sublist("Solver").get<std::string>("nonLinSolver");
	std::string lineSearchName = pl->sublist("Solver").get<std::string>( "linSearch" );

	S tolNOX   = pl->sublist("Solver").get<S>("tolNOX");
	int maxIter = pl->sublist("Solver").get<int>("maxIter");
	int initZero = pl->sublist("Solver").get<int>("initZero");
	S tolBelos = pl->sublist("Solver").get<S>("tolBelos");
	S tolInnerBelos = pl->sublist("Solver").get<S>("tolInnerBelos");

	int withprec=pl->sublist("Solver").get<int>("withprec");
	int withoutput=pl->sublist("Solver").get<int>("withoutput");

	int refinement=pl->sublist("Solver").get<int>("refinement");

	std::string linSolName  = pl->sublist("Solver").get<std::string>( "linSolver" );
	if( linSolName.size()< 2 )
		linSolName = (withprec>0)?"Block GMRES":"GMRES";

	auto space = Pimpact::createSpace<S,O,3,4>( Teuchos::rcpFromRef( pl->sublist("Space", true) ) );


//	int flow = pl->get<int>("flow");
	int force = pl->get<int>("forcing");


	/////////////////////////////////////////// end of set up parameters ///////////////////////////
	///////////////////////////////////////////  set up initial stuff //////////////////////////////
	
	// init vectors
	auto x = Pimpact::create<MF>( space );
//		Pimpact::createMultiField(
//			Pimpact::createCompoundField(
//				Pimpact::createMultiHarmonic< Pimpact::VectorField<SpaceT> >( space ),
//				Pimpact::createMultiHarmonic< Pimpact::ScalarField<SpaceT> >( space )) );

	// init Fields
	x->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->initField( Pimpact::EVectorField( pl->get<int>("baseflow") ), 1. );

	if( 0==initZero )
		x->init(0.);
	else if( 1==initZero ) {
		x->getFieldPtr(0)->getVFieldPtr()->random();
		x->getFieldPtr(0)->getSFieldPtr()->init(0.);
		x->scale(1.e-32);
	}
	else if( 2==initZero )
		x->random();

	////////////////////////////////////// end of set up initial stuff //////////////////////////////

	//tolBelos*=l1*l2/n1/n2*(nfe-1)/nfs;

  /******************************************************************************************/
	for( int refine=0; refine<refinement; ++refine ) {

		std::string rl = "";
		if( refinement>1 )
			rl = std::to_string( (long long)refine ); // long long needed on brutus(intel)

		auto fu   = x->clone( Pimpact::ShallowCopy );
		if( 0==force )
			fu->init( 0. );
		else {
//			S re = space->getDomainSize()->getRe();
			fu->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->getFieldPtr(Pimpact::U)->initField( Pimpact::FPoint,  -2. );
			fu->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->getFieldPtr(Pimpact::V)->initField( Pimpact::FPoint, pl->get<S>( "lambda", 1. ) );
			//		fu->getFieldPtr(0)->getVFieldPtr()->getSFieldPtr(0)->getFieldPtr(Pimpact::W)->initField( Pimpact::FPoint, 1. );
		}
		//	fu->getFieldPtr(0)->getVFieldPtr()->write( 700,true );


		auto opV2V = Pimpact::createMultiDtConvectionDiffusionOp( space );
		auto opS2V = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::GradOp>( space ) );
		auto opV2S = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::DivOp>( space ) );

		auto op =
			Pimpact::createMultiOperatorBase(
					Pimpact::createCompoundOpWrap(
						opV2V,
						opS2V,
						opV2S )
					);


    Teuchos::RCP<BOp> jop;
		jop = op;

    Teuchos::RCP<BOp> lprec;

		/*** init preconditioner ****************************************************************************/
		if( withprec>0 ) {

			auto mgSpaces = Pimpact::createMGSpaces<FSpaceT,CSpaceT,CS>( space, pl->sublist("Multi Grid").get<int>("maxGrids") );

      auto pls = Teuchos::parameterList();
			pls->set( "numCycles", pl->sublist("Multi Grid").get<int>("numCycles") );
//			pls->sublist("Smoother").set( "omega", 0.5 );
//      pls->sublist("Smoother").set<int>( "Ordering", 1 );
//			pls->sublist("Smoother").set( "numIters", 10 );
//			pls->sublist("Smoother").set( "numIters", 10 );

			pls->sublist("Smoother").set<std::string>("Solver name", "GMRES" );
			pls->sublist("Smoother").set( "Solver",
					*Pimpact::createLinSolverParameter( "GMRES", 1.e-16, -1,
						Teuchos::rcp<std::ostream>( new Teuchos::oblackholestream() ), 4,
						"Smoother" ) );
			pls->sublist("Coarse Grid Solver").set<std::string>( "Solver name", "GMRES" );
			pls->sublist("Coarse Grid Solver").set( "Solver",
					*Pimpact::createLinSolverParameter( "GMRES", 1.e-3, -1,
//						Teuchos::rcp<std::ostream>( &std::cout, false ), 100,
						Teuchos::rcp<std::ostream>( new Teuchos::oblackholestream() ), 100,
						"Coarse Grid Solver" )
					);

//			auto mg = Pimpact::createMGFields<FT>( mgSpaces );
//			auto mgop = Pimpact::createMGOperators<DTConvDiffOpT,DTConvDiffOpT>(mgSpaces);
//			auto mgsop = Pimpact::createMGSmoothers<MOP>( mgop, pl );
//
//			auto mgtra = Pimpact::createMGTransfers<TransF,RestrF,InterF>(mgSpaces);

			auto mg =
				Pimpact::createMultiGrid<
					FT,
					TransF,
					RestrF,
					InterF,
					DTConvDiffOpT,
					DTConvDiffOpT,
					MOP,
					MOP > ( mgSpaces, pls ) ;

			mg->print();

			lprec = Pimpact::createMultiOperatorBase( mg );

		}

		/*** end of init preconditioner ****************************************************************************/


		auto para = Pimpact::createLinSolverParameter( linSolName, tolBelos, 1, Pimpact::createOstream( "stats_linSolve"+rl+".txt",  space->rankST() ), 100, op->getLabel() );

		para->set( "Num Blocks",         5   );

		para->set( "Verbosity",	
				Belos::Errors +
				Belos::Warnings +
				Belos::IterationDetails +
				Belos::OrthoDetails +
				Belos::FinalSummary +
				Belos::StatusTestDetails +
				Belos::Debug
				);

    auto lp_ = Pimpact::createLinearProblem<MF>( jop, x->clone(), fu->clone(), para, linSolName );
		if( withprec ) lp_->setRightPrec( lprec );
//		if( withprec ) lp_->setLeftPrec( lprec );


    auto lp = Pimpact::createInverseOperatorBase( lp_ );


    auto inter = NOX::Pimpact::createInterface( fu, op, lp );

    auto nx = NOX::Pimpact::createVector(x);

    auto bla = Teuchos::parameterList();

    auto group = NOX::Pimpact::createGroup<Inter>( bla, inter, nx );

    // Set up the status tests
    auto statusTest = NOX::Pimpact::createStatusTest( maxIter, tolNOX, 1.e-16 );

    // Create the list of solver parameters
    auto solverParametersPtr =
        NOX::Pimpact::createNOXSolverParameter( nonLinSolName, lineSearchName );

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver =
        NOX::Solver::buildSolver( group, statusTest, solverParametersPtr);

		if( 0==space->rankST() ) std::cout << "\n\t--- Nf: "<< space->nGlo(3) <<"\tdof: "<<x->getLength(true)<<"\t---\n";

    // Solve the nonlinear system
		{
			Teuchos::TimeMonitor LocalTimer(*Teuchos::TimeMonitor::getNewCounter("Pimpact:: Solving Time"));
			try{
				solver->solve();
			}
			catch( std::logic_error & e ) {
				std::cout << e.what() << "\n";
			}
		}

		Teuchos::TimeMonitor::summarize();

    // Get the answer
    *group = solver->getSolutionGroup();

		if( withoutput ) {
			x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr();
			Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr()->write( 800 );
			Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr()->getFieldPtr(0)->getVFieldPtr()->write( 400, true );
			//		Teuchos::rcp_dynamic_cast<const NV>( group->getFPtr() )->getConstFieldPtr()->write(900);
		}

		// spectral refinement of x, fu
		if( refinement>1 ) {
			auto spaceF =
				Pimpact::RefinementStrategy<SpaceT>::createRefinedSpace(
//						space, Teuchos::tuple( false,false,false,true ) );
				space, Teuchos::tuple(1,1,1,1) );

			auto refineOp =
				Teuchos::rcp(
						new Pimpact::TransferCompoundOp<
							Pimpact::TransferMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp<SpaceT> > >,
							Pimpact::TransferMultiHarmonicOp< Pimpact::InterpolationOp<SpaceT> >
							>( space, spaceF ) );
			//		refineOp->print();

			auto xf = Pimpact::create<CF>( spaceF );
			// init Fields
			//		boundaries
			xf->getVFieldPtr()->get0FieldPtr()->initField( Pimpact::EVectorField(baseflow), 1. );
			auto temp = Pimpact::create<CF>( spaceF );

			refineOp->apply( x->getField(0), *temp );

			xf->add( 1., *temp, 0., *temp );

			xf->write();

			x = Pimpact::createMultiField( xf );
			space = spaceF;
		}
  } 
	/******************************************************************************************/


	Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );

  MPI_Finalize();
  return( 0 );

}
