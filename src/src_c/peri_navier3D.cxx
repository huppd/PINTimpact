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
#include "Teuchos_oblackholestream.hpp"

#include "pimpact.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_DomainSize.hpp"
#include "Pimpact_Space.hpp"
#include "Pimpact_ProcGridSize.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"

#include "Pimpact_LinSolverParameter.hpp"

#include "Pimpact_MultiGrid.hpp"

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

	// Space parameters
  // physical constants
	S re = 200.;

	S alpha2 = 64.;

  // domain type
  int domain = 1;

  // domain size
  int dim = 3;
  S lx = 8.;
  S ly = 2.;
  S lz = 4.;

  // grid size
	O nx = 65;
	O ny = 33;
	O nz = 65;
	O nf = 1;

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
  int forcing = 0;

	// solver parameters
  std::string nonLinSolName = "Newton";

  std::string lineSearchName = "Backtrack";

  std::string linSolName = "";

	int withprec=2;

  int maxIter = 10;

  int initZero = 1;

  S tolBelos = 1.e-1;
  S tolInnerBelos = 1.e-3;
  S tolNOX   = 1.e-3;

  // multi grid parameters
  int maxGrids = 10;
  int numCycles = 2;

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
  my_CLP.setOption("withprec", &withprec, "0: no preconditioner, 1: block triangular, 2: block traingular with inner MG multigrid, 3: block triangular with mg solver");
  my_CLP.setOption("maxIter", &maxIter, "bla");
  my_CLP.setOption("initZero", &initZero, "bla");
  my_CLP.setOption("tolBelos", &tolBelos, "bla");
  my_CLP.setOption("tolInnerBelos", &tolInnerBelos, "bla");
  my_CLP.setOption("tolNOX", &tolNOX, "bla");

	// multi grid parameters
  my_CLP.setOption("maxGrids", &maxGrids, "bla");
  my_CLP.setOption("numCycles", &numCycles, "bla");


  my_CLP.recogniseAllOptions(true);
  my_CLP.throwExceptions(true);

  my_CLP.parse(argi,argv);


	// Set all the valid parameters and their default values.
	Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList("Pimpact");


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

typedef Pimpact::CoarsenStrategyGlobal<FSpaceT,CSpaceT,5> CS;

typedef Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> > VF;
typedef Pimpact::MultiHarmonicField< Pimpact::ScalarField<SpaceT> > SF;

typedef Pimpact::MultiField<VF> MVF;
typedef Pimpact::MultiField<SF> MSF;

typedef Pimpact::CompoundField< VF, SF> CF;
typedef Pimpact::MultiField<CF> MF;

typedef Pimpact::OperatorBase<MF> BOp;

typedef NOX::Pimpact::Interface<MF> Inter;
typedef NOX::Pimpact::Vector<typename Inter::Field> NV;


template<class T1,class T2> using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T> using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionOp<T> >;
template<class T> using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;

template<class T> using ConvDiffOpT = Pimpact::ConvectionVOp<Pimpact::ConvectionDiffusionSOp<T> >;
template<class T> using ConvDiffSORT = Pimpact::ConvectionVSmoother<T,Pimpact::ConvectionDiffusionSORSmoother >;


template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;






int main(int argi, char** argv ) {


  // intialize MPI
  MPI_Init( &argi, &argv );


	Teuchos::RCP<Teuchos::ParameterList> pl = getSpaceParametersFromCL<S,O>( argi, argv  ) ;

	pl->print();

  // solver stuff
  std::string nonLinSolName = pl->sublist("Solver").get<std::string>("nonLinSolver");
  std::string lineSearchName = pl->sublist("Solver").get<std::string>( "linSearch" );

	S tolNOX   = pl->sublist("Solver").get<S>("tolNOX");
  int maxIter = pl->sublist("Solver").get<int>("maxIter");
  int initZero = pl->sublist("Solver").get<int>("initZero");
  S tolBelos = pl->sublist("Solver").get<S>("tolBelos");
  S tolInnerBelos = pl->sublist("Solver").get<S>("tolInnerBelos");

	int withprec=pl->sublist("Solver").get<int>("withprec");

  std::string linSolName  = pl->sublist("Solver").get<std::string>( "linSolver" );
	if( linSolName.size()< 2 )
		linSolName = (withprec>0)?"Block GMRES":"GMRES";

	auto space = Pimpact::createSpace<S,O,3,4>( Teuchos::rcpFromRef( pl->sublist("Space", true) ) );
//  auto space = Pimpact::createSpace<S,O,3,4>( pl );

	int baseflow = pl->get<int>("baseflow");
	int flow = pl->get<int>("flow");
	int force = pl->get<int>("forcing");

  // outputs
  Teuchos::RCP<std::ostream> outPar;
  Teuchos::RCP<std::ostream> outLinSolve;
  Teuchos::RCP<std::ostream> outPrec;
  //  Teuchos::RCP<std::ostream> outLap2;
	 Teuchos::RCP<std::ostream> outSchur;

  if( space->rankST()==0 ) {
    outPar   = Teuchos::rcp( new std::ofstream("para_case.txt") );
    outLinSolve  = Teuchos::rcp( new std::ofstream("stats_linSolve.txt") );
    outPrec  = Teuchos::rcp( new std::ofstream("stats_DtConvDif.txt") );
    outSchur     = Teuchos::rcp( new std::ofstream("stats_DivGrad.txt") );
  } else
    outPar = Teuchos::rcp( new Teuchos::oblackholestream() ) ;


  if( space->rankST()==0 ) {
    Teuchos::rcp_static_cast<std::ofstream>(outPar)->close();
  }
  outPar = Teuchos::null;

	// init vectors
	auto x = Pimpact::createMultiField(
			Pimpact::createCompoundField(
				Pimpact::createMultiHarmonic< Pimpact::VectorField<SpaceT> >( space ),
				Pimpact::createMultiHarmonic< Pimpact::ScalarField<SpaceT> >( space )) );
  auto fu   = x->clone();

  // init Fields, init and rhs
	x->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->initField( Pimpact::EVectorField(baseflow), 1. );
//		if( 0 != flow ) {
//			x->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->initField( Pimpact::EVectorField(flow), 0., 1., 0., 0.1 );
//			x->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->getFieldPtr(Pimpact::U)->initField();
//		}

	auto sol = x->getFieldPtr(0)->getVFieldPtr()->clone( Pimpact::DeepCopy );
	//x->write(700);

	if( 0==initZero )
		x->init(0.);
	else if( 1==initZero ) {
		x->getFieldPtr(0)->getSFieldPtr()->init(0.);
		x->getFieldPtr(0)->getVFieldPtr()->random();
//		x->random();
		S re = space->getDomain()->getDomainSize()->getRe();
		x->scale(1.e-16);
	}
	if( 0==force )
		fu->init( 0. );
	else {
		S re = space->getDomain()->getDomainSize()->getRe();
		fu->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->getFieldPtr(Pimpact::U)->initField( Pimpact::FPoint,  -2. );
		fu->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->getFieldPtr(Pimpact::V)->initField( Pimpact::FPoint, 1. );
		fu->getFieldPtr(0)->getVFieldPtr()->getSFieldPtr(0)->getFieldPtr(Pimpact::W)->initField( Pimpact::FPoint, 1. );
	}
//	fu->getFieldPtr(0)->getVFieldPtr()->write( 700,true );


	//tolBelos*=l1*l2/n1/n2*(nfe-1)/nfs;

  /******1***********************************************************************************/
	{

		auto para = Pimpact::createLinSolverParameter( linSolName, tolBelos, -1, outLinSolve );
		para->set( "Num Blocks",         10  );
		para->set( "Maximum Iterations", 100 );
		para->set( "Maximum Restarts",   5	 );

    para->set( "Timer Label",	"Linear Problem");
		para->set( "Verbosity",	
				Belos::Errors +
				Belos::Warnings +
				Belos::IterationDetails +
				Belos::OrthoDetails +
				Belos::FinalSummary +
				Belos::StatusTestDetails +
				Belos::Debug
				);


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
		if( withprec>0 ) {

			// create Mulit space
			auto mgSpaces = Pimpact::createMGSpaces<FSpaceT,CSpaceT,CS>( space, pl->sublist("Multi Grid").get<int>("maxGrids") );

			// create Hinv
			auto v2v_para = Pimpact::createLinSolverParameter( (2==withprec)?"Block GMRES":"GMRES", tolInnerBelos, -1, outPrec );
			v2v_para->set( "Num Blocks",		     10  );
			v2v_para->set( "Maximum Iterations", 100 );
			v2v_para->set( "Maximum Restarts",	 10  );
			v2v_para->set( "Timer Label",	"F_inv");
			auto opV2Vprob =
					Pimpact::createLinearProblem<MVF>(
							Pimpact::createMultiOperatorBase(
									opV2V ),
									Teuchos::null,
									Teuchos::null,
									v2v_para,
									(2==withprec)?"Block GMRES":"GMRES" );

			// creat Hinv prec
      auto pls = Teuchos::parameterList();
			pls->set( "numCycles", pl->sublist("Multi Grid").get<int>("numCycles") );
			pls->sublist("Smoother").set( "omega", 1. );
			pls->sublist("Smoother").set( "numIters", 1 );
      pls->sublist("Smoother").set<int>( "Ordering", 1 );

			auto mgConvDiff =
						Pimpact::createMultiGrid<
						Pimpact::VectorField,
						TransVF,
						RestrVF,
						InterVF,
						ConvDiffOpT,
						ConvDiffOpT,
						ConvDiffSORT,
						MOP > ( mgSpaces, pls ) ;

			auto zeroOp = Pimpact::create<ConvDiffOpT>( space );
			auto zeroInv = Pimpact::create<MOP>( zeroOp );

			auto bla = Pimpact::createLinSolverParameter( "GMRES", tolInnerBelos/10., -1, Teuchos::rcp( new Teuchos::oblackholestream() ) );
//			auto bla = Pimpact::createLinSolverParameter( "GMRES", 1.e-6, -1, Teuchos::rcp( new Teuchos::oblackholestream() ) );
//			auto bla = Pimpact::createLinSolverParameter( "GMRES", tolInnerBelos/10., -1, Teuchos::rcp( &std::cout, false  ) );
			bla->set( "Num Blocks",				 	10   );
			bla->set( "Maximum Iterations", 100 );
			bla->set( "Maximum Restarts",	  10  );
			bla->set( "Timer Label",	"F_zero_inv");

			zeroInv->getOperatorPtr()->getLinearProblem()->setParameters( bla );

			zeroInv->getOperatorPtr()->getLinearProblem()->setRightPrec( Pimpact::createMultiOperatorBase(mgConvDiff) );

			Teuchos::RCP<Pimpact::OperatorBase<Pimpact::MultiField<Pimpact::MultiHarmonicField<Pimpact::VectorField<SpaceT> > > > >
				opV2Vprec = Teuchos::null;
				opV2Vprec = 
					Pimpact::createMultiOperatorBase(
							Pimpact::createMultiHarmonicDiagOp(zeroInv) );
			
			if( 2==withprec )
				opV2Vprob->setRightPrec( opV2Vprec );
			
			auto opV2Vinv =
				Pimpact::createInverseOperatorBase( opV2Vprob );


			// schurcomplement approximator ...
			auto opSchur =
					Pimpact::createMultiOperatorBase(
							Pimpact::createTripleCompositionOp(
									opV2S,
									opV2V,
									opS2V
							)
					);
			////--- inverse DivGrad
			auto divGradOp =
				Pimpact::createMultiOperatorBase(
//						Pimpact::createMultiHarmonicOpWrap(
							Pimpact::createCompositionOp(
								opV2S->getOperatorPtr(),
								opS2V->getOperatorPtr()
								)
//							)
						);

			auto pl_divGrad =
				Pimpact::createLinSolverParameter( (withprec==2||withprec==3)?"Block GMRES":"GMRES", tolInnerBelos/100., -1, outSchur );
			//				Pimpact::createLinSolverParameter( (withprec==2||withprec==3)?"Block GMRES":"GMRES", tolInnerBelos, -1, outSchur );
			pl_divGrad->set( "Timer Label",	"DivGrad");
			//			pl_divGrad->set( "Block Size", std::max( space->nGlo(3)/2, 1) );
			//			pl_divGrad->set( "Block Size", space->nGlo(3) );
			pl_divGrad->set( "Block Size", space->nGlo(3)*2+1 );
			pl_divGrad->set( "Adaptive Block Size", false );		
			//			pl_divGrad->set( "Block Size", 1 );

			auto divGradProb =
					Pimpact::createLinearProblem<Pimpact::MultiField<Pimpact::ScalarField<SpaceT> > >(
							divGradOp ,
							Teuchos::null,
							Teuchos::null,
							pl_divGrad,
							(withprec==2||withprec==3)?"Block GMRES":"GMRES" );

			/// init multigrid divgrad
			auto plmg = Teuchos::parameterList("MultiGrid");
			plmg->set<int>("numCycles",pl->sublist("Multi Grid").get<int>("numCycles"));
//			plmg->sublist("Smoother").set<S>("omega",0.8);
//			plmg->sublist("Smoother").set<int>("numIters",4);
		
			auto mg_divGrad =
				Pimpact::createMultiOperatorBase(
//						Pimpact::createMultiHarmonicOpWrap(
							Pimpact::createMultiGrid<
								Pimpact::ScalarField,
								Pimpact::TransferOp,
								Pimpact::RestrictionOp,
								Pimpact::InterpolationOp,
								Pimpact::DivGradOp,
								Pimpact::DivGradO2Op,
								Pimpact::DivGradO2JSmoother,
								MOP>( mgSpaces, plmg )
//									)
						);

			if( 2==withprec or 3==withprec )
				divGradProb->setRightPrec( mg_divGrad );

			auto divGradInv =
				Pimpact::createMultiOperatorBase(
						Pimpact::createMultiHarmonicMultiOpWrap( // fuse
//						Pimpact::createMultiHarmonicOpWrap( // fuse
//							Pimpact::createMultiOpUnWrap(
								Pimpact::createInverseOperatorBase( divGradProb )
								)
//							)
						);

			auto opS2Sinv =
					Pimpact::createOperatorBase(
							Pimpact::createTripleCompositionOp(
									divGradInv,
									opSchur,
									divGradInv
							)
					);

			auto invTriangOp = Pimpact::createInverseTriangularOp(
					opV2Vinv,
					opS2V,
					opS2Sinv );

			lprec =
				Pimpact::createMultiOperatorBase(
						invTriangOp );
		}


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

		if( 0==space->rankST() ) std::cout << "\n\t--- Nf: 0\tdof: "<<x->getLength(true)<<"\t---\n";

    // Solve the nonlinear system
		{
			Teuchos::TimeMonitor LocalTimer(*Teuchos::TimeMonitor::getNewCounter("Pimpact:: Solving Time"));
			solver->solve();
		}

		Teuchos::TimeMonitor::summarize();

    // Get the answer
    *group = solver->getSolutionGroup();

    x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr();
    Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr()->write( 800 );
    Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr()->getFieldPtr(0)->getVFieldPtr()->write( 400, true );
//		Teuchos::rcp_dynamic_cast<const NV>( group->getFPtr() )->getConstFieldPtr()->write(900);
  }
	/******************************************************************************************/
//	x->write(800);
//	auto er = x->clone(Pimpact::ShallowCopy);
	sol->add( 1., *sol, -1, *x->getFieldPtr(0)->getVFieldPtr());

	auto er = sol->norm();
	if( 0==space->rankST() )
		std::cout << "error: " << er << "\n";


  MPI_Finalize();
  return( 0 );

}
