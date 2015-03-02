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
	S re = 100.;

	S alpha2 = 125.;

  // domain type
  int domain = 2;

  // domain size
  int dim = 2;
  S lx = 1.;
  S ly = 1.;
  S lz = 1.;

  // grid size
	O nx = 33;
	O ny = 33;
	O nz = 2;
	O nf = 4;

	O nfs = 4;
	O nfe = 4;

  // processor grid size
  O npx = 1;
  O npy = 1;
  O npz = 1.;
  O npf = 1.;

	// flow and forcing type
  int flow = 7;
  int forcing = 0;

	// solver parameters
  std::string nonLinSolName = "Newton";

  std::string lineSearchName = "Backtrack";

  std::string linSolName = "GMRES";

	int withprec=0;

  int maxIter = 10;

  S tolBelos = 1.e-1;
  S tolInnerBelos = 1.e-3;
  S tolNOX   = 1.e-3;

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

	// flow and forcint
  my_CLP.setOption( "flow", &flow,
      "Flow type" );
  my_CLP.setOption( "force", &forcing,
      "forcing, ja?" );

	// solver parameters
  my_CLP.setOption("nonLinSolver", &nonLinSolName, "Newton");
  my_CLP.setOption("linSearch", &lineSearchName, "Backtrack");
  my_CLP.setOption("linSolver", &linSolName, "bla");
  my_CLP.setOption("withprec", &withprec, "withprec");
  my_CLP.setOption("maxIter", &maxIter, "bla");
  my_CLP.setOption("tolBelos", &tolBelos, "bla");
  my_CLP.setOption("tolInnerBelos", &tolInnerBelos, "bla");
  my_CLP.setOption("tolNOX", &tolNOX, "bla");


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
	pl->sublist("Space").set<O>("npy", npx, "amount of processors in y-direction" );
	pl->sublist("Space").set<O>("npz", npz, "amount of processors in z-direction" );
	pl->sublist("Space").set<O>("npf", npf, "amount of processors in f-direction" );

	// flow and forcint
	pl->set<int>( "flow", flow, "Flow type: depending on main" );
	pl->set<int>( "forcing", forcing, "forcing, ja?" );

	// solver parameters
  pl->sublist("Solver").set("nonLinSolver", nonLinSolName );
  pl->sublist("Solver").set("linSearch", lineSearchName );
  pl->sublist("Solver").set("linSolver", linSolName );
  pl->sublist("Solver").set("withprec", withprec );
  pl->sublist("Solver").set("maxIter", maxIter );
  pl->sublist("Solver").set("tolBelos", tolBelos );
  pl->sublist("Solver").set("tolInnerBelos", tolInnerBelos );
  pl->sublist("Solver").set("tolNOX", tolNOX );

	return( pl );
}


typedef double S;
typedef int O;

typedef Pimpact::Space<S,O,3,4> SpaceT;

typedef Pimpact::Space<S,O,3,4> FSpaceT;
typedef Pimpact::Space<S,O,3,2> CSpaceT;

typedef Pimpact::CoarsenStrategy<FSpaceT,CSpaceT> CS;

typedef Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> > VF;
typedef Pimpact::MultiHarmonicField< Pimpact::ScalarField<SpaceT> > SF;

typedef Pimpact::MultiField<VF> MVF;
typedef Pimpact::MultiField<SF> MSF;

typedef Pimpact::CompoundField< VF, SF> CF;
typedef Pimpact::MultiField<CF> MF;

typedef Pimpact::OperatorBase<MF> BOp;

typedef NOX::Pimpact::Interface<MF> Inter;
typedef NOX::Pimpact::Vector<typename Inter::Field> NV;


template<class T> using ConvDiffOpT = Pimpact::ConvectionVOp<Pimpact::ConvectionDiffusionSOp<T> >;

template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;


//template<class ST> using BOPF = Pimpact::MultiOpWrap< Pimpact::DivGradOp<ST> >;
//template<class ST> using BOPC = Pimpact::MultiOpWrap< Pimpact::DivGradO2Op<ST> >;
//template<class ST> using BSM = Pimpact::MultiOpWrap< Pimpact::DivGradO2JSmoother<ST> >;


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
  S tolBelos = pl->sublist("Solver").get<S>("tolBelos");
  S tolInnerBelos = pl->sublist("Solver").get<S>("tolInnerBelos");

	int withprec=pl->sublist("Solver").get<int>("withprec");

//  std::string linSolName = withprec?"Block GMRES":"GMRES";
//  std::string linSolName = false?"Block GMRES":"GMRES";
  std::string linSolName  = pl->sublist("Solver").get<std::string>( "linSolver" );

	auto space = Pimpact::createSpace<S,O,3,4>( Teuchos::rcpFromRef( pl->sublist("Space", true) ) );
//  auto space = Pimpact::createSpace<S,O,3,4>( pl );

	int flow = pl->get<int>("flow");
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
  auto x    = Pimpact::createMultiField(
			Pimpact::createCompoundField(
      	Pimpact::createMultiHarmonic< Pimpact::VectorField<SpaceT> >( space ),
      	Pimpact::createMultiHarmonic< Pimpact::ScalarField<SpaceT> >( space )) );
  auto fu   = x->clone();

  // init Fields, init and rhs
  x->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->initField( Pimpact::EFlowField(flow), 1 );

	x->init(0.);
  fu->init( 0. );


	//tolBelos*=l1*l2/n1/n2*(nfe-1)/nfs;

  /******************************************************************************************/
	{

    auto para = Pimpact::createLinSolverParameter( linSolName, tolBelos, -1 );
		para->set( "Maximum Iterations", 1000 );
		para->set( "Implicit Residual Scaling", "Norm of RHS" );
		para->set( "Explicit Residual Scaling", "Norm of RHS" );
		para->set( "Flexible Gmres", true );
		para->set( "Output Frequency", withprec?1:100 );
		para->set( "Verbosity",	
				Belos::Errors +
				Belos::Warnings +
				Belos::IterationDetails +
				Belos::OrthoDetails +
				Belos::FinalSummary +
				Belos::TimingDetails +
				Belos::StatusTestDetails +
				Belos::Debug
				);


    auto opV2V =
				Pimpact::createMultiDtConvectionDiffusionOp( space );
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
		if( withprec ) {

			// create Mulit space
			auto mgSpaces = Pimpact::createMGSpaces<FSpaceT,CSpaceT,CS>( space, 3 );

			// create Hinv
			auto opV2Vprob =
					Pimpact::createLinearProblem<MVF>(
							Pimpact::createMultiOperatorBase(
									opV2V ),
									Teuchos::null,
									Teuchos::null,
									Pimpact::createLinSolverParameter( "Block GMRES", tolInnerBelos, 100, outPrec ), "Block GMRES" );
//										Pimpact::createLinSolverParameter( "GMRES", tolInnerBelos, 100, outPrec ), "GMRES" );

			// creat Hinv prec
			auto zeroOp = Pimpact::create<ConvDiffOpT>( space );
			auto zeroInv = Pimpact::create<MOP>( zeroOp );

			auto opV2Vprec = Pimpact::createMultiOperatorBase(
					Pimpact::createMultiHarmonicDiagOp(zeroInv) );
			
			opV2Vprob->setRightPrec( opV2Vprec );
			
			auto opV2Vinv =
				Pimpact::createInverseOperatorBase( opV2Vprob );

//			opV2Vinv = opV2Vprec;
//			auto eye = x->getFieldPtr(0)->getVFieldPtr()->clone();
//			for( O i=0; i<nf; ++i ) {
//				eye->getCFieldPtr(i)->initField( Pimpact::ConstFlow, 1., 1., 1. );
//				eye->getSFieldPtr(i)->initField( Pimpact::ConstFlow, 1., 1., 1. );
//			}
//
//			auto opV2Vinv = 
//				Pimpact::createMultiOperatorBase(
//						Pimpact::createForcingOp( eye ) );

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
						Pimpact::createCompositionOp(
							opV2S,
							opS2V
							)
						);


			auto divGradProb =
					Pimpact::createLinearProblem<MSF>(
							divGradOp ,
							Teuchos::null,
							Teuchos::null,
//							Pimpact::createLinSolverParameter( "Block GMRES", tolInnerBelos, -1, outSchur ), "Block GMRES" );
							Pimpact::createLinSolverParameter( "GMRES", tolInnerBelos, -1, outSchur ), "GMRES" );
//							Pimpact::createLinSolverParameter( "CG", tolInnerBelos, -1, outSchur ), "CG" );

			auto mg_divGrad =
				Pimpact::createMultiOperatorBase(
						Pimpact::createMultiHarmonicOpWrap(
							Pimpact::createMultiGrid<
								Pimpact::ScalarField,
								Pimpact::TransferOp,
								Pimpact::RestrictionOp,
								Pimpact::InterpolationOp,
								Pimpact::DivGradOp,
								Pimpact::DivGradO2Op,
								Pimpact::DivGradO2JSmoother,
								MOP>( mgSpaces ) ) );

//			divGradProb->setRightPrec( mg_divGrad );

			auto divGradInv = Pimpact::createInverseOperatorBase( divGradProb );
//			divGradInv = mg_divGrad;

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
    auto statusTest = NOX::Pimpact::createStatusTest( maxIter, tolNOX, tolBelos/100 );

    // Create the list of solver parameters
    auto solverParametersPtr =
        NOX::Pimpact::createNOXSolverParameter( nonLinSolName, lineSearchName );

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver =
        NOX::Solver::buildSolver( group, statusTest, solverParametersPtr);

    // Solve the nonlinear system
		solver->solve();


    // Get the answer
    *group = solver->getSolutionGroup();

    x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr();
  }
  /******************************************************************************************/
  x->write(800);

  MPI_Finalize();
  return( 0 );

}
