#include <iostream>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_VectorFieldOpWrap.hpp"
#include "Pimpact_MultiGrid.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Operator.hpp"

#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"

#include "Pimpact_InterpolationTimeOp.hpp"
#include "Pimpact_RestrictionTimeOp.hpp"
#include "Pimpact_CoarsenStrategy.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"


#include "Pimpact_IntResCompoundOp.hpp"
#include "Pimpact_TransferCompoundOp.hpp"
#include "Pimpact_TransferTimeOp.hpp"
#include "Pimpact_TimeStokesBSmoother.hpp"

namespace {

bool output = false;

using ST = double;
using OT = int;
const int sd = 3;

using SpaceT = Pimpact::Space<ST,OT,sd,4,4>; 

using FSpace3T = Pimpact::Space<ST,OT,sd,3,4>;
using FSpace4T = Pimpact::Space<ST,OT,sd,4,4>; 

using CSpace3T = Pimpact::Space<ST,OT,sd,3,2>;
using CSpace4T = Pimpact::Space<ST,OT,sd,4,2>; 

using CS3L = Pimpact::CoarsenStrategy<FSpace3T,CSpace3T>;
using CS3G = Pimpact::CoarsenStrategyGlobal<FSpace3T,CSpace3T,5>;

using CS4L = Pimpact::CoarsenStrategy<FSpace4T,CSpace4T>;
using CS4G = Pimpact::CoarsenStrategyGlobal<FSpace4T,CSpace4T,5>;

template<class ST> using BSF = Pimpact::MultiField< Pimpact::ScalarField<ST> >;
//template<class T> using BVF = Pimpact::MultiField< Pimpact::VectorField<T> >;

template<class ST> using TSF = Pimpact::TimeField< Pimpact::ScalarField<ST> >;
template<class ST> using TOPF = Pimpact::TimeOpWrap< Pimpact::DivGradOp<ST> >;
template<class ST> using TOPC = Pimpact::TimeOpWrap< Pimpact::DivGradO2Op<ST> >;


template<class ST> using BOPF = Pimpact::MultiOpWrap< Pimpact::DivGradOp<ST> >;
template<class ST> using BOPC = Pimpact::MultiOpWrap< Pimpact::DivGradO2Op<ST> >;
template<class ST> using BSM = Pimpact::MultiOpWrap< Pimpact::DivGradO2JSmoother<ST> >;

template<class T> using ConvDiffOpT = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<T> >;

//template<class T> using ConvDiffOpT = Pimpact::ConvectionVOp<Pimpact::ConvectionDiffusionSOp<T> >;

template<class T> using ConvDiffSORT = Pimpact::NonlinearSmoother<T,Pimpact::ConvectionDiffusionSORSmoother >;

template<class T> using ConvDiffJT = Pimpact::NonlinearSmoother<T,Pimpact::ConvectionDiffusionJSmoother >;

bool testMpi = true;
double eps = 1e-6;

int domain = 0;
int ftype = 0;

int fs = 0;
int fe = 4;

int npx = 1;
int npy = 1;
int npz = 1;
int npf = 1;

int nx = 17;
int ny = 17;
int nz = 17;
int nf = 16;

int rankbla = -1;

int maxGrids = 10;

auto pl = Teuchos::parameterList();


TEUCHOS_STATIC_SETUP() {
	Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
	clp.addOutputSetupOptions(true);
	clp.setOption(
			"test-mpi", "test-serial", &testMpi,
			"Test MPI (if available) or force test of serial.  In a serial build,"
			" this option is ignored and a serial comm is always used." );
	clp.setOption(
			"error-tol-slack", &eps,
			"Slack off of machine epsilon used to check test results" );
	clp.setOption(
			"domain", &domain,
			"Slack off of machine epsilon used to check test results" );
	clp.setOption(
			"ftype", &ftype,
			"Slack off of machine epsilon used to check test results" );
	clp.setOption(
			"fs", &fs,
			"Slack off of machine epsilon used to check test results" );
	clp.setOption(
			"fe", &fe,
			"Slack off of machine epsilon used to check test results" );
	clp.setOption( "npx", &npx, "" );
	clp.setOption( "npy", &npy, "" );
	clp.setOption( "npz", &npz, "" );
	clp.setOption( "npf", &npf, "" );
	clp.setOption( "nx", &nx, "" );
	clp.setOption( "ny", &ny, "" );
	clp.setOption( "nz", &nz, "" );
	clp.setOption( "nf", &nf, "" );
	clp.setOption( "rank", &rankbla, "" );
	clp.setOption( "maxGrids", &maxGrids, "" );
	clp.setOption(
			"output", "noutput", &output,
			"Test MPI (if available) or force test of serial.  In a serial build,"
			" this option is ignored and a serial comm is always used." );


	pl->set( "dim", 3 );
	Pimpact::setBoundaryConditions( pl, domain );

	pl->set( "lx", 2. );
	pl->set( "ly", 2. );
	pl->set( "lz", 2. );

	pl->set( "Re", 1. );
	pl->set( "alpha2", 1. );
}


template<class SpaceT> using CVF = Pimpact::CompoundField<Pimpact::TimeField<Pimpact::VectorField<SpaceT> >,
Pimpact::TimeField<Pimpact::ScalarField<SpaceT> > >;

template<class SpaceT> using INT = Pimpact::IntResCompoundOp<
	Pimpact::InterpolationTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<SpaceT> > >,
	Pimpact::InterpolationTimeOp<                           Pimpact::InterpolationOp<SpaceT> > >;

template<class SpaceT> using RES = Pimpact::IntResCompoundOp<
Pimpact::RestrictionTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<SpaceT> > >,
	Pimpact::RestrictionTimeOp<                           Pimpact::RestrictionSFOp<SpaceT> > >;

template<class SpaceT1, class SpaceT2> using TCO = Pimpact::TransferCompoundOp<
Pimpact::TransferTimeOp<Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<SpaceT1, SpaceT2> > >,
	Pimpact::TransferTimeOp<                           Pimpact::TransferOp<SpaceT1, SpaceT2> > >; 

template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, MG, CS ) {

	pl->set("nx", nx );
	pl->set("ny", ny );
	pl->set("nz", nz );
	pl->set("nf", nf );

	pl->set("npx", npx );
	pl->set("npy", npy );
	pl->set("npz", npz );
	pl->set("npf", npf );

	Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl ); 

	auto mgSpaces = Pimpact::createMGSpaces<FSpace4T,CSpace4T,CS>( space, maxGrids );


	auto mgPL = Teuchos::parameterList();
	mgPL->sublist("Smoother").set<std::string>("Solver name", "GMRES" );
	mgPL->sublist("Smoother").set( "Solver",
			*Pimpact::createLinSolverParameter( "GMRES", 1.e-16, -1,
				Teuchos::rcp<std::ostream>( new Teuchos::oblackholestream() ), 4 ) );

	mgPL->sublist("Smoother").set<int>( "numIters", 4 );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<int>( "Maximum Iterations", 1000 );
	mgPL->sublist("Coarse Grid Solver").set<std::string>("Solver name", "GMRES" );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<std::string>("Timer Label", "Coarse Grid Solver" );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<ST>("Convergence Tolerance" , 1.e-1 );

	auto mg = Pimpact::createMultiGrid<
		CVF,
		TCO,
		RES,
		INT,
		Pimpact::TimeNSOp,
		Pimpact::TimeNSOp,								
		Pimpact::TimeNS4DBSmoother,
		//									Pimpact::TimeStokesBSmoother
		MOP
			> ( mgSpaces, mgPL );

	//	mg->print();

	auto x = Pimpact::createCompoundField( Pimpact::createTimeField< Pimpact::VectorField<FSpace4T> >( space ),
			Pimpact::createTimeField< Pimpact::ScalarField<FSpace4T> >( space ));

	using OpT = Pimpact::TimeStokesOp<FSpace4T>;
	auto op = Pimpact::create<OpT>( space );

	auto op_true = Pimpact::create<OpT>( space );

	double p = 1;
	double alpha = std::sqrt(pl->get<double>("alpha2"));

	auto b = x->clone();
	auto b_bc = x->clone();
	auto Ax = x->clone();
	auto true_sol = x->clone();
	auto err = x->clone();
	err->init(1);

	Pimpact::initVectorTimeField( true_sol->getVFieldPtr(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha );

	op->assignField(*true_sol);

	op_true->assignField(*true_sol);

	//auto true_sol2 = true_sol->clone();

	op->apply(*true_sol,*b);

	// consistency  (is ok but in the corsest level x=0 initial guess)
	//mg->apply( *b, *true_sol2);
	//err->add(-1.,*true_sol,1.,*true_sol2);
	//std::cout << "\n" << "consistency err: " << err->norm()/std::sqrt( err->getLength() );

	//err->write();

	// put BC in the RHS
	true_sol->init(0);
	op->apply(*true_sol,*b_bc);

	b->add(1.,*b,-1.,*b_bc);

	/////////

	//consistency
	//true_sol2 = true_sol->clone();
	//mg->apply( *b, *true_sol2 );
	//err->add(-1.,*true_sol,1.,*true_sol2);
	//std::cout << "\n" << "consistency err with bc in RHS: " << err->norm()/std::sqrt( err->getLength() );

	///////

	Pimpact::initVectorTimeField( true_sol->getVFieldPtr(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha );

	// for consistency, ok
	//x = true_sol->clone();

	x->random();
	//x->getSField().init(0.);
	x->scale(10);
	x->getSField().level();

	//op->assignField(*x);

	err->add( -1, *x, 1., *true_sol );
	std::cout << "\n" << "err: " << err->norm()/std::sqrt( err->getLength() )<< std :: endl;

	if (output)
		err->write();

	op->apply(*x,*Ax);

	err->add( -1, *Ax, 1., *b );
	std::cout << err->norm()/std::sqrt( err->getLength() ) << "\n";

	std::cout << "Start MG" << std::endl;

	for( int i=0; i<50; ++i ) {
		mg->apply( *b, *x );

		x->level();

		err->add( -1, *x, 1., *true_sol );

		if (output)
			err->write((i+1)*100);

		std::cout << "err: " << err->norm()/std::sqrt( err->getLength() );

		op->apply(*x,*Ax);

		err->add( -1, *Ax, 1., *b );
		std::cout <<  "       res " << err->norm()/std::sqrt( err->getLength() ) << "\n";
	}
	//x->write();

	TEST_EQUALITY( err->norm()<1.e-5, true );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, MG, CS4L )
//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, MG, CS4G )

} // end of namespae
