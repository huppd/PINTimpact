#include <ostream>
#include <fstream>
#include <mpi.h>

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
#include "Pimpact_GridSize.hpp"
#include "Pimpact_ProcGridSize.hpp"
//#include "Pimpact_FieldSpace.hpp"
//#include "Pimpact_IndexSpace.hpp"
//#include "Pimpact_ScalarField.hpp"
//#include "Pimpact_VectorField.hpp"
//#include "Pimpact_ModeField.hpp"
//#include "Pimpact_MultiField.hpp"
#include "Pimpact_FieldFactory.hpp"
#include "Pimpact_OperatorMV.hpp"
#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "BelosPimpactAdapter.hpp"




int main(int argi, char** argv ) {

	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing

	typedef double Scalar;
	typedef int Ordinal;
	typedef Pimpact::ModeField<Pimpact::VectorField<Scalar,Ordinal> > VF;
	typedef Pimpact::ModeField<Pimpact::ScalarField<Scalar,Ordinal> >  SF;
	typedef Pimpact::MultiField<VF> MVF;
	typedef Pimpact::MultiField<SF> MSF;
	typedef Pimpact::OperatorMV< Pimpact::DtL<Scalar,Ordinal> >  Lap;
	typedef Pimpact::OperatorMV< Pimpact::Div_DtLinv_Grad<Scalar,Ordinal> >  Schur;
	typedef Pimpact::OperatorMV< Pimpact::Grad<Scalar,Ordinal> >  G;

	// intialize MPI
	MPI_Init( &argi, &argv );

	//// get problem values form Comand line
	Teuchos::CommandLineProcessor my_CLP;

	// physical constants
	Scalar re = 1.;
	my_CLP.setOption( "re", &re, "Reynolds number" );

	Scalar omega = 1.;
	my_CLP.setOption( "omega", &omega, "introduced frequency" );

	Scalar px = 1.;
	my_CLP.setOption( "px", &px, "pressure gradient(only necessary for pulsatile flows)" );

	// flow type
	int flow = 1;
	my_CLP.setOption( "flow", &flow,
			"Flow type: 0=zero flow, 1=2D Poiseuille flow in x, 2=2D Poiseuille flow in y, 3=2D pulsatile flow in x, 4=2D pulsatile flow in, 5=2D streaming" );

	// domain type
	int domain = 1;
	my_CLP.setOption( "domain", &domain,
			"Domain type: 0:all dirichlet, 1:dirichlet 2d channel, 2: periodic 2d channel" );

	// domain size
	Scalar l1 = 2.;
	my_CLP.setOption( "lx", &l1, "length in x-direction" );

	Scalar l2 = 2.;
	my_CLP.setOption( "ly", &l2, "length in y-direction" );

	Scalar l3 = 2.;
	my_CLP.setOption( "lz", &l3, "length in z-direction" );

	// grid size
	Ordinal n1 = 33;
	my_CLP.setOption( "nx", &n1, "amount of grid points in x-direction: a*2**q+1" );

	Ordinal n2 = 33;
	my_CLP.setOption( "ny", &n2, "amount of grid points in y-direction: a*2**q+1" );

	Ordinal n3 = 2.;
	my_CLP.setOption( "nz", &n3, "amount of grid points in z-direction: a*2**q+1" );

	// processor grid size
	Ordinal np1 = 2;
	my_CLP.setOption( "npx", &np1, "amount of processors in x-direction" );

	Ordinal np2 = 2;
	my_CLP.setOption( "npy", &np2, "amount of processors in y-direction" );

	Ordinal np3 = 1.;
	my_CLP.setOption( "npz", &np3, "amount of processors in z-direction" );

	// solver name
	std::string solver_name_1 = "GMRES";
	my_CLP.setOption( "solver1", &solver_name_1, "name of the solver for H" );

	std::string solver_name_2 = "GMRES";
	my_CLP.setOption( "solver2", &solver_name_2, "name of the solver for Schur complement" );

	my_CLP.recogniseAllOptions(true);
	my_CLP.throwExceptions(true);

	my_CLP.parse(argi,argv);
	// end of parsing

	// starting with ininializing
	int rank = Pimpact::init_impact_pre();

	// outputs
	Teuchos::RCP<std::ostream> outPar;
	Teuchos::RCP<std::ostream> outLap1;
	Teuchos::RCP<std::ostream> outLap2;
	Teuchos::RCP<std::ostream> outSchur;

	if(rank==0) {
		outPar   = Teuchos::rcp( new std::ofstream("para_case.txt") );
		outLap1  = Teuchos::rcp( new std::ofstream("stats_solvLap1.txt") );
		outLap2  = Teuchos::rcp( new std::ofstream("stats_solvLap2.txt") );
		outSchur = Teuchos::rcp( new std::ofstream("stats_solvSchur.txt") );
	}	else
//		outPar = Teuchos::rcp( &blackhole, false) ;
		outPar = Teuchos::rcp( new Teuchos::oblackholestream() ) ;

	*outPar << " \tflow=" << flow << "\n";
	*outPar << " \tdomain" << domain << "\n";
	*outPar << " \tre=" << re << "\n";
	*outPar << " \tpx=" << px << "\n";
	*outPar << " \tomega=" << omega << "\n";

	auto ds = Pimpact::createDomainSize<Scalar>(l1,l2,l3);
	ds->set_Impact();
	ds->print( *outPar );

	auto bc = Pimpact::createBC( Pimpact::EDomainType(domain) );
	bc->set_Impact();

	auto gs = Pimpact::createGridSize<Ordinal>(n1,n2,n3);
	gs->set_Impact();
	gs->print( *outPar );

	auto pgs = Pimpact::createProcGridSize<Ordinal>(np1,np2,np3);
	pgs->set_Impact();
	pgs->print( *outPar );

	if(rank==0) {
		Teuchos::rcp_static_cast<std::ofstream>(outPar)->close();
	}
	outPar = Teuchos::null;


	// init IMPACT
	Pimpact::init_impact_post();


	// init Spaces
	auto fS = Pimpact::createFieldSpace<Ordinal>();

	auto iIS = Pimpact::createInnerFieldIndexSpaces<Ordinal>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<Ordinal>();


	// init vectors
	auto p     = Pimpact::createInitMSF<Scalar,Ordinal>( fS );
	auto temps = Pimpact::createInitMSF<Scalar,Ordinal>( fS );

	auto u     = Pimpact::createInitMVF<Scalar,Ordinal>( Pimpact::EFlowType(flow), fS, iIS, fIS );
	auto tempv = Pimpact::createInitMVF<Scalar,Ordinal>( Pimpact::EFlowType(flow), fS, iIS, fIS );
	auto f     = Pimpact::createInitMVF<Scalar,Ordinal>( Pimpact::EFlowType(flow), fS, iIS, fIS );

	p->init(0.);
	u->init(0);
	tempv->init(0);


	// init operators
	auto lap  = Pimpact::createDtL<Scalar,Ordinal>( omega, 0., 1./re );
	auto div  = Pimpact::createOperatorMV<Pimpact::Div<double,int> >();
	auto grad = Pimpact::createOperatorMV<Pimpact::Grad<double,int> >();

	// create parameter for linsovlers
	auto solverParams = Pimpact::createLinSolverParameter( solver_name_1, 1.e-6*l1*l2/n1/n2 );

	solverParams->get()->set ("Output Stream", outLap1 );
//
	Teuchos::writeParameterListToXmlFile( *solverParams->get(), "para_solver.xml" );


	// create preconditioner
//	auto lprec = Pimpact::createDtL<Scalar,Ordinal>( -1./omega, 0., 0. );

	// create problems/solvers
	auto lap_problem = Pimpact::createLinearProblem<Scalar,MVF,Lap>( lap, u, f, solverParams->get(), solver_name_1 );
//	lap_problem->setLeftPrec( lprec );

	auto schur = Pimpact::createDivDtLinvGrad<double,int>( u, lap_problem );

	solverParams = Pimpact::createLinSolverParameter( solver_name_2, 1.e-6*l1*l2/n1/n2  );
	solverParams->get()->set( "Output Stream", outSchur );
	auto schur_prob = Pimpact::createLinearProblem<Scalar,MSF,Schur>( schur, p,temps, solverParams->get(), solver_name_2);

	// solve stationary stokes
	lap_problem->solve( tempv, f );
	//  tempv->write( 2000 );

	div->apply( *tempv, *temps );
	//temps->write( 2002 );

	solverParams = Pimpact::createLinSolverParameter( solver_name_1, 1.e-6*l1*l2/n1/n2 );
	solverParams->get()->set( "Output Stream", outLap2 );
	solverParams->get()->set ("Verbosity", int( Belos::Errors) );
	lap_problem->setParameters( solverParams->get() );
	schur_prob->solve( p, temps );
	p->write();

	grad->apply( *p, *tempv );
	//tempv->write(2006);

	tempv->add( -1., *tempv, 1., *f );

	solverParams->get()->set ("Verbosity",  Belos::Errors + Belos::Warnings + Belos::IterationDetails +
			Belos::OrthoDetails + Belos::FinalSummary +	Belos::TimingDetails +
			Belos::StatusTestDetails + Belos::Debug );
	solverParams->get()->set( "Output Stream", outLap2 );
	lap_problem->setParameters( solverParams->get() );

	lap_problem->solve(u,tempv);

	u->write();

	if(rank==0) {
		Teuchos::rcp_static_cast<std::ofstream>( outLap1)->close();
		Teuchos::rcp_static_cast<std::ofstream>( outLap2)->close();
		Teuchos::rcp_static_cast<std::ofstream>(outSchur)->close();
	}


	MPI_Finalize();
	return( 0 );
}
