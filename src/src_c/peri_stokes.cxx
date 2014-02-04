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
#include "Pimpact_GridSize.hpp"
#include "Pimpact_ProcGridSize.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiField.hpp"
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

	Pimpact::init_impact_post();

	// init Spaces
	auto fS = Pimpact::createFieldSpace<int>();

	auto iIS = Pimpact::createInnerFieldIndexSpaces<int>();
	auto fIS = Pimpact::createFullFieldIndexSpaces<int>();

	// init vectors
	auto scac = Pimpact::createScalarField<double,int>(fS);
	auto scas = Pimpact::createScalarField<double,int>(fS);

	auto sca = Pimpact::createModeField( scac, scas );

	auto p     = Pimpact::createMultiField< SF, double, int >(*sca,1);
	auto temps = Pimpact::createMultiField< SF, double, int >(*sca,1);

	scac=Teuchos::null;
	scas=Teuchos::null;
	sca=Teuchos::null;
//	vel=Teuchos::null;
//	velc=Teuchos::null;
//	vels=Teuchos::null;

	p->Init(0.);

	/// \todo move this to a factory
	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );

	auto u     = Pimpact::createMultiField< VF, double, int >(*vel,1);
	auto f     = Pimpact::createMultiField< VF, double, int >(*vel,1);
	auto tempv = Pimpact::createMultiField< VF, double, int >(*vel,1);

	switch( Pimpact::EFlowType(flow) ) {
	case Pimpact::ZeroFLow :
		u->GetVec(0).getFieldC()->init_field( Pimpact::ZeroProf );
		u->GetVec(0).getFieldS()->init_field( Pimpact::ZeroProf );

//		tempv->GetVec(0).getFieldC()->init_field( Pimpact::ZeroProf );
//		tempv->GetVec(0).getFieldS()->init_field( Pimpact::ZeroProf );
		break;
	case Pimpact::Poiseuille_inX :
		u->GetVec(0).getFieldC()->init_field( Pimpact::Poiseuille2D_inX );
		u->GetVec(0).getFieldS()->init_field( Pimpact::ZeroProf );

//		tempv->GetVec(0).getFieldC()->init_field( Pimpact::Poiseuille2D_inX );
//		tempv->GetVec(0).getFieldS()->init_field( Pimpact::ZeroProf );
		break;
	case Pimpact::Poiseuille_inY :
		u->GetVec(0).getFieldC()->init_field( Pimpact::Poiseuille2D_inY );
		u->GetVec(0).getFieldS()->init_field( Pimpact::ZeroProf );

//		tempv->GetVec(0).getFieldC()->init_field( Pimpact::Poiseuille2D_inY );
//		tempv->GetVec(0).getFieldS()->init_field( Pimpact::ZeroProf );
		break;
	case Pimpact::Pulsatile_inX :
		u->GetVec(0).getFieldC()->init_field( Pimpact::Pulsatile2D_inXC, re, omega, px );
		u->GetVec(0).getFieldS()->init_field( Pimpact::Pulsatile2D_inXS, re, omega, px );

//		tempv->GetVec(0).getFieldC()->init_field( Pimpact::Pulsatile2D_inXC, re, omega, px );
//		tempv->GetVec(0).getFieldS()->init_field( Pimpact::Pulsatile2D_inXS, re, omega, px );
		break;
	case Pimpact::Pulsatile_inY :
		u->GetVec(0).getFieldC()->init_field( Pimpact::Pulsatile2D_inYC, re, omega, px );
		u->GetVec(0).getFieldS()->init_field( Pimpact::Pulsatile2D_inYS, re, omega, px );

//		tempv->GetVec(0).getFieldC()->init_field( Pimpact::Pulsatile2D_inYC, re, omega, px );
//		tempv->GetVec(0).getFieldS()->init_field( Pimpact::Pulsatile2D_inYS, re, omega, px );
		break;
	case Pimpact::Streaming2DFlow :
		u->GetVec(0).getFieldC()->init_field( Pimpact::Streaming2D, re, omega, px );
		u->GetVec(0).getFieldS()->init_field( Pimpact::ZeroProf, re, omega, px );

//		tempv->GetVec(0).getFieldC()->init_field( Pimpact::Streaming2D, re, omega, px );
//		tempv->GetVec(0).getFieldS()->init_field( Pimpact::ZeroProf, re, omega, px );
		break;
	}

	u->Init(0);
	tempv->Init(0);

	f->GetVec(0).getFieldC()->init_field( Pimpact::ZeroProf );
	f->GetVec(0).getFieldS()->init_field( Pimpact::ZeroProf );

	// init operators
	auto lap  = Pimpact::createDtL<Scalar,Ordinal>( omega, 0., 1./re );
	auto div  = Pimpact::createOperatorMV<Pimpact::Div<double,int> >();
	auto grad = Pimpact::createOperatorMV<Pimpact::Grad<double,int> >();


//	RCP<ParameterList> solverParams = parameterList();
	auto solverParams = Pimpact::createLinSolverParameter( solver_name_1, 1.e-6*l1*l2/n1/n2 );

//	solverParams->set ("Num Blocks", 100);
//	solverParams->set ("Maximum Iterations", 2000);
//	solverParams->set ("Maximum Restarts", 20);
//	solverParams->set ("Convergence Tolerance", 1.0e-6*l1*l2/n1/n2);
//	solverParams->set ("Implicit Residual Scaling","None");
//	solverParams->set ("Explicit Residual Scaling","None");
//	solverParams->set ("Output Frequency", 10 );
//	solverParams->set ("Output Style", 1 );
	solverParams->get()->set ("Output Stream", outLap1 );
//	solverParams->set ("Verbosity",  int(Belos::Errors + Belos::Warnings + Belos::IterationDetails +
//			Belos::OrthoDetails + Belos::FinalSummary +	Belos::TimingDetails +
//			Belos::StatusTestDetails + Belos::Debug) );
//
	Teuchos::writeParameterListToXmlFile( *solverParams->get(), "para_solver.xml" );

//	solverParams = Teuchos::getParametersFromXmlFile("solver1.xml");

	// create problems/solvers
	auto lap_problem = Pimpact::createLinearProblem<Scalar,MVF,Lap>( lap, u, f, solverParams->get(), solver_name_1 );

	auto schur = Pimpact::createDivDtLinvGrad<double,int>( u, lap_problem );

//	solverParams = parameterList();
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

	tempv->Add( -1., *tempv, 1., *f );

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
	return 0;
}
