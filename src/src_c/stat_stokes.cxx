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
#include "Pimpact_MultiField.hpp"
#include "Pimpact_OperatorMV.hpp"
#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_Operator.hpp"
#include "BelosPimpactAdapter.hpp"




int main(int argi, char** argv ) {

	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing

	typedef double Scalar;
	typedef int Ordinal;
	typedef Pimpact::VectorField<Scalar,Ordinal> VF;
	typedef Pimpact::ScalarField<Scalar,Ordinal> SF;
	typedef Pimpact::MultiField<VF> MVF;
	typedef Pimpact::MultiField<SF> MSF;
	typedef Pimpact::OperatorMV< Pimpact::Helmholtz<Scalar,Ordinal> >  Lap;
	typedef Pimpact::OperatorMV< Pimpact::Div_Hinv_Grad<Scalar,Ordinal> >  Schur;
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

	// flow type \todo + Boundary conditions
	int flow = 1;
	my_CLP.setOption( "flow", &flow, "Flow type: 0=zero flow, 1=2D Poiseuille flow in x, 2=2D Poiseuille flow in y" );

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
	*outPar << " \tre=" << re << "\n";
	*outPar << " \tomega=" << omega << "\n";

	auto ds = Pimpact::createDomainSize<Scalar>(l1,l2,l3);
	ds->set_Impact();
	ds->print( *outPar );

//	auto bc = Pimpact::createPeriodicChannelBC2D();
//	auto bc = Pimpact::createPeriodicChannelBC2D();
	auto bc = Pimpact::createAllDirichletBC2D();
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
	auto sca = Pimpact::createScalarField<double,int>(fS);
	auto vel = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto p     = Pimpact::createMultiField<SF,double,int>(*sca,1);
	auto temps = Pimpact::createMultiField<SF,double,int>(*sca,1);
	auto u     = Pimpact::createMultiField<VF,double,int>(*vel,1);
	auto f     = Pimpact::createMultiField<VF,double,int>(*vel,1);
	auto tempv = Pimpact::createMultiField<VF,double,int>(*vel,1);

	sca = Teuchos::null;
	vel = Teuchos::null;

	p->Init(0.);
	u->GetVec(0).init_field( Pimpact::EFlowProfile(flow) );
	u->Init(0);

	tempv->GetVec(0).init_field( Pimpact::EFlowProfile(flow) );
	f->GetVec(0).init_field( Pimpact::ZeroProf );

	// init operators
	auto lap  = Pimpact::createHelmholtz<Scalar,Ordinal>( 0., 1./re );
	auto div  = Pimpact::createOperatorMV<Pimpact::Div<double,int> >();
	auto grad = Pimpact::createOperatorMV<Pimpact::Grad<double,int> >();


	// solve parameter for GMRES
	RCP<ParameterList> solveParaGMRES = parameterList();

	solveParaGMRES->set ("Num Blocks", 100);
	solveParaGMRES->set ("Maximum Iterations", 2000);
	solveParaGMRES->set ("Maximum Restarts", 20);
	solveParaGMRES->set ("Convergence Tolerance", 1.0e-16);
	solveParaGMRES->set ("Output Frequency", 10 );
	solveParaGMRES->set ("Output Style", Belos::General );
	solveParaGMRES->set ("Output Stream", outLap1 );
	solveParaGMRES->set ("Verbosity",  Belos::Errors + Belos::Warnings + Belos::IterationDetails +
			Belos::OrthoDetails + Belos::FinalSummary +	Belos::TimingDetails +
			Belos::StatusTestDetails );

	Teuchos::writeParameterListToXmlFile( *solveParaGMRES, "para_solverGMRES.xml" );

	// solve parameter for CG
	RCP<ParameterList> solveParaCG = parameterList();

	solveParaCG->set ("Convergence Tolerance", 1.0e-16);
	solveParaCG->set ("Maximum Iterations", 1000);
	solveParaCG->set ("Output Frequency", 10 );
	solveParaCG->set ("Output Style", int(Belos::General) );
	solveParaCG->set ("Output Stream", outLap1 );
	solveParaCG->set ("Verbosity",  int(Belos::Errors + Belos::Warnings + Belos::IterationDetails +
			Belos::OrthoDetails + Belos::FinalSummary +	Belos::TimingDetails +
			Belos::StatusTestDetails ) );

	Teuchos::writeParameterListToXmlFile( *solveParaCG, "para_solverCG.xml" );

//	solveParaGMRES = Teuchos::getParametersFromXmlFile("solver1.xml");

	// create problems/solvers
	auto lap_problem = Pimpact::createLinearProblem<Scalar,MVF,Lap>( lap, u, f, solveParaCG, "CG" );

	auto schur = Pimpact::createDivHinvGrad<double,int>( u, lap_problem );

	solveParaGMRES->set( "Output Stream", outSchur );
	auto schur_prob = Pimpact::createLinearProblem<Scalar,MSF,Schur>( schur, p, temps, solveParaGMRES, "GMRES" );


	// solve stationary stokes
	lap_problem->solve( tempv, f );
	//	tempv->write( 2000 );

	div->apply( *tempv, *temps );
	//	temps->write( 2001 );

	solveParaCG->set( "Output Stream", outLap2 );
	solveParaCG->set ("Verbosity",  int(Belos::Errors) );
	lap_problem->setParameters( solveParaCG );
	schur_prob->solve( p, temps );
	p->write();

	grad->apply( *p, *tempv );
	//	tempv->write(2003);

	tempv->Add( -1., *tempv, 1., *f );

	solveParaCG->set ("Verbosity",  int(Belos::Errors + Belos::Warnings + Belos::IterationDetails +
			Belos::OrthoDetails + Belos::FinalSummary +	Belos::TimingDetails +
			Belos::StatusTestDetails ) );
	solveParaCG->set( "Output Stream", outLap2 );
	lap_problem->setParameters( solveParaCG );

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
