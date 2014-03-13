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
#include "BelosPimpactAdapter.hpp"




int main(int argi, char** argv ) {

	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp; // Save some typing

	typedef double													Scalar;
	typedef int															Ordinal;
	typedef double					                S;
	typedef int   						              O;
	typedef Pimpact::ScalarField<S,O>       SF;
	typedef Pimpact::VectorField<S,O>       VF;
	typedef Pimpact::ModeField<SF>          MSF;
	typedef Pimpact::ModeField<VF>          MVF;
	typedef Pimpact::CompoundField<MVF,MSF> CF;
	typedef Pimpact::MultiField<CF>         MV;
	typedef Pimpact::OperatorMV<Pimpact::CompoundStokes<S,O> > OP;
//	typedef Pimpact::OperatorMV< Pimpact::DtL<Scalar,Ordinal> >  Lap;
//	typedef Pimpact::OperatorMV< Pimpact::Div_DtLinv_Grad<Scalar,Ordinal> >  Schur;
//	typedef Pimpact::OperatorMV< Pimpact::Grad<Scalar,Ordinal> >  G;

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
	my_CLP.setOption( "px", &px, "pressure gradient" );

	// flow type \todo + Boundary conditions
	int flow = 1;
	my_CLP.setOption( "flow", &flow,
			"Flow type: 0=zero flow, 1=2D Poiseuille flow in x, 2=2D Poiseuille flow in y, 3=2D pulsatile flow in x, 4=2D pulsatile flow in y" );

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
	*outPar << " \tpx=" << px << "\n";
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
	auto scac = Pimpact::createScalarField<double,int>(fS);
	auto scas = Pimpact::createScalarField<double,int>(fS);

	auto velc = Pimpact::createVectorField<double,int>(fS,iIS,fIS);
	auto vels = Pimpact::createVectorField<double,int>(fS,iIS,fIS);

	auto vel = Pimpact::createModeField( velc, vels );
	auto sca = Pimpact::createModeField( scac, scas );

	auto q_ 	 = Pimpact::createCompoundField( vel, sca );

	auto q     = Pimpact::createMultiField< CF >(*q_,1);
	auto f     = Pimpact::createMultiField< CF >(*q_,1);

	scac=Teuchos::null;
	scas=Teuchos::null;
	sca=Teuchos::null;
	vel=Teuchos::null;
//	velc=Teuchos::null;
	vels=Teuchos::null;


	switch( Pimpact::EFlowType(flow) ) {
	case Pimpact::Zero2DFlow :
		q->getField(0).getVField()->getCFieldPtr()->initField( Pimpact::ZeroProf );
		q->getField(0).getVField()->getSFieldPtr()->initField( Pimpact::ZeroProf );
		break;
	case Pimpact::Poiseuille_inX :
		q->getField(0).getVField()->getCFieldPtr()->initField( Pimpact::Poiseuille2D_inX );
		q->getField(0).getVField()->getSFieldPtr()->initField( Pimpact::ZeroProf );
		break;
	case Pimpact::Poiseuille_inY :
		q->getField(0).getVField()->getCFieldPtr()->initField( Pimpact::Poiseuille2D_inY );
		q->getField(0).getVField()->getSFieldPtr()->initField( Pimpact::ZeroProf );
		break;
	case Pimpact::Pulsatile_inX :
		q->getField(0).getVField()->getCFieldPtr()->initField( Pimpact::Pulsatile2D_inXC, re, omega, px );
		q->getField(0).getVField()->getSFieldPtr()->initField( Pimpact::Pulsatile2D_inXS, re, omega, px );
		break;
	case Pimpact::Pulsatile_inY :
		q->getField(0).getVField()->getCFieldPtr()->initField( Pimpact::Pulsatile2D_inYC, re, omega, px );
		q->getField(0).getVField()->getSFieldPtr()->initField( Pimpact::Pulsatile2D_inYS, re, omega, px );
		break;
	}

	q->init(0);

	f->getField(0).getVField()->getCFieldPtr()->initField( Pimpact::ZeroProf );
	f->getField(0).getVField()->getSFieldPtr()->initField( Pimpact::ZeroProf );

	// init operators
	auto op = Pimpact::createOperatorMV(
	    Pimpact::createCompoundStokes( omega, 0., 1./re, velc ) );


	RCP<ParameterList> solverParams = parameterList();

	solverParams->set ("Num Blocks", 1000);
	solverParams->set ("Maximum Iterations", 20000);
	solverParams->set ("Maximum Restarts", 2000);
	solverParams->set ("Convergence Tolerance", 1.0e-16*l1*l2/n1/n2);
	solverParams->set ("Implicit Residual Scaling","None");
	solverParams->set ("Explicit Residual Scaling","None");
	solverParams->set ("Output Frequency", 10 );
	solverParams->set ("Output Style", 1 );
	solverParams->set ("Output Stream", outLap1 );
	solverParams->set ("Verbosity",  int(Belos::Errors + Belos::Warnings + Belos::IterationDetails +
			Belos::OrthoDetails + Belos::FinalSummary +	Belos::TimingDetails +
			Belos::StatusTestDetails + Belos::Debug) );

	Teuchos::writeParameterListToXmlFile( *solverParams, "para_solver.xml" );

//	solverParams = Teuchos::getParametersFromXmlFile("solver1.xml");

	// create problems/solvers
	auto problem = Pimpact::createLinearProblem<Scalar,MV,OP>( op, q, f, solverParams, "GMRES" );

	// solve stationary stokes
	problem->solve( q, f );

	q->write();


	if(rank==0) {
		Teuchos::rcp_static_cast<std::ofstream>( outLap1)->close();
		Teuchos::rcp_static_cast<std::ofstream>( outLap2)->close();
		Teuchos::rcp_static_cast<std::ofstream>(outSchur)->close();
	}


	MPI_Finalize();
	return( 0 );
}
