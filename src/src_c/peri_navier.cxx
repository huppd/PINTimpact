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


//template<class ST> using BOPF = Pimpact::MultiOpWrap< Pimpact::DivGradOp<ST> >;
//template<class ST> using BOPC = Pimpact::MultiOpWrap< Pimpact::DivGradO2Op<ST> >;
//template<class ST> using BSM = Pimpact::MultiOpWrap< Pimpact::DivGradO2JSmoother<ST> >;

template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;


int main(int argi, char** argv ) {


  // intialize MPI
  MPI_Init( &argi, &argv );

  //// get problem values form Comand line
  Teuchos::CommandLineProcessor my_CLP;

  // physical constants
  S re = 1.e2;

  S alpha2 = 1.e3;

  // flow type
  int flow = 7;

  // domain type
  int domain = 2;

  // domain size
  int dim = 2;
  S l1 = 1.;
  S l2 = 1.;
  S l3 = 1.;

  // grid size
  O n1 = 65;
  O n2 = 65;
  O n3 = 17;
  //O n3 = 2.;

  O nf = 4;

  O nfs = 4;

  O nfe = 4;

  // processor grid size
  O np1 = 1;
  O np2 = 1;
  O np3 = 1.;

  // solver stuff
  std::string linSolName = "Block GMRES";
  std::string nonLinSolName = "Newton";

  std::string lineSearchName = "Backtrack";

	bool withprec=true;

  int maxIter = 10;

  S tolBelos = 2.e-1;
  S tolNOX   = 1.e-1;
  S tolNF    = 1.e-1;

  // end of parsing

  if( nfs==nfe ) {
    nfs=nf;
    nfe=nf+1;
  }
  else {
    nf=nfe;
  }
  // starting with ininializing
  auto pl = Teuchos::parameterList();

  pl->set( "Re", re );
  pl->set( "alpha2", alpha2 );
  pl->set( "domain", domain );

  pl->set( "lx", l1 );
  pl->set( "ly", l2 );
  pl->set( "lz", l3 );

  pl->set( "dim", dim );

  pl->set("nx", n1 );
  pl->set("ny", n2 );
  pl->set("nz", n3 );

  pl->set("nf", nf );
  pl->set("nfs", nfs );
  pl->set("nfe", nfe );

  // processor grid size
  pl->set("npx", np1 );
  pl->set("npy", np2 );
  pl->set("npz", np3 );

  auto space = Pimpact::createSpace<S,O,3,4>( pl );

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

  *outPar << " \tflow=" << flow << "\n";
  *outPar << " \tdomain=" << domain << "\n";



  if( space->rankST()==0 ) {
    Teuchos::rcp_static_cast<std::ofstream>(outPar)->close();
  }
  outPar = Teuchos::null;

  // init vectors
  auto x    = Pimpact::createMultiField(
			Pimpact::createCompoundField(
      	Pimpact::createMultiHarmonicVectorField( space, nfs ),
      	Pimpact::createMultiHarmonicScalarField( space, nfs )) );
  auto fu   = x->clone();

  // init Fields, init and rhs
  x->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->initField( Pimpact::EFlowField(flow), 1 );

	x->init(0.);
  fu->init( 0. );


	//tolBelos*=l1*l2/n1/n2*(nfe-1)/nfs;

  /******************************************************************************************/
  for( nf=nfs; nf<nfe; nf*=2) {

	tolBelos/=2;

    if( nf!=nfs ) {
      S toltemp = x->getConstFieldPtr(0)->getConstVFieldPtr()->getConstFieldPtr(x->getConstFieldPtr(0)->getConstVFieldPtr()->getNumberModes()-1)->norm()/std::sqrt(l1*l2/n1/n2);
      if( 0==space->rankST() ) std::cout << "\n\t--- ||u_Nf||: "<<toltemp<<"\t---\n";
      if( toltemp < tolNF ) {
        if( 0==space->rankST() ) std::cout << "\n\t--- Nf: "<<x->getConstFieldPtr(0)->getConstVFieldPtr()->getNumberModes()<<"\tdof: "<<x->getLength(true)<<"\t---\n";
        break;
      }
      do {
        x->getFieldPtr(0)->getVFieldPtr()->push_back();
        x->getFieldPtr(0)->getSFieldPtr()->push_back();
        fu->getFieldPtr(0)->getVFieldPtr()->push_back();
        fu->getFieldPtr(0)->getSFieldPtr()->push_back();
      } while( x->getConstFieldPtr(0)->getConstVFieldPtr()->getNumberModes() < nf );
      tolNOX /= 10;
    }

    auto para = Pimpact::createLinSolverParameter( linSolName, tolBelos, -1 );
		para->set( "Maximum Iterations", 1000 );
		para->set( "Implicit Residual Scaling", "Norm of RHS" );
		para->set( "Explicit Residual Scaling", "Norm of RHS" );
		para->set( "Flexible Gmres", true );
		para->set( "Output Frequency", 1 );
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
				Pimpact::createMultiDtConvectionDiffusionOp( space, nf );
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

			auto mgSpaces = Pimpact::createMGSpaces<FSpaceT,CSpaceT,CS>( space, 2 );

			auto opV2Vprob =
					Pimpact::createLinearProblem<MVF>(
							Pimpact::createMultiOperatorBase(
									opV2V ),
									Teuchos::null,
									Teuchos::null,
									Pimpact::createLinSolverParameter( "GMRES", 1.e-6, -1, outPrec ),
									"GMRES" );

			auto opV2Vinv = Pimpact::createInverseOperatorBase( opV2Vprob );

			// opS2Sinv alias schurcomplement ...
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
							Pimpact::createLinSolverParameter( "Block GMRES", 1.e-9, -1, outSchur ),
							"Block GMRES" );

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

			divGradProb->setRightPrec( mg_divGrad );

			auto divGradInv = Pimpact::createInverseOperatorBase( divGradProb );

			auto opS2Sinv =
					Pimpact::createOperatorBase(
							Pimpact::createTripleCompositionOp(
									divGradInv,
									opSchur,
									divGradInv
							)
					);

			auto invSchur = Pimpact::createInverseTriangularOp(
					opV2Vinv,
					opS2V,
					opS2Sinv );

			lprec =
				Pimpact::createMultiOperatorBase(
						invSchur );
		}


    auto lp_ = Pimpact::createLinearProblem<MF>( jop, x->clone(), fu->clone(), para, linSolName );
		if( withprec ) lp_->setRightPrec( lprec );

    auto lp = Pimpact::createInverseOperatorBase( lp_ );

    auto inter = NOX::Pimpact::createInterface( fu, op, lp );

    auto nx = NOX::Pimpact::createVector(x);

    auto bla = Teuchos::parameterList();

    auto group = NOX::Pimpact::createGroup<Inter>( bla, inter, nx );

    // Set up the status tests
    auto statusTest = NOX::Pimpact::createStatusTest( maxIter, tolNOX, tolBelos );

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
