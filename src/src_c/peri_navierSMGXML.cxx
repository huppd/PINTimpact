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
#include "Pimpact_CompoundSmoother.hpp"

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

template<class OpT>
using Smoother =
	Pimpact::CompoundSmoother<
		OpT,
		MOP,
//		Pimpact::MultiHarmonicOpWrap< Pimpact::DivGradO2JSmoother<Pimpact::DivGradO2Op<CSpaceT> > >//
		Pimpact::MultiHarmonicOpWrap< MOP<Pimpact::DivGradO2Op<CSpaceT> > >//
	>;


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

  std::string xmlFilename = "parameterSMG.xml";
  my_CLP.setOption("filename", &xmlFilename, "file name of the input xml paramerterlist");

  my_CLP.recogniseAllOptions(true);
  my_CLP.throwExceptions(true);

  my_CLP.parse(argi,argv);


	auto pl = Teuchos::getParametersFromXmlFile( xmlFilename );
	pl->print();

	/////////////////////////////////////////// end of set up parameters ///////////////////////////

	///////////////////////////////////////////  set up initial stuff //////////////////////////////
	
	int initZero = pl->sublist("Solver").get<int>("initZero");
	Pimpact::EVectorField baseflow = (Pimpact::EVectorField)pl->get<int>( "baseflow", 1 );
	
	auto space = Pimpact::createSpace<S,O,3,4>( Teuchos::rcpFromRef( pl->sublist("Space",true) ) );

	// init vectors
	auto x = Pimpact::create<MF>( space );

	// init Fields
	x->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->initField( baseflow, 1. );

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

  /******************************************************************************************/
	int refinement=pl->sublist("Solver").get<int>("refinement");
	int withprec=pl->sublist("Solver").get<int>("withprec");

	for( int refine=0; refine<refinement; ++refine ) {

		std::string rl = "";
		if( refinement>1 )
			rl = std::to_string( (long long)refine ); // long long needed on brutus(intel)

		auto fu   = x->clone( Pimpact::ShallowCopy );
		if( 0==pl->get<int>("forcing", 1) )
			fu->init( 0. );
		else {
			//			S re = space->getDomainSize()->getRe();
			fu->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->getFieldPtr(Pimpact::U)->initField(  Pimpact::FPoint, pl->get<S>( "lambda0x", -2. ) );
			fu->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->getFieldPtr(Pimpact::V)->initField( Pimpact::FPoint, pl->get<S>( "lambdaCy",  1. ) );
			fu->getFieldPtr(0)->getVFieldPtr()->getCFieldPtr(0)->getFieldPtr(Pimpact::W)->initField( Pimpact::FPoint, pl->get<S>( "lambdaCz",  0. ) );
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

      auto pls = Teuchos::rcpFromRef( pl->sublist("Multi Grid") );

			auto mg =
				Pimpact::createMultiGrid<
					FT,
					TransF,
					RestrF,
					InterF,
					DTConvDiffOpT,
					DTConvDiffOpT,
//					MOP,
					Smoother,
					MOP > ( mgSpaces, pls ) ;

			mg->print();

			lprec = Pimpact::createMultiOperatorBase( mg );

		}

		/*** end of init preconditioner ****************************************************************************/



		pl->sublist("Belos Solver").sublist("Solver").set( "Output Stream", Pimpact::createOstream("stats_linSolve"+rl+".txt",space->rankST()) );

    auto lp_ = Pimpact::createLinearProblem<MF>(
				jop,
				x->clone(),
				fu->clone(),
				Teuchos::rcpFromRef(pl->sublist("Belos Solver").sublist("Solver") ),
				pl->sublist("Belos Solver").get<std::string>("Solver name") );

		if( withprec ) lp_->setRightPrec( lprec );

    auto lp = Pimpact::createInverseOperatorBase( lp_ );

    auto inter = NOX::Pimpact::createInterface( fu, op, lp );

    auto nx = NOX::Pimpact::createVector(x);

    auto bla = Teuchos::parameterList();

    auto group = NOX::Pimpact::createGroup<Inter>( bla, inter, nx );

    // Set up the status tests
		auto statusTest =
			NOX::StatusTest::buildStatusTests( pl->sublist("NOX Solver").sublist("Status Test"), NOX::Utils() );


    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver =
        NOX::Solver::buildSolver( group, statusTest, Teuchos::rcpFromRef( pl->sublist("NOX Solver") ) );

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

		if( pl->sublist("Solver").get<int>("withoutput") ) {
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
			xf->getVFieldPtr()->get0FieldPtr()->initField( baseflow, 1. );
			auto temp = Pimpact::create<CF>( spaceF );

			refineOp->apply( x->getField(0), *temp );

			xf->add( 1., *temp, 0., *temp );

			x = Pimpact::createMultiField( xf );
			space = spaceF;
		}
  } 
	/******************************************************************************************/

	Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );

  MPI_Finalize();
  return( 0 );

}
