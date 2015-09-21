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

#include "Pimpact_MultiOpSmoother.hpp"
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


template<class T1,class T2> using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T> using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionHWOp<T> >;
template<class T> using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;


template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;


//template<class T> using PrecS = Pimpact::MultiOpSmoother< Pimpact::DivGradO2JSmoother<T> >;
template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2JSmoother >;
template<class T> using POP2  = Pimpact::PrecInverseOp< T, ConvDiffJT >;
template<class T> using POP3  = Pimpact::PrecInverseOp< T, ConvDiffSORT >;
//template<class T> using POP2 = Pimpact::PrecInverseOP< T, Pimpact::NonlinearSmoother<T,Pimpact::ConvectionDiffusionSORSmoother > >;
//template<class T> using POP2 = Pimpact::PrecInverseOP< T, ConvDiffJT >;







int main(int argi, char** argv ) {


/////////////////////////////////////////// set up parameters ///////////////////////////
  // intialize MPI
  MPI_Init( &argi, &argv );

	Teuchos::CommandLineProcessor my_CLP;

  std::string xmlFilename = "parameter3D.xml";
  my_CLP.setOption("filename", &xmlFilename, "file name of the input xml paramerterlist");

  my_CLP.recogniseAllOptions(true);
  my_CLP.throwExceptions(true);

  my_CLP.parse(argi,argv);


	auto pl = Teuchos::getParametersFromXmlFile( xmlFilename );
	pl->print();

/////////////////////////////////////////// end of set up parameters ///////////////////////////


///////////////////////////////////////////  set up initial stuff //////////////////////////////

	int initZero = pl->sublist("Solver").get<int>( "initZero", 0);

	int withprec=pl->sublist("Solver").get<int>( "withprec", 2 );
	int withoutput=pl->sublist("Solver").get<int>( "withoutput", 1 );

	int refinement     = pl->sublist("Solver").get<int>( "refinement level", 1     );
	S   refinementTol  = pl->sublist("Solver").get<S>(   "refinement tol",   1.e-6 );
	int refinementStep = pl->sublist("Solver").get<int>( "refinement step",  2     );

	auto space = Pimpact::createSpace<S,O,3,4>( Teuchos::rcpFromRef( pl->sublist("Space",true) ) );

	int baseflow = pl->get<int>("baseflow");
//	int flow = pl->get<int>("flow");
	int force = pl->get<int>("forcing");


	// init vectors
	auto x = Pimpact::create<MF>( space );

	// init Fields
	x->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->initField( Pimpact::EVectorField(baseflow), 1. );

	if( 0==initZero )
		x->init(0.);
	else if( 1==initZero ) {
		x->getFieldPtr(0)->getVFieldPtr()->random();
		x->getFieldPtr(0)->getSFieldPtr()->init(0.);
		x->scale(1.e-32);
	}
	else if( 2==initZero )
		x->random();
	else if( 3==initZero ) {
		x->getFieldPtr(0)->getSFieldPtr()->get0FieldPtr()->initField( Pimpact::Grad2D_inX, -2./space->getDomain()->getDomainSize()->getRe() );
	}



  /******************************************************************************************/
	for( int refine=0; refine<refinement; ++refine ) {

		std::string rl = "";
		if( refinement>1 )
			rl = std::to_string( (long long)refine ); // long long needed on brutus(intel)

		auto fu   = x->clone( Pimpact::ShallowCopy );
		if( 0==force )
			fu->init( 0. );
		else {
//			S re = space->getDomain()->getDomainSize()->getRe();
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

			// create Mulit space
			auto mgSpaces = Pimpact::createMGSpaces<FSpaceT,CSpaceT,CS>( space, pl->sublist("Multi Grid").get<int>("maxGrids") );



			Teuchos::RCP< Pimpact::OperatorBase<MVF> > opV2Vinv = Teuchos::null;

			if( withprec>1 ) {

				// creat H0-inv prec
				auto zeroOp = Pimpact::create<ConvDiffOpT>( space );

				pl->sublist("ConvDiff").sublist("Solver").set(
						"Output Stream",
						Pimpact::createOstream( zeroOp->getLabel()+rl+".txt", space->rankST() ));

				auto zeroInv =
					Pimpact::create<MOP>( zeroOp, Teuchos::rcpFromRef( pl->sublist("ConvDiff") ) );

				auto modeOp =
					Teuchos::rcp( new Pimpact::ModeNonlinearOp< ConvDiffOpT<SpaceT> >( zeroOp ) );

				pl->sublist("M_ConvDiff").sublist("Solver").set(
						"Output Stream",
						Pimpact::createOstream( modeOp->getLabel()+rl+".txt", space->rankST() ) );

				auto modeInv =
					Pimpact::create<MOP>(
							modeOp,
							Teuchos::rcpFromRef( pl->sublist("M_ConvDiff") ) );


				auto mgConvDiff =
					Pimpact::createMultiGrid<
					Pimpact::VectorField,
					TransVF,
					RestrVF,
					InterVF,
					ConvDiffOpT,
					ConvDiffOpT,
//					ConvDiffSORT,
					ConvDiffJT,
					//					MOP
					POP2
						> ( mgSpaces, Teuchos::rcpFromRef( pl->sublist("ConvDiff").sublist("Multi Grid") ) ) ;

				if( 0==space->rankST() )
					mgConvDiff->print();

				if( withprec>2 ) {
					zeroInv->getOperatorPtr()->getLinearProblem()->setRightPrec( Pimpact::createMultiOperatorBase(mgConvDiff) );
					modeInv->getOperatorPtr()->getLinearProblem()->setRightPrec( Pimpact::createMultiOperatorBase( Pimpact::create<Pimpact::EddyPrec>( zeroInv ) ) );
				}

				// create Hinv prec
				Teuchos::RCP<Pimpact::OperatorBase<MVF> >
					opV2Vprec = 
					Pimpact::createMultiOperatorBase(
							Pimpact::createMultiHarmonicDiagOp(
								zeroInv, 
								modeInv
								)
							);


				if( withprec>3 )
					opV2Vinv = opV2Vprec;
				else {
					pl->sublist("MH_ConvDiff").sublist("Solver").set(
							"Output Stream",
							Pimpact::createOstream( opV2V->getLabel()+rl+".txt", space->rankST() ) );

					auto opV2Vprob =
						Pimpact::createLinearProblem<MVF>(
								Pimpact::createMultiOperatorBase( opV2V ),
								Teuchos::null,
								Teuchos::null,
								Teuchos::rcpFromRef( pl->sublist("MH_ConvDiff").sublist("Solver") ),
								pl->sublist("MH_ConvDiff").get<std::string>("Solver name","GMRES") );

					opV2Vprob->setRightPrec( opV2Vprec );
					opV2Vinv = Pimpact::createInverseOperatorBase( opV2Vprob );
				}

			}


			/////////////////////////////////////////end of opv2v///////////////////////////////////////////////////////////////////

			////--- inverse DivGrad
			pl->sublist("DivGrad").sublist("Solver").set(
					"Output Stream",
					Pimpact::createOstream( "DivGrad"+rl+".txt", space->rankST() ) );

			auto divGradInv2 =
						Pimpact::createInverseOp( 
							Pimpact::createDivGradOp(
								opV2S->getOperatorPtr(),
								opS2V->getOperatorPtr()
								)
							, Teuchos::rcpFromRef( pl->sublist("DivGrad") )
							);

			{ // init multigrid divgrad

				auto mgDivGrad = 
					Pimpact::createMultiGrid<
					Pimpact::ScalarField,
					Pimpact::TransferOp,
					Pimpact::RestrictionHWOp,
					Pimpact::InterpolationOp,
					Pimpact::DivGradOp,
					Pimpact::DivGradO2Op,
					Pimpact::DivGradO2JSmoother,
					POP
					>( mgSpaces, Teuchos::rcpFromRef( pl->sublist("DivGrad").sublist("Multi Grid") ) );

				if( 0==space->rankST() )
					mgDivGrad->print();

				divGradInv2->setRightPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );
			}

			auto divGradInv =
				Pimpact::createMultiHarmonicMultiOpWrap(
						divGradInv2
						);

			// schurcomplement approximator ...
			auto opSchur =
				Pimpact::createTripleCompositionOp(
						opV2S,
						opV2V,
						opS2V
						);

			auto opS2Sinv =
				Pimpact::createTripleCompositionOp(
						divGradInv,
						opSchur,
						divGradInv
						);

			if( space->rankST()==0 )
				std::cout << opS2Sinv->getLabel() << "\n";

			auto invTriangOp =
				Pimpact::createInverseTriangularOp(
						Pimpact::createMultiOpUnWrap( opV2Vinv ),
						opS2V,
						opS2Sinv );

			lprec = Pimpact::createMultiOperatorBase( invTriangOp );
			
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
//		if( withprec ) lp_->setLeftPrec( lprec );

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


    // Get the answer
    *group = solver->getSolutionGroup();

	x = Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr();
	if( withoutput ) {
		Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr()->write(refine*1000 );
		Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr()->getFieldPtr(0)->getVFieldPtr()->write( 500+refine*1000, true );
		//		Teuchos::rcp_dynamic_cast<const NV>( group->getFPtr() )->getConstFieldPtr()->write(900);
	}

	{
		auto out = Pimpact::createOstream( "energy_dis"+rl+".txt", space->rankST() );
		*out << 0 << "\t" << x->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->norm() << "\t" << std::sqrt(x->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->getLength()) << "\n";
		for( int i=0; i<space->nGlo(3); ++i )
			*out << i+1 << "\t" << x->getFieldPtr(0)->getVFieldPtr()->getFieldPtr(i)->norm() << "\t" << std::sqrt( (S)x->getFieldPtr(0)->getVFieldPtr()->getFieldPtr(i)->getLength()) << "\n";
	}

	// spectral refinement of x, fu
	//std::cout<< "bla bool: "<<       (             x->getFieldPtr(0)->getVFieldPtr()->getFieldPtr(space->nGlo(3)-1)->norm() /
			//std::sqrt( (S)x->getFieldPtr(0)->getVFieldPtr()->getFieldPtr(space->nGlo(3)-1)->getLength() ) <	refinementTol) << "\n";
	//std::cout<< "bla e: "<<                    x->getFieldPtr(0)->getVFieldPtr()->getFieldPtr(space->nGlo(3)-1)->norm()  << "\n";
	//std::cout<< "bla sqrt: "<< std::sqrt( (S)x->getFieldPtr(0)->getVFieldPtr()->getFieldPtr(space->nGlo(3)-1)->getLength() ) << "\n";

	if(                   x->getFieldPtr(0)->getVFieldPtr()->getFieldPtr(space->nGlo(3)-1)->norm() /
			std::sqrt( (S)x->getFieldPtr(0)->getVFieldPtr()->getFieldPtr(space->nGlo(3)-1)->getLength() ) <	refinementTol ) break;

	if( refinement>1 ) {
		auto spaceF =
			Pimpact::RefinementStrategy<SpaceT>::createRefinedSpace(
					space, Teuchos::tuple<int>( 0, 0, 0, refinementStep ) );
		//				space, Teuchos::tuple<int>( 1, 1, 1, refinementStep ) );

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

		//			if( withoutput )  xf->write();

		x = Pimpact::createMultiField( xf );
		space = spaceF;
	}
	} 
	/******************************************************************************************/

	Teuchos::TimeMonitor::summarize();

	Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );

  MPI_Finalize();
  return( 0 );

}
