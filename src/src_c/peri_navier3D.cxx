#include <fstream>
#include <ostream>

#include <mpi.h>

#include <Teuchos_Array.hpp>
#include "Teuchos_CommandLineProcessor.hpp"
#include <Teuchos_CommHelpers.hpp>
#include "Teuchos_Range1D.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Tuple.hpp>
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

#include "BelosOutputManager.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosTypes.hpp"

#include "NOX.H"

#include "BelosPimpactAdapter.hpp"
#include "NOX_Pimpact_Group.hpp"
#include "NOX_Pimpact_Interface.hpp"
#include "NOX_Pimpact_Vector.hpp"
#include "NOX_Pimpact_StatusTest.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"
#include "Pimpact_CoarsenStrategy.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_FieldFactory.hpp"
#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_ModeNonlinearOp.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_MultiOpSmoother.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorFactory.hpp"
#include "Pimpact_RefinementStrategy.hpp"
#include "Pimpact_TransferCompoundOp.hpp"
#include "Pimpact_TransferMultiHarmonicOp.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_Space.hpp"





using S = double;
using O = int;

const int dNC = 4;

using SpaceT = Pimpact::Space<S,O,4,dNC>;

using FSpaceT = Pimpact::Space<S,O,4,dNC>;
using CSpaceT = Pimpact::Space<S,O,4,2>;

using CS = Pimpact::CoarsenStrategyGlobal<FSpaceT,CSpaceT>;

using VF = Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> >;
using SF = Pimpact::MultiHarmonicField< Pimpact::ScalarField<SpaceT> >;

using MVF = Pimpact::MultiField<VF>;
using MSF = Pimpact::MultiField<SF>;

using CF = Pimpact::CompoundField< VF, SF>;
using MF = Pimpact::MultiField<CF>;

using BOp = Pimpact::OperatorBase<MF>;

using NV = NOX::Pimpact::Vector<MF>;


template<class T1,class T2>
using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T>
using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionHWOp<T> >;
template<class T>
using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;


//template<class T> using PrecS = Pimpact::MultiOpSmoother< Pimpact::DivGradO2JSmoother<T> >;
//template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2SORSmoother >;
template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2JSmoother >;
//template<class T> using POP   = Pimpact::PrecInverseOp< T, Pimpact::Chebyshev >;
template<class T> using POP2  = Pimpact::PrecInverseOp< T, ConvDiffJT >;
template<class T> using POP3  = Pimpact::PrecInverseOp< T, ConvDiffSORT >;
//template<class T> using POP2 = Pimpact::PrecInverseOP< T, Pimpact::NonlinearSmoother<T,Pimpact::ConvectionDiffusionSORSmoother > >;
//template<class T> using POP2 = Pimpact::PrecInverseOP< T, ConvDiffJT >;







int main( int argi, char** argv ) {


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

	/////////////////////////////////////////// end of set up parameters /////////////////////////


	///////////////////////////////////////////  set up initial stuff ////////////////////////////

	std::string initGuess = pl->sublist("Solver").get<std::string>( "initial guess", "zero" );

	int withoutput=pl->sublist("Solver").get<int>( "withoutput", 1 );

	int refinement     = pl->sublist("Solver").get<int>( "refinement level", 1     );
	S   refinementTol  = pl->sublist("Solver").get<S>(   "refinement tol",   1.e-6 );
	int refinementStep = pl->sublist("Solver").get<int>( "refinement step",  2     );

	Teuchos::RCP<const Pimpact::Space<S,O,4,dNC> > space =
		Pimpact::createSpace<S,O,4,dNC>(
				Teuchos::rcpFromRef( pl->sublist( "Space", true ) ) );


	// init vectors
	Teuchos::RCP<MF> x = createMultiField(
			createCompoundField(
				Teuchos::rcp( new VF(space,true) ),
				Teuchos::rcp( new SF(space) ) ) ) ;

	// init Fields
	x->getFieldPtr(0)->getVFieldPtr()->initField( pl->sublist("Base flow") );


	if( "zero"==initGuess )
		x->init( 0. );
	else if( "almost zero"==initGuess ) {
		x->random();
		x->scale(1.e-12);
	}
	else if( "random"==initGuess )
		x->random();
	else if( "exact"==initGuess ) {
		x->getFieldPtr(0)->getSFieldPtr()->get0FieldPtr()->initField( Pimpact::Grad2D_inX, -2./space->getDomainSize()->getRe() );
	}
	x->getFieldPtr(0)->getVFieldPtr()->changed();
	x->getFieldPtr(0)->getSFieldPtr()->changed();

	/******************************************************************************************/
	for( int refine=0; refine<refinement; ++refine ) {

		std::string rl = "";
		if( refinement>1 )
			rl = std::to_string( static_cast<long long>(refine) ); // long long needed on brutus(intel)

		auto fu = x->clone( Pimpact::ShallowCopy );
		fu->getFieldPtr(0)->getVFieldPtr()->initField( pl->sublist("Force") );

		//if( withoutput )
			//fu->write( 90000 );


		auto opV2V = Pimpact::createMultiDtConvectionDiffusionOp( space );
		auto opS2V = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::GradOp>( space ) );
		auto opV2S = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::DivOp>( space ) );

		auto op = Pimpact::createCompoundOpWrap(
					opV2V,
					opS2V,
					opV2S );

		pl->sublist("Picard Solver").sublist("Solver").set( "Output Stream", Pimpact::createOstream("Picard"+rl+".txt", space->rankST() ) );

		auto opInv = Pimpact::createInverseOp( op, Teuchos::rcpFromRef( pl->sublist("Picard Solver") ) );

		/*** init preconditioner ****************************************************************************/

		std::string picardPrecString = pl->sublist("Picard Solver").get<std::string>( "preconditioner", "none" );

		if( "none" != picardPrecString ) { 

			// create Multi space
			auto mgSpaces = Pimpact::createMGSpaces<FSpaceT,CSpaceT,CS>( space, pl->sublist("Multi Grid").get<int>("maxGrids") );
			//auto mgSpaces = Pimpact::createMGSpaces<FSpaceT,FSpaceT,CS>( space, pl->sublist("Multi Grid").get<int>("maxGrids") );

			// creat H0-inv prec
			auto zeroOp = Pimpact::create<ConvDiffOpT>( space );

			pl->sublist("ConvDiff").sublist("Solver").set( "Output Stream",
					Pimpact::createOstream( zeroOp->getLabel()+rl+".txt", space->rankST() ) );

			auto zeroInv = Pimpact::createMultiOpUnWrap( Pimpact::createInverseOp(
						zeroOp, Teuchos::rcpFromRef( pl->sublist("ConvDiff") ) ) );

			auto modeOp = Teuchos::rcp( new
					Pimpact::ModeNonlinearOp<ConvDiffOpT<SpaceT> >( zeroOp ) );

			pl->sublist("M_ConvDiff").sublist("Solver").set(
					"Output Stream",
					Pimpact::createOstream( modeOp->getLabel()+rl+".txt", space->rankST() ) );

			auto modeInv = Pimpact::createMultiOpUnWrap( Pimpact::createInverseOp(
						modeOp, Teuchos::rcpFromRef( pl->sublist("M_ConvDiff") ) ) );

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
				POP2
					> ( mgSpaces, Teuchos::rcpFromRef( pl->sublist("ConvDiff").sublist("Multi Grid") ) ) ;

			//if( 0==space->rankST() )
				//mgConvDiff->print();

			zeroInv->getOperatorPtr()->setRightPrec( Pimpact::createMultiOperatorBase(mgConvDiff) );
			modeInv->getOperatorPtr()->setRightPrec( Pimpact::createMultiOperatorBase( Pimpact::create<Pimpact::EddyPrec>(zeroInv) ) );

			// create Hinv prec
			Teuchos::RCP<Pimpact::OperatorBase<MVF> > opV2Vprec = 
				Pimpact::createMultiOperatorBase(
						Pimpact::createMultiHarmonicDiagOp(
							zeroInv, 
							modeInv ) );


			pl->sublist("MH_ConvDiff").sublist("Solver").set(
					"Output Stream",
					Pimpact::createOstream( opV2V->getLabel()+rl+".txt", space->rankST() ) );


			auto opV2Vinv = Pimpact::createInverseOp( opV2V, Teuchos::rcpFromRef(
						pl->sublist("MH_ConvDiff") ) );

			opV2Vinv->setRightPrec( opV2Vprec);

			/////////////////////////////////////////end of opv2v//////////////////////////////////////
			////--- inverse DivGrad

			pl->sublist("DivGrad").sublist("Solver").set( "Output Stream",
					Pimpact::createOstream( "DivGrad"+rl+".txt", space->rankST() ) );

			auto divGradOp =
				Pimpact::createDivGradOp(
						opV2S->getOperatorPtr(),
						opS2V->getOperatorPtr() );

			auto divGradInv2 =
				Pimpact::createInverseOp( 
						divGradOp,
						Teuchos::rcpFromRef( pl->sublist("DivGrad") ) );

			std::string divGradScalString =
				pl->sublist("DivGrad").get<std::string>("scaling","none");

			if( "none" != divGradScalString ) { 
				if( "left" != divGradScalString )
					divGradInv2->setLeftPrec( Pimpact::createMultiOperatorBase(
								Pimpact::createInvDiagonal( divGradOp ) ) );
				if( "right" != divGradScalString )
					divGradInv2->setRightPrec( Pimpact::createMultiOperatorBase(
								Pimpact::createInvDiagonal( divGradOp ) ) );
			}

			std::string divGradPrecString =
				pl->sublist("DivGrad").get<std::string>("preconditioner","none");

			if( "none" != divGradPrecString ) { // init multigrid divgrad

				auto mgDivGrad = Pimpact::createMultiGrid<
					Pimpact::ScalarField,
					Pimpact::TransferOp,
					Pimpact::RestrictionHWOp,
					Pimpact::InterpolationOp,
					Pimpact::DivGradOp,
					//Pimpact::DivGradOp,
					Pimpact::DivGradO2Op,
					Pimpact::DivGradO2JSmoother,
					//Pimpact::Chebyshev,
					//Pimpact::DivGradO2SORSmoother,
					//POP
					//Pimpact::Chebyshev
					//Pimpact::DivGradO2Inv
					//Pimpact::DivGradO2SORSmoother
					Pimpact::DivGradO2JSmoother
						>( mgSpaces, Teuchos::rcpFromRef( pl->sublist("DivGrad").sublist("Multi Grid") ) );

				if( 0==space->rankST() )
					mgDivGrad->print();

				if( "right" == divGradPrecString ) 
					divGradInv2->setRightPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );
				if( "left" == divGradPrecString ) 
					divGradInv2->setLeftPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );
			}

			auto divGradInv =
				Pimpact::createMultiHarmonicMultiOpWrap(
						divGradInv2 );

			// schurcomplement approximator ...
			auto opSchur =
				Pimpact::createTripleCompositionOp(
						opV2S,
						opV2V,
						opS2V );

			auto opS2Sinv =
				Pimpact::createTripleCompositionOp(
						divGradInv,
						opSchur,
						divGradInv );

			if( space->rankST()==0 )
				std::cout << opS2Sinv->getLabel() << "\n";

			auto invTriangOp =
				Pimpact::createInverseTriangularOp(
						Pimpact::createMultiOpUnWrap( opV2Vinv ),
						opS2V,
						opS2Sinv );

			if( "right" == picardPrecString ) 
				opInv->setRightPrec( Pimpact::createMultiOperatorBase( invTriangOp ) );
			if( "left" == picardPrecString ) 
				opInv->setLeftPrec( Pimpact::createMultiOperatorBase( invTriangOp ) );
		}
		/*** end of init preconditioner ****************************************************************************/

		auto inter = NOX::Pimpact::createInterface( fu, createMultiOpWrap(op), opInv );

		auto nx = NOX::Pimpact::createVector( x );

		auto bla = Teuchos::parameterList();

		auto group = NOX::Pimpact::createGroup( bla, inter, nx );

		// Set up the status tests
		auto statusTest =
			NOX::StatusTest::buildStatusTests( pl->sublist("NOX Solver").sublist("Status Test"), NOX::Utils() );

		// Create the solver
		Teuchos::RCP<NOX::Solver::Generic> solver =
			NOX::Solver::buildSolver( group, statusTest, Teuchos::rcpFromRef( pl->sublist("NOX Solver") ) );


		Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );
		if( 0==space->rankST() )
			std::cout << "\n\t--- Nf: "<< space->nGlo(3) <<"\tdof: "<<x->getLength()<<"\t---\n";

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
			Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr()->level();
			Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr()->write( refine*1000 );
			Teuchos::rcp_const_cast<NV>(Teuchos::rcp_dynamic_cast<const NV>( group->getXPtr() ))->getFieldPtr()->getFieldPtr(0)->getVFieldPtr()->write( 500+refine*1000, true );
		}

		{
			auto out = Pimpact::createOstream( "energy_dis"+rl+".txt", space->rankST() );

			*out << 0 << "\t" << x->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->norm() << "\t" << std::sqrt( static_cast<S>( x->getFieldPtr(0)->getVFieldPtr()->get0FieldPtr()->getLength() ) ) << "\n";
			for( int i=1; i<=space->nGlo(3); ++i )
				*out << i << "\t" << x->getFieldPtr(0)->getVFieldPtr()->getFieldPtr(i)->norm() << "\t" << std::sqrt( static_cast<S>( x->getFieldPtr(0)->getVFieldPtr()->getFieldPtr(i)->getLength() ) ) << "\n";
		}

		// spectral refinement criterion
		if( refinement>1 ) {

			x->getFieldPtr(0)->getVFieldPtr()->exchange();
			S u_nf = x->getFieldPtr(0)->getVFieldPtr()->getFieldPtr(space->nGlo(3))->norm();
			S u_1  = x->getFieldPtr(0)->getVFieldPtr()->getFieldPtr(1             )->norm();
			S truncError = 1.;
			if( u_nf != u_1 ) // just in case u_1=u_nf=0
			 truncError = u_nf / u_1 ;
			if( truncError < refinementTol ) {
				if( 0==space->rankST() )
					std::cout << "\n||u[nf]||/||u[1]|| = " << truncError << " < " << refinementTol << "\n\n"; 
				break;
			}	
			else
				if( 0==space->rankST() )
					std::cout << "\n||u[nf]||/||u[1]|| = " << truncError << " >= " << refinementTol << "\n\n"; 

			auto spaceF =
				Pimpact::RefinementStrategy<SpaceT>::createRefinedSpace(
						space, Teuchos::tuple<int>( 0, 0, 0, refinementStep ) );

			auto refineOp =
				Teuchos::rcp(
						new Pimpact::TransferCompoundOp<
						Pimpact::TransferMultiHarmonicOp< Pimpact::VectorFieldOpWrap< Pimpact::InterpolationOp<SpaceT> > >,
						Pimpact::TransferMultiHarmonicOp< Pimpact::InterpolationOp<SpaceT> >
						>( space, spaceF ) );
			// refineOp->print();

			// init Fields for fine Boundary conditions
			auto xf = Pimpact::create<CF>( spaceF );
			xf->getVFieldPtr()->initField( pl->sublist("Base flow") );

			auto temp = Pimpact::create<CF>( spaceF );

			refineOp->apply( x->getField(0), *temp );

			xf->add( 1., *temp, 0., *temp );

			x = Pimpact::createMultiField( xf );
			space = spaceF;
		}
	} // end of for( int refine=0; refine<refinement; ++refine ) {
	/******************************************************************************************/

	Teuchos::TimeMonitor::summarize();

	Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );

	MPI_Finalize();
	return( 0 );

}
