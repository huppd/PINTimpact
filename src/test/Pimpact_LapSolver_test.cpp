#include <fstream>
#include <ostream>

#include <mpi.h>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"
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




namespace {


using ST = double;
using OT = int;

const int dimension = 4;
const int dNC = 4;
const int sd = 3;
//const int dNC = 2;

using SpaceT = Pimpact::Space<ST,OT,sd,dimension,4>;

using FSpaceT = SpaceT;
using CSpaceT = Pimpact::Space<ST,OT,sd,dimension,2>;

using CS = Pimpact::CoarsenStrategyGlobal<FSpaceT,CSpaceT>;

using MGSpacesT = Pimpact::MGSpaces<FSpaceT,CSpaceT>;


template<class T>
	using BSF = Pimpact::MultiField< Pimpact::ScalarField<T> >;

template<class T> using BOPF = Pimpact::MultiOpWrap< Pimpact::DivGradOp<T> >;
template<class T> using BOPC = Pimpact::MultiOpWrap< Pimpact::DivGradO2Op<T> >;
template<class T> using BSM = Pimpact::MultiOpWrap< Pimpact::DivGradO2JSmoother<T> >;

template<class T> using ConvDiffOpT =
	Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<T> >;

template<class T> using ConvDiffSORT =
	Pimpact::NonlinearSmoother<T,Pimpact::ConvectionDiffusionSORSmoother >;

template<class T> using ConvDiffJT =
	Pimpact::NonlinearSmoother<T,Pimpact::ConvectionDiffusionJSmoother >;

template<class T1,class T2> using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T> using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> >;
template<class T> using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;

template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;
template<class T> using POP = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2JSmoother >;
template<class T> using POP2 = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2SORSmoother >;


using DGJMGT = Pimpact::MultiGrid<
	MGSpacesT,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	//Pimpact::DivGradO2Op,
	Pimpact::DivGradO2Op,
	Pimpact::DivGradO2JSmoother,
	Pimpact::DivGradO2Inv >;

using DGSORMGT = Pimpact::MultiGrid<
	MGSpacesT,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::DivGradO2SORSmoother,
	Pimpact::DivGradO2Inv >;

using DGLMGT = Pimpact::MultiGrid<
	MGSpacesT,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::DivGradO2LSmoother,
	Pimpact::DivGradO2Inv >;

using DGCMGT = Pimpact::MultiGrid<
	MGSpacesT,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::Chebyshev,
	Pimpact::DivGradO2Inv >;


// blup
using VF = Pimpact::MultiHarmonicField< Pimpact::VectorField<SpaceT> >;
using SF = Pimpact::MultiHarmonicField< Pimpact::ScalarField<SpaceT> >;

using MVF = Pimpact::MultiField<VF>;
using MSF = Pimpact::MultiField<SF>;

using CF = Pimpact::CompoundField< VF, SF>;


template<class T1,class T2>
	using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1,T2> >;
template<class T>
	using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> >;
template<class T>
	using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;


template<class T>
	using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;
template<class T>
	using POP = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2JSmoother >;



std::string xmlFilename = "../../XML/parameterLAP_PF.xml";


TEUCHOS_STATIC_SETUP() {

	Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
	clp.addOutputSetupOptions(true);

	clp.setOption("filename", &xmlFilename, "file name of the input xml paramerterlist");

}




TEUCHOS_UNIT_TEST( bla, bla  ) {

	/////////////////////////////////////////// set up parameters ///////////////////////////
	Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::getParametersFromXmlFile( xmlFilename );
	pl->print();
	/////////////////////////////////////// end of set up parameters /////////////////////////


	///////////////////////////////////////////  set up initial stuff ////////////////////////////
	Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>(
			Teuchos::rcpFromRef( pl->sublist( "Space", true ) ) );


	// init Fields
	auto x = Pimpact::create<CF>( space );
	auto rhs = x->clone( Pimpact::ECopy::Shallow );
	auto rhs2 = x->clone( Pimpact::ECopy::Shallow );

	x->getVFieldPtr()->initField( pl->sublist("Base flow") );
	x->getSFieldPtr()->get0FieldPtr()->initField( Pimpact::Grad2D_inX, -2./space->getDomainSize()->getRe() );
	x->getVFieldPtr()->write();


	// construct NS Operators
	auto opV2V = Pimpact::createMultiDtConvectionDiffusionOp( space );
	auto opS2V = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::GradOp>( space ) );
	auto opV2S = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::DivOp>( space ) );
	opV2V->assignField( x->getVField() );
	auto divGradOp =
		Pimpact::createDivGradOp(
				opV2S->getOperatorPtr(),
				opS2V->getOperatorPtr() );

	// check solution for divergence freeness
	{
		opV2S->apply( x->getVField(), rhs->getSField() );
		ST div = rhs->getSFieldPtr()->norm();
		if( 0==space->rankST() )
			std::cout << "\ndiv(U) = " << div << "\n";
	}

	// generate RHS
	opV2V->apply( x->getVField(), rhs->getVField() );
	opV2S->apply( rhs->getVField(), rhs->getSField() );
	divGradOp->apply( x->getSFieldPtr()->get0Field(), rhs2->getSFieldPtr()->get0Field() );

	{ // check pressure solution connsistent
		auto res = rhs->clone();
		opS2V->apply( x->getSField(), res->getVField() );
		res->getVFieldPtr()->add( 1., rhs->getVField(), 1., res->getVField() );

		//rhs->getVFieldPtr()->write(100);
		//res->getVFieldPtr()->write(200);

		ST residual = res->getVFieldPtr()->norm();
		res->getVFieldPtr()->get0FieldPtr()->write(300);
		if( 0==space->rankST() ) 
			std::cout << "|| rhs_u - grad(p) ||: " << residual << "\n";
	}
	{ // check consistency of RHS
		auto temp = rhs2->getSFieldPtr()->get0FieldPtr()->clone();
		rhs->getSFieldPtr()->get0FieldPtr()->write(  100 );
		rhs2->getSFieldPtr()->get0FieldPtr()->write( 200 );
		temp->add(
				1., rhs2->getSFieldPtr()->get0Field(),
				1., rhs->getSFieldPtr()->get0Field() );
		temp->write(300);
		std::cout << "difference rhs: " << temp->norm( Belos::InfNorm ) << "\n";
	}

	// initial guess
	//x->getSFieldPtr()->init( 0. );
	//x->getSFieldPtr()->random();

	// --- inverse DivGrad

	auto divGradInv2 =
		Pimpact::createInverseOp( 
				divGradOp,
				Teuchos::rcpFromRef( pl->sublist("DivGrad") ) );

	// --- inverse DivGrad^T
	auto divGradInvTT =
		Pimpact::createInverseOp( 
				Pimpact::createTranspose(
					divGradOp ),
				Teuchos::rcpFromRef( pl->sublist("DivGrad") ) );

	//  Multigrid DivGrad
	auto mgSpaces = Pimpact::createMGSpaces<FSpaceT,CSpaceT,CS>( space,
			pl->sublist("Multi Grid").get<int>("maxGrids") );

	auto mgDivGrad = 
		Pimpact::createMultiGrid<
		Pimpact::ScalarField, Pimpact::TransferOp, Pimpact::RestrictionSFOp,
		Pimpact::InterpolationOp,
		Pimpact::DivGradOp,
		Pimpact::DivGradO2Op,
		//Pimpact::DivGradO2JSmoother,
		//Pimpact::DivGradO2SORSmoother,
		//Pimpact::DivGradO2LSmoother,
		Pimpact::Chebyshev,
		Pimpact::DivGradO2Inv
			>( mgSpaces, Teuchos::rcpFromRef( pl->sublist("DivGrad").sublist("Multi Grid") ) );

	if( 0==space->rankST() )
		mgDivGrad->print();

	// set preconditioners
	divGradInv2->setLeftPrec( Pimpact::createMultiOperatorBase(
				Pimpact::createInvDiagonal( divGradOp ) ) );
	//divGradInv2->setRightPrec( Pimpact::createMultiOperatorBase(
	//Pimpact::createInvDiagonal( divGradOp ) ) );
	divGradInv2->setRightPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );
	//divGradInv2->setLeftPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );
	//divGradInvTT->setLeftPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );
	//divGradInvTT->setLeftPrec( Pimpact::createMultiOperatorBase(

	auto divGradInv =
		Pimpact::createMultiHarmonicMultiOpWrap(
				divGradInv2 );
	auto divGradInvT =
		Pimpact::createMultiHarmonicMultiOpWrap(
				divGradInvTT );


	//for( int i=0; i<10; ++i ) 
	//mgDivGrad->apply(
	//rhs->getSFieldPtr()->get0Field(),
	//x->getSFieldPtr()->get0Field() );

	//{ // compute contribution to nullspace
	//auto fSlevel = rhs->getSFieldPtr()->clone( Pimpact::ECopy::Deep );
	////rhs->getSFieldPtr()->level();
	//fSlevel->level();
	//fSlevel->add( 1., rhs->getSField(), -1., *fSlevel );
	//std::cout << " contribution nullspace: " << fSlevel->norm() << "\n";
	//}
	//rhs->getSFieldPtr()->level();
	//
	// solve nullspace IMPACT style
	/*auto nullspace = rhs->getSFieldPtr()->clone();*/
	//{
	//nullspace->init( 1. );
	////nullspace->random();
	//auto ones = nullspace->clone();
	//ones->init( 0. );
	////ones->random();
	//divGradInvT->apply( *ones, *nullspace );
	//nullspace->write(888);
	//}

	// --------------------- solve with plaine rhs
	x->getSFieldPtr()->init( 0. );
	//x->getSFieldPtr()->random();
	divGradInv->apply( rhs2->getSField(), x->getSField() );
	x->getSFieldPtr()->init( 0. );
	//x->getSFieldPtr()->random();
	divGradInv->apply( rhs->getSField(), x->getSField() );

	// write solution
	x->getSFieldPtr()->get0FieldPtr()->write(400);

	// compute residual NS
	{
		//auto res = rhs->clone();
		//opS2V->apply( x->getSField(), res->getVField() );
		//res->getVFieldPtr()->add( 1., rhs->getVField(), -1., res->getVField() );

		//rhs->getVFieldPtr()->write(100);
		//res->getVFieldPtr()->write(200);

		//ST residual = res->getVFieldPtr()->norm();
		//res->getVFieldPtr()->write(300);
		//if( 0==space->rankST() ) 
		//std::cout << "|| rhs_u - grad(p) ||: " << residual << "\n";
	}


	////// compute Residual I
	//

	// --------------------- solve with nullspaced rhs
	//auto xl = x->getSFieldPtr()->clone();
	//xl->init(0.);
	////rhs->getSField().level();
	//auto f_ = rhs->getSFieldPtr()->clone( Pimpact::ECopy::Deep );
	//f_->get0FieldPtr()->add(
	//1., rhs->getSFieldPtr()->get0Field(),
	//-nullspace->get0FieldPtr()->dot( rhs->getSFieldPtr()->get0Field() )/nullspace->get0FieldPtr()->dot( nullspace->get0Field() ), nullspace->get0Field() );
	//std::cout << "(rhs,nullspace): " << f_->get0FieldPtr()->dot( nullspace->get0Field() ) << "\n";
	//divGradInv->apply( *f_, *xl );

	//xl->level();
	//xl->write( 10 );
	//// diff between solution
	//x->getSFieldPtr()->add( 1., x->getSField(), -1., *xl );
	//x->getSFieldPtr()->write(666);
	//std::cout << "||x-xl||: " << x->getSFieldPtr()->norm() << "\n";

	// compute resudial II
	//opS2V->apply( xl->getSField(), x->getVField() );
	//x->getVField().add( 1., xl->getVField(), 1., rhs->getVField() );
	//std::cout << "residual: " << x->getVFieldPtr()->norm() << "\n";

	Teuchos::TimeMonitor::summarize();

	Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );


} // end of TEUCHOS_UNIT_TEST 


} // end of namespace
