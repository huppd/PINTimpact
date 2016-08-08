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
//const int dNC = 2;

using SpaceT = Pimpact::Space<ST,OT,dimension,4>;

using FSpaceT = SpaceT;
using CSpaceT = Pimpact::Space<ST,OT,dimension,2>;

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
template<class T> using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionHWOp<T> >;
template<class T> using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;

template<class T> using MOP = Pimpact::MultiOpUnWrap<Pimpact::InverseOp< Pimpact::MultiOpWrap< T > > >;
template<class T> using POP = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2JSmoother >;
template<class T> using POP2 = Pimpact::PrecInverseOp< T, Pimpact::DivGradO2SORSmoother >;


using DGJMGT = Pimpact::MultiGrid<
	MGSpacesT,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionHWOp,
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
	Pimpact::RestrictionHWOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::DivGradO2SORSmoother,
	Pimpact::DivGradO2Inv >;

using DGLMGT = Pimpact::MultiGrid<
	MGSpacesT,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionHWOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::DivGradO2LSmoother,
	Pimpact::DivGradO2Inv >;

using DGCMGT = Pimpact::MultiGrid<
	MGSpacesT,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionHWOp,
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
using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionHWOp<T> >;
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
	Teuchos::RCP<const SpaceT> space = Pimpact::createSpace<ST,OT,dimension,dNC>(
			Teuchos::rcpFromRef( pl->sublist( "Space", true ) ) );


	// init vectors
	auto x = Pimpact::create<CF>( space );
	auto f = x->clone( Pimpact::ShallowCopy );

	// init Fields
	x->getVFieldPtr()->initField( pl->sublist("Base flow") );
	x->getVFieldPtr()->write();


	auto opV2V = Pimpact::createMultiDtConvectionDiffusionOp( space );
	auto opS2V = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::GradOp>( space ) );
	auto opV2S = Pimpact::createMultiHarmonicOpWrap( Pimpact::create<Pimpact::DivOp>( space ) );
	opV2V->assignField( x->getVField() );

	opV2S->apply( x->getVField(), f->getSField() );
	if( 0==space->rankST() )
		std::cout << "\ndiv(U) = " << f->getSFieldPtr()->norm() << "\n";

	// create Multi space
	auto mgSpaces =
		Pimpact::createMGSpaces<FSpaceT,CSpaceT,CS>(
				space,
				pl->sublist("Multi Grid").get<int>("maxGrids") );


	// --- inverse DivGrad
	auto divGradOp =
		Pimpact::createDivGradOp( opV2S->getOperatorPtr(),
				opS2V->getOperatorPtr() );

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


	auto mgDivGrad = 
		Pimpact::createMultiGrid<
		Pimpact::ScalarField, Pimpact::TransferOp, Pimpact::RestrictionHWOp,
		Pimpact::InterpolationOp,
		Pimpact::DivGradOp,
		Pimpact::DivGradO2Op,
		//Pimpact::DivGradO2LSmoother,
		Pimpact::Chebyshev,
		Pimpact::DivGradO2Inv
			>( mgSpaces, Teuchos::rcpFromRef( pl->sublist("DivGrad").sublist("Multi Grid") ) );

	if( 0==space->rankST() )
		mgDivGrad->print();

	{ // without prec

		divGradInv2->setLeftPrec( Pimpact::createMultiOperatorBase(
					Pimpact::createInvDiagonal( divGradOp ) ) );
		//divGradInv2->setRightPrec( Pimpact::createMultiOperatorBase(
					//Pimpact::createInvDiagonal( divGradOp ) ) );

		//divGradInv2->setRightPrec( Pimpact::createMultiOperatorBase( mgDivGrad ) );

		//divGradInvTT->setLeftPrec( Pimpact::createMultiOperatorBase( mgDivGrad )
		//);
		
		//divGradInvTT->setLeftPrec( Pimpact::createMultiOperatorBase(
					//Pimpact::createInvDiagonal( divGradOp ) ) );

		auto divGradInv =
			Pimpact::createMultiHarmonicMultiOpWrap(
					divGradInv2 );
		auto divGradInvT =
			Pimpact::createMultiHarmonicMultiOpWrap(
					divGradInvTT );

		opV2V->apply( x->getVField(), f->getVField() );
		opV2S->apply( f->getVField(), f->getSField() );

		f->getSFieldPtr()->write( 99 );
		{ // compute contribution to nullspace
			auto fSlevel = f->getSFieldPtr()->clone( Pimpact::DeepCopy );
			//f->getSFieldPtr()->level();
			fSlevel->level();
			fSlevel->add( 1., f->getSField(), -1., *fSlevel );
			std::cout << " contribution nullspace: " << fSlevel->norm() << "\n";
		}

		//x->getSFieldPtr()->random();
		//x->getSFieldPtr()->get0FieldPtr()->initField( Pimpact::Grad2D_inX, -2./space->getDomainSize()->getRe() ); //for( int i=0; i<10; ++i ) //mgDivGrad->apply( f->getSFieldPtr()->get0Field(), x->getSFieldPtr()->get0Field() );

		//f->getSFieldPtr()->level();
		//
		// solve nullspace IMPACT style
		auto nullspace = f->getSFieldPtr()->clone();
		{
		nullspace->init( 1. );
		//nullspace->random();
		auto ones = nullspace->clone();
		ones->init( 0. );
		//ones->random();
		divGradInvT->apply( *ones, *nullspace );
		nullspace->write(888);
		}

		// --------------------- solve with plaine rhs
		x->getSFieldPtr()->init( 0. );
		divGradInv->apply( f->getSField(), x->getSField() );

		//// compute Residual I
		auto res =  Pimpact::createMultiField( f->getSFieldPtr()->get0FieldPtr()->clone() );
		//res->push_back();
		//res->push_back();
		divGradInv->getOperatorPtr()->getLinearProblem()->getProblem()->getOperator()->apply(
				*Pimpact::createMultiField( x->getSFieldPtr()->get0FieldPtr() ),
				*res );
		res->getFieldPtr( 0 )->add( 1., f->getSFieldPtr()->get0Field(), -1., res->getField( 0 ) );
		std::cout << "||res||: " << res->norm() << "\n";
		res->write( 777 );
		//res->push_back();
		//divGradInv->getOperatorPtr()->getLinearProblem()->getProblem()->computeCurrResVec( res.getRawPtr() );


		x->level();
		x->getSFieldPtr()->write();

		// --------------------- solve with nullspaced rhs
		auto xl = x->getSFieldPtr()->clone();
		xl->init(0.);
		//f->getSField().level();
		auto f_ = f->getSFieldPtr()->clone( Pimpact::DeepCopy );
		f_->get0FieldPtr()->add(
				1., f->getSFieldPtr()->get0Field(),
				-nullspace->get0FieldPtr()->dot( f->getSFieldPtr()->get0Field() )/nullspace->get0FieldPtr()->dot( nullspace->get0Field() ), nullspace->get0Field() );
		std::cout << "(rhs,nullspace): " << f_->get0FieldPtr()->dot( nullspace->get0Field() ) << "\n";
		divGradInv->apply( *f_, *xl );

		xl->level();
		xl->write( 10 );
		// diff between solution
		x->getSFieldPtr()->add( 1., x->getSField(), -1., *xl );
		x->getSFieldPtr()->write(666);
		std::cout << "||x-xl||: " << x->getSFieldPtr()->norm() << "\n";

		// compute resudial II
		//opS2V->apply( xl->getSField(), x->getVField() );
		//x->getVField().add( 1., xl->getVField(), 1., f->getVField() );
		//std::cout << "residual: " << x->getVFieldPtr()->norm() << "\n";

	}
	Teuchos::TimeMonitor::summarize();

	Teuchos::writeParameterListToXmlFile( *pl, "parameterOut.xml" );

	//MPI_Finalize();
	//return( 0 );

} // end of TEUCHOS_UNIT_TEST 

} // end of namespace
