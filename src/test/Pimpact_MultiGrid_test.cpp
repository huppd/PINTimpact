#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Pimpact_CoarsenStrategy.hpp"
#include "Pimpact_CoarsenStrategyGlobal.hpp"
#include "Pimpact_Fields.hpp"
#include "Pimpact_LinearProblem.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_MultiGrid.hpp"
#include "Pimpact_Operator.hpp"

#include "Pimpact_Test.hpp"




namespace {



using CSpace2DT = Pimpact::Space<ST,OT,2,d,2>;
using CSpace3DT = Pimpact::Space<ST,OT,3,d,2>;

using CSL2D = Pimpact::CoarsenStrategy<D2,CSpace2DT>;
using CSG2D = Pimpact::CoarsenStrategyGlobal<D2,CSpace2DT>;

using CSL3D = Pimpact::CoarsenStrategy<D3,CSpace3DT>;
using CSG3D = Pimpact::CoarsenStrategyGlobal<D3,CSpace3DT>;

using Spaces2D = Pimpact::MGSpaces<D2,CSpace2DT>;
using Spaces3D = Pimpact::MGSpaces<D3,CSpace3DT>;



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


using DGJMGT2D = Pimpact::MultiGrid<
	Spaces2D,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	//Pimpact::DivGradO2Op,
	Pimpact::DivGradO2Op,
	Pimpact::DivGradO2JSmoother,
	//MOP >;
	//POP >;
	//Pimpact::DivGradO2JSmoother >;
	Pimpact::DivGradO2Inv >;

using DGJMGT3D = Pimpact::MultiGrid<
	Spaces3D,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	//Pimpact::DivGradO2Op,
	Pimpact::DivGradO2Op,
	Pimpact::DivGradO2JSmoother,
	//MOP >;
	//POP >;
	//Pimpact::DivGradO2JSmoother >;
	Pimpact::DivGradO2Inv >;



using DGLMGT2D = Pimpact::MultiGrid<
	Spaces2D,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::DivGradO2LSmoother,
	MOP >;
//Pimpact::DivGradO2JSmoother >;
//Pimpact::DivGradO2LSmoother >;
//Pimpact::DivGradO2Inv >;

using DGLMGT3D = Pimpact::MultiGrid<
	Spaces3D,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::DivGradO2LSmoother,
	MOP >;
//Pimpact::DivGradO2JSmoother >;
//Pimpact::DivGradO2LSmoother >;
//Pimpact::DivGradO2Inv >;

using DGCMGT2D = Pimpact::MultiGrid<
	Spaces2D,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::Chebyshev,
	MOP >;
//Pimpact::Chebyshev >;
//Pimpact::DivGradO2Inv >;

using DGCMGT3D = Pimpact::MultiGrid<
	Spaces3D,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::Chebyshev,
	MOP >;
//Pimpact::Chebyshev >;
//Pimpact::DivGradO2Inv >;



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGSpaces, constructor, CS ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	int rank = space->rankST();
	if( print ) space->print();

	{
		Teuchos::RCP<std::ostream> fstream = Pimpact::createOstream( "coord_rank_"+std::to_string(rank)+".txt", 0 );

		space->getCoordinatesGlobal()->print( *fstream );
	}

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids);

	std::cout << "rank: " << space->rankST() << "\tnGridLevels: " << mgSpaces->getNGrids() << "\n";

	if( space->rankST()==0 && print )
		mgSpaces->print();

	for( int level=0; level<mgSpaces->getNGrids(); ++level ) {
		Teuchos::RCP<std::ostream> fstream = Pimpact::createOstream( "coord_l"+std::to_string(level)+"rank_"+std::to_string(rank)+".txt",0 );
		mgSpaces->get(level)->getCoordinatesGlobal()->print( *fstream );
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSpaces, constructor, CSG2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSpaces, constructor, CSG3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGFields, SFconstructor, CS ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids );

	std::cout << "nGridLevels: " << mgSpaces->getNGrids() << "\n";
	if( space->rankST()==0 && print )
		mgSpaces->print();

	auto mgFields = Pimpact::createMGFields<Pimpact::ScalarField>( mgSpaces );

	auto field = mgFields->get( -1 );
	if( mgSpaces->participating(-1) ){

		field->init(5.);
		TEST_FLOATING_EQUALITY( std::sqrt( std::pow(5.,2)*field->getLength() ), field->norm(Belos::TwoNorm), eps );

		if( write ) field->write(0);
	}

}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, SFconstructor, CSL )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, SFconstructor, CSG2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, SFconstructor, CSG3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGFields, VFconstructor, CS ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids );

	std::cout << "nGridLevels: " << mgSpaces->getNGrids() << "\n";
	if( space->rankST()==0 && print )
		mgSpaces->print();

	auto mgFields = Pimpact::createMGFields<Pimpact::VectorField>( mgSpaces );

	auto field = mgFields->get( -1 );
	if( mgSpaces->participating(-1) ){

		field->init(5.);
		TEST_FLOATING_EQUALITY( std::sqrt( std::pow(5.,2)*field->getLength() ), field->norm(Belos::TwoNorm), eps );

		if( write ) field->write(0);
	}

}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, VFconstructor, CSL )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, VFconstructor, CSG2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGFields, VFconstructor, CSG3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGOperators, SFconstructor, CS ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids );

	auto mgOps = Pimpact::createMGOperators<Pimpact::DivGradO2Op>( mgSpaces );

	auto mgOps2 = Pimpact::createMGOperators<Pimpact::DivGradOp,Pimpact::DivGradO2Op>( mgSpaces );

	for( int i=0; i<3; ++i ) {
		auto op = mgOps->get( i );
		if( 0==space->rankST() )
			std::cout << "\n\n --- level: " << i << " ---\n\n";

		if( mgSpaces->participating(i) && print )
			op->print();
	}

}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, SFconstructor, CSL )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, SFconstructor, CSG2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, SFconstructor, CSG3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGOperators, VFconstructor, CS ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids );

	auto mgOps = Pimpact::createMGOperators<ConvDiffOpT>( mgSpaces );


	auto op = mgOps->get( -1 );

	if( mgSpaces->participating(-1) && print )
		op->print();

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, VFconstructor, CSG2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGOperators, VFconstructor, CSG3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGSmoothers, SFconstructor, CS ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids );


	auto mgOps = Pimpact::createMGOperators<Pimpact::DivGradOp,Pimpact::DivGradO2Op>( mgSpaces );

	auto mgSmoother = Pimpact::createMGSmoothers<Pimpact::DivGradO2JSmoother>( mgOps );

	auto op = mgSmoother->get( -2 );

	if( mgSpaces->participating(-2) && print )
		op->print();

}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, SFconstructor, CSL )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, SFconstructor, CSG2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, SFconstructor, CSG3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGSmoothers, VFconstructor, CS ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids );

	auto mgOps = Pimpact::createMGOperators<ConvDiffOpT>( mgSpaces );

	auto mgSmoother = Pimpact::createMGSmoothers<ConvDiffSORT>( mgOps );

	auto op = mgSmoother->get( -2 );

	if( mgSpaces->participating(-2) && print )
		op->print();

}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, VFconstructor, CSL )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, VFconstructor, CSG2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGSmoothers, VFconstructor, CSG3D )



TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MGTransfers, Restrictor, CS, RestrictorType ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids );

	for( int level=1; level<mgSpaces->getNGrids(); ++level ) {
		if( 0==space->rankST() ) {
			std::cout << "\n\n\t--- level: " << level-1 << "---\n";
		}
		Teuchos::RCP<const RestrictorType > op = 
			Teuchos::rcp( new RestrictorType(
						mgSpaces->get(level-1),
						mgSpaces->get(level),
						mgSpaces->get()->getProcGrid()->getNP() ) );

		if( mgSpaces->participating(level-1) && print ) op->print();

		std::vector<Pimpact::F> types;
		if( "Restriction SF"==op->getLabel() )
			types.push_back( Pimpact::F::S );
		else if( "Restriction VF"==op->getLabel() ) {
			types.push_back( Pimpact::F::U );
			types.push_back( Pimpact::F::V );
			types.push_back( Pimpact::F::W );
		}

		for( auto type=types.begin(); type!=types.end(); ++type ) {
			if( 2==CS::SpaceT::sdim && *type==Pimpact::F::W ) break;
			if( 0==space->rankST() )
				std::cout << " --- ftype: " << Pimpact::toString(*type) << " ---\n";

			Teuchos::RCP< Pimpact::ScalarField<typename CS::CSpaceT> > fieldf;
			Teuchos::RCP< Pimpact::ScalarField<typename CS::CSpaceT> > fieldc;
			Teuchos::RCP< Pimpact::ScalarField<typename CS::CSpaceT> > sol;
			Teuchos::RCP< Pimpact::ScalarField<typename CS::CSpaceT> > er;

			fieldf = Pimpact::createScalarField( mgSpaces->get( level-1 ), *type );

			fieldc = Pimpact::createScalarField( mgSpaces->get( level ), *type );
			sol = fieldc->clone();
			er = fieldc->clone();

			// the zero test
			fieldf->initField( Pimpact::ConstField, 0. );
			fieldc->init( 1. );

			if( mgSpaces->participating(level-1) ) op->apply( *fieldf, *fieldc );

			if( mgSpaces->participating(level-1) )
				TEST_FLOATING_EQUALITY( 0., fieldf->norm(), eps );

			if( mgSpaces->participating(level) ) 
				TEST_FLOATING_EQUALITY( 0., fieldc->norm(), eps );

			// the random test
			fieldf->random();

			if( mgSpaces->participating(level-1) ) op->apply( *fieldf, *fieldc );

			if( mgSpaces->participating(level) )
				TEST_INEQUALITY( 0., fieldc->norm() );

			// the const test
			fieldf->initField( Pimpact::ConstField, 1. );
			fieldc->initField( Pimpact::ConstField, 0. );
			sol->initField( Pimpact::ConstField, 1. );

			if( mgSpaces->participating(level-1) ) op->apply( *fieldf, *fieldc );

			if( mgSpaces->participating(level) ) {

				er->add( 1., *sol, -1., *fieldc, Pimpact::B::Y );
				ST errInf = er->norm(Belos::InfNorm, Pimpact::B::Y);
				if( 0==space->rankST() )
					std::cout << "error Const: " << errInf << " ("<< op->getDD() << ")\n";
				//if( i>0 )
				TEST_EQUALITY( errInf<eps, true ); // boundaries?
				if( errInf>=eps )
					if( write ) er->write(0);
			}

			// the hard test
			for( int dir=1; dir<=3; ++ dir ) {

				Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);
				fieldf->initField( type );
				fieldc->initField( Pimpact::ConstField, 0. );
				sol->initField( type );

				if( mgSpaces->participating(level-1) ) op->apply( *fieldf, *fieldc );

				if( mgSpaces->participating(level) ) {
					er->add( 1., *sol, -1., *fieldc, Pimpact::B::Y );
					double errInf = er->norm(Belos::InfNorm, Pimpact::B::Y );
					if( 0==space->rankST() )
						std::cout << "error ("<< toString(type) << "): " << errInf << " ("<< op->getDD() << ")\n";
					TEST_EQUALITY( errInf<eps, true );
					if( errInf>=eps ) {
						if( print ) er->print();
						if( write ) {
							er->write(dir);
							fieldc->write(10*(dir));
							sol->write(100*(dir));
						}
					}
				}
			}
		}
	} // end of for( )

} // end of TEUCHOS_UNIT_TEST_TEMPLATE


using ResSF2D = Pimpact::RestrictionSFOp<CSpace2DT>;
using ResVF2D = Pimpact::RestrictionVFOp<CSpace2DT>;
using ResSF3D = Pimpact::RestrictionSFOp<CSpace3DT>;
using ResVF3D = Pimpact::RestrictionVFOp<CSpace3DT>;

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MGTransfers, Restrictor, CSG2D, ResSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MGTransfers, Restrictor, CSG2D, ResVF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MGTransfers, Restrictor, CSG3D, ResSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MGTransfers, Restrictor, CSG3D, ResVF3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGTransfers, Interpolator, CS ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids );

	std::cout << "\nrank: " << space->rankST() << "\tnGridLevels: " << mgSpaces->getNGrids() << "\n";

	auto mgTransfers = Pimpact::createMGTransfers<Pimpact::TransferOp,Pimpact::RestrictionSFOp,Pimpact::InterpolationOp>( mgSpaces );

	for( int level=0; level<mgSpaces->getNGrids()-1; ++level ) {

		Teuchos::RCP< const Pimpact::InterpolationOp<typename CS::CSpaceT> > op = 
			Teuchos::rcp( new Pimpact::InterpolationOp<typename CS::CSpaceT>(
						mgSpaces->get(level+1),
						mgSpaces->get(level),
						mgSpaces->get()->getProcGrid()->getNP() ) );
		if( 0==space->rankST() ) {
			std::cout << "\n\n\t--- level: " << level << "---\n";
		}
		if( space->rankST()==rankbla && print ) {
			op->print();
		}

		Teuchos::Tuple<Pimpact::F,4> type =
			Teuchos::tuple(
					Pimpact::F::S,
					Pimpact::F::U,
					Pimpact::F::V,
					Pimpact::F::W );

		for( int i=fs; i<fe; ++i ) {
			if( 0==space->rankST() )
				std::cout << "field type: " << i << "\n";

			auto fieldf = Pimpact::createScalarField( mgSpaces->get( level   ), type[i] );
			auto fieldc = Pimpact::createScalarField( mgSpaces->get( level+1 ), type[i] );
			auto sol = fieldf->clone();
			auto er = fieldf->clone();


			// the zero test
			fieldf->init( 1. );
			fieldc->initField( Pimpact::ConstField, 0. );

			if( mgSpaces->participating(level) )
				op->apply( *fieldc, *fieldf );

			if( mgSpaces->participating(level) )
				TEST_EQUALITY( eps>fieldf->norm(Belos::InfNorm), true );
			if( mgSpaces->participating(level+1) )
				TEST_EQUALITY( eps>fieldc->norm(Belos::InfNorm), true );

			// the random test
			fieldc->random();
			fieldf->init(0.);

			if( mgSpaces->participating(level+1) )
				TEST_INEQUALITY( 0., fieldc->norm() );

			if( mgSpaces->participating(level) )
				op->apply( *fieldc, *fieldf );

			if( mgSpaces->participating(level+1) )
				TEST_INEQUALITY( 0., fieldc->norm() );


			// the stronger init test
			fieldc->initField( Pimpact::ConstField, 1. );
			fieldf->initField( Pimpact::ConstField, 0. );
			sol->initField( Pimpact::ConstField, 1. );
			er->random();

			if( mgSpaces->participating(level) )
				op->apply( *fieldc, *fieldf );


			er->add( 1., *sol, -1., *fieldf );
			if( mgSpaces->participating(level) ) {
				ST errInf = er->norm( Belos::InfNorm );
				if( 0==space->rankST() )
					std::cout << "error Const: " << errInf << "\n";
				TEST_EQUALITY( errInf < eps, true  );
				if( errInf>=eps )
					if( write ) er->write(0);
			}
			if( er->normLoc(Belos::InfNorm)>=eps && space->rankST()==rankbla && print ){
				//			std::cout << "rank: " << space->rankST() << "\n";
				er->print();
			}


			// --- hardcore test ---
			for( int dir=1; dir<=3; ++dir ) {

				Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

				fieldc->initField( type );
				fieldf->initField( Pimpact::ConstField, 0. );
				sol->initField( type );
				er->random();

				if( mgSpaces->participating(level) ){
					op->apply( *fieldc, *fieldf );

					er->add( 1., *sol, -1., *fieldf );

					ST errInf = er->norm( Belos::InfNorm );
					if( 0==space->rankST() )
						std::cout << "error (" << Pimpact::toString(type) << "): " << errInf << "\n";
					TEST_EQUALITY( errInf<eps, true  );
					if( errInf>=eps ) {
						int i = dir;
						if( write ) {
							er->write( i );
							fieldf->write( 10*i );
							sol->write( 100*i );
						}
						if( print ) er->print();
					}
				}
			}
		}
	}

}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, Interpolator, CSL )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, Interpolator, CSG2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, Interpolator, CSG3D )




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGTransfers, MGTransfersSF, CS ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids );

	auto mgTransfers = Pimpact::createMGTransfers<
		Pimpact::TransferOp,Pimpact::RestrictionSFOp,Pimpact::InterpolationOp>(
				mgSpaces );

	auto x = Pimpact::createMGFields<Pimpact::ScalarField>( mgSpaces );

	auto fieldf = x->get();
	auto fieldc = x->get(-1);

	ST errInf = 0.;

	// interpolation 
	{
		auto sol = fieldf->clone( Pimpact::ECopy::Shallow );
		auto er = fieldf->clone( Pimpact::ECopy::Shallow );

		// the zero test
		fieldf->init( 1. );
		fieldc->initField();

		mgTransfers->interpolation( x );

		if( mgSpaces->participating(-1) ) {
			errInf = fieldc->norm( Belos::InfNorm );
			TEST_EQUALITY( errInf<eps, true );
		}
		if( mgSpaces->participating(0) ) {
			errInf = fieldf->norm( Belos::InfNorm );
			TEST_EQUALITY( errInf<eps, true );
			if( 0==space->rankST() )
				std::cout << "\ninterpolation error zero: " << errInf << "\n";
		}

		// the random test
		fieldc->random();
		fieldf->init(0.);

		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm( Belos::InfNorm ) );

		mgTransfers->interpolation( x );

		if( mgSpaces->participating(0) )
			TEST_INEQUALITY( 0., fieldf->norm( Belos::InfNorm ) );

		// the Const test
		fieldc->initField( Pimpact::ConstField, 1. );
		fieldf->initField();
		sol->initField( Pimpact::ConstField, 1. );

		mgTransfers->interpolation( x );

		er->add( 1., *sol, -1., *fieldf );

		if( mgSpaces->participating(0) ) {
			errInf = er->norm( Belos::InfNorm );
			if( 0==space->rankST() )
				std::cout << "interpolation error Const: " << errInf << "\n";
			TEST_EQUALITY( errInf<eps, true  );
			if( errInf>=eps ) {
				int i = 0;
				if( write ) er->write( std::abs(i) );
				std::string r = std::to_string( static_cast<long long>( space->rankST() ) ); // long long needed on brutus(intel)
				if( print ) er->print( *Pimpact::createOstream( "int_error_c_r"+r+".txt" ) );
			}
		}


		// the Grad test
		for( int dir=1; dir<=3; ++dir ) {

			Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

			fieldc->initField( type );
			fieldf->initField();
			sol->initField( type );

			mgTransfers->interpolation( x );

			er->add( 1., *sol, -1., *fieldf );

			if( mgSpaces->participating(0) ) {
				errInf = er->norm( Belos::InfNorm );
				if( 0==space->rankST() )
					std::cout << "interpolation error (" << Pimpact::toString(type) << "): " << errInf << "\n";
				TEST_EQUALITY( errInf<eps, true  );
				if( errInf>=eps ) {
					int i = dir-2;
					if( write ) er->write( std::abs(i) );
					std::string d = std::to_string( static_cast<long long>( dir-1 ) ); // long long needed on brutus(intel)
					std::string r = std::to_string( static_cast<long long>( space->rankST() ) ); // long long needed on brutus(intel)
					if( print ) er->print( *Pimpact::createOstream( "int_error_g"+d+"_r"+r+".txt" ) );
				}
			}
		}

	}

	// restriction 
	{
		auto sol = fieldc->clone( Pimpact::ECopy::Shallow );
		auto er = fieldc->clone( Pimpact::ECopy::Shallow );


		// the zero test
		fieldf->initField();
		fieldc->init( 1. );

		mgTransfers->restriction( x );

		if( mgSpaces->participating(0) )
			TEST_EQUALITY( fieldf->norm()<eps, true );
		if( mgSpaces->participating(-1) )
			TEST_EQUALITY( fieldc->norm()<eps, true );

		// the random test
		fieldf->random();
		fieldc->init(0.);

		if( mgSpaces->participating(0) )
			TEST_INEQUALITY( 0., fieldf->norm() );

		mgTransfers->restriction( x );

		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm() );

		// the const test
		fieldf->initField( Pimpact::ConstField, 1. );
		fieldc->initField();
		sol->initField( Pimpact::ConstField, 1. );

		mgTransfers->restriction( x );

		er->add( 1., *sol, -1., *fieldc );

		if( mgSpaces->participating(-1) ) {
			errInf = er->norm( Belos::InfNorm );
			TEST_EQUALITY( errInf<eps, true  );
		}

		if( 0==space->rankST() )
			std::cout << "\nrestriction error Const: " << errInf << "\n";


		// hardcore Grad test
		for( int dir=1; dir<=3; ++ dir ) {

			Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

			fieldf->initField( type );
			fieldc->initField();

			sol->initField( type );

			mgTransfers->restriction( x );

			er->add( 1., *sol, -1., *fieldc );
			if( mgSpaces->participating(-1) ) {

				errInf = er->norm( Belos::InfNorm );
				TEST_EQUALITY( errInf<eps, true  );
				if( errInf>=eps )
					if( write ) er->write( std::abs(1*(dir-2)) );
			}
			if( 0==space->rankST() )
				std::cout << "restriction error (" << Pimpact::toString(type) << "): " << errInf << "\n";
		}
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, MGTransfersSF, CSG2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, MGTransfersSF, CSG3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MGTransfers, MGTransfersVF, CS ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids );

	auto mgTransfers = Pimpact::createMGTransfers<
		TransVF,RestrVF,InterVF>( mgSpaces );

	auto x = Pimpact::createMGFields<Pimpact::VectorField>( mgSpaces );

	auto fieldf = x->get();
	auto fieldc = x->get(-1);

	// interpolation 
	{
		auto sol = fieldf->clone( Pimpact::ECopy::Shallow );
		auto er = fieldf->clone( Pimpact::ECopy::Shallow );


		// the zero test
		fieldf->init( 1. );
		fieldc->initField();

		mgTransfers->interpolation( x );

		if( mgSpaces->participating(0) )
			TEST_EQUALITY( fieldf->norm()<eps, true );
		if( mgSpaces->participating(-1) )
			TEST_EQUALITY( fieldc->norm()<eps, true );


		// the random test
		fieldc->random();
		fieldf->init(0.);

		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm() );

		mgTransfers->interpolation( x );

		if( mgSpaces->participating(0) )
			TEST_INEQUALITY( 0., fieldf->norm() );

		// the stronger init test
		(*fieldc)(Pimpact::F::U).initField( Pimpact::ConstField, 1. );
		(*fieldc)(Pimpact::F::V).initField( Pimpact::ConstField, 1. );
		(*fieldc)(Pimpact::F::W).initField( Pimpact::ConstField, 1. );
		fieldf->init(0.);
		(*sol)(Pimpact::F::U).initField( Pimpact::ConstField, 1. );
		(*sol)(Pimpact::F::V).initField( Pimpact::ConstField, 1. );
		(*sol)(Pimpact::F::W).initField( Pimpact::ConstField, 1. );

		mgTransfers->interpolation( x );

		er->add( 1., *sol, -1., *fieldf );

		if( mgSpaces->participating(0) ) {
			ST rel_error = er->norm( Belos::InfNorm );
			if( 0==space->rankST() )
				std::cout << "\nint. error Const: " << rel_error << "\n";
			TEST_EQUALITY( rel_error < eps, true  );
			if( rel_error>=eps || isnan(rel_error) ) {
				if( write ) er->write(0);
			}
		}


		//// hardcore test init test in X
		//fieldc->getField( Pimpact::F::U ).initField( Pimpact::Grad2D_inX );
		//fieldc->getField( Pimpact::F::V ).initField( Pimpact::Grad2D_inX );
		//fieldc->getField( Pimpact::F::W ).initField( Pimpact::Grad2D_inX );
		//fieldf->initField( Pimpact::ConstFlow, 0., 0., 0. );

		//sol->getField( Pimpact::F::U )->initField( Pimpact::Grad2D_inX );
		//sol->getField( Pimpact::F::V )->initField( Pimpact::Grad2D_inX );
		//sol->getField( Pimpact::F::W )->initField( Pimpact::Grad2D_inX );

		//mgTransfers->interpolation( x );

		//er->add( 1., *sol, -1., *fieldf );

		//if( mgSpaces->participating(0) ) {
		//S rel_error = er->norm( Belos::InfNorm );
		//if( 0==space->rankST() )
		//std::cout << "int. error GradX: " << rel_error << "\n";
		//TEST_EQUALITY( rel_error<eps, true  );
		//if( rel_error>=eps || isnan(rel_error) ) {
		//er->write(1);
		//fieldf->write(10);
		//sol->write(100);
		//}
		//}

		//// hardcore test init test in Y
		//fieldc->getField( Pimpact::F::U ).initField( Pimpact::Grad2D_inY );
		//fieldc->getField( Pimpact::F::V ).initField( Pimpact::Grad2D_inY );
		//fieldc->getField( Pimpact::F::W ).initField( Pimpact::Grad2D_inY );
		//fieldf->initField( Pimpact::ConstFlow, 0., 0., 0. );

		//sol->getField( Pimpact::F::U ).initField( Pimpact::Grad2D_inY );
		//sol->getField( Pimpact::F::V ).initField( Pimpact::Grad2D_inY );
		//sol->getField( Pimpact::F::W ).initField( Pimpact::Grad2D_inY );

		//mgTransfers->interpolation( x );

		//er->add( 1., *sol, -1., *fieldf );
		//if( mgSpaces->participating(0) ) {
		//S rel_error = er->norm( Belos::InfNorm );
		//if( 0==space->rankST() )
		//std::cout << "int. error GradY: " << rel_error << "\n";
		//TEST_EQUALITY( rel_error<eps, true  );
		//if( rel_error>=eps || isnan(rel_error) ) {
		//er->write(2);
		//fieldf->write(20);
		//sol->write(200);
		//}
		//}

		// hardcore Grad test
		for( int dir=1; dir<=3; ++ dir ) {

			Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

			(*fieldc)( Pimpact::F::U ).initField( type );
			(*fieldc)( Pimpact::F::V ).initField( type );
			(*fieldc)( Pimpact::F::W ).initField( type );
			fieldf->initField();

			(*sol)( Pimpact::F::U ).initField( type );
			(*sol)( Pimpact::F::V ).initField( type );
			(*sol)( Pimpact::F::W ).initField( type );

			mgTransfers->interpolation( x );

			er->add( 1., *sol, -1., *fieldf );

			if( mgSpaces->participating(0) ) {
				ST rel_error = er->norm( Belos::InfNorm );
				if( 0==space->rankST() )
					std::cout << "interpolation error (" << Pimpact::toString(type) << "): " << rel_error << "\n";
				TEST_EQUALITY( rel_error<eps, true  );
				if( rel_error>=eps || isnan(rel_error) ) {
					if( write ) er->write(1*(dir-2));
					if( write ) fieldf->write(10*(dir-2));
					if( write ) sol->write(100*(dir-2));
				}
			}
		}

	}
	// restriction 
	{
		auto sol = fieldc->clone( Pimpact::ECopy::Shallow );
		auto er = fieldc->clone( Pimpact::ECopy::Shallow );

		// the zero test
		fieldf->initField();
		fieldc->init( 1. );

		mgTransfers->restriction( x );

		if( mgSpaces->participating(0) )
			TEST_EQUALITY( fieldf->norm()<eps, true );
		if( mgSpaces->participating(-1) )
			TEST_EQUALITY( fieldc->norm()<eps, true );


		// the random test
		fieldf->random();
		fieldc->init(0.);

		if( mgSpaces->participating(0) )
			TEST_INEQUALITY( 0., fieldf->norm() );

		mgTransfers->restriction( x );

		if( mgSpaces->participating(-1) )
			TEST_INEQUALITY( 0., fieldc->norm() );


		// the stronger init test
		(*fieldf)(Pimpact::F::U).initField( Pimpact::ConstField, 1. );
		(*fieldf)(Pimpact::F::V).initField( Pimpact::ConstField, 1. );
		(*fieldf)(Pimpact::F::W).initField( Pimpact::ConstField, 1. );
		fieldc->init(0.);
		(*sol)(Pimpact::F::U).initField( Pimpact::ConstField, 1. );
		(*sol)(Pimpact::F::V).initField( Pimpact::ConstField, 1. );
		(*sol)(Pimpact::F::W).initField( Pimpact::ConstField, 1. );

		mgTransfers->restriction( x );

		er->add( 1., *sol, -1., *fieldc );
		if( mgSpaces->participating(-1) ) {
			ST rel_error = er->norm()/std::sqrt( (ST)er->getLength() );
			if( 0==space->rankST() )
				std::cout << "res. error Const: " << rel_error << "\n";
			TEST_EQUALITY( rel_error<eps, true  );
			if( rel_error>eps ) {
				if( write ) er->write(0);
			}
		}

		// hardcore grad test
		for( int dir=1; dir<=3; ++ dir ) {

			Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

			(*fieldf)( Pimpact::F::U ).initField( type );
			(*fieldf)( Pimpact::F::V ).initField( type );
			(*fieldf)( Pimpact::F::W ).initField( type );
			fieldc->initField();

			(*sol)( Pimpact::F::U ).initField( type );
			(*sol)( Pimpact::F::V ).initField( type );
			(*sol)( Pimpact::F::W ).initField( type );

			mgTransfers->restriction( x );

			er->add( 1., *sol, -1., *fieldc );
			if( mgSpaces->participating(-1) ) {
				ST rel_error = er->norm()/std::sqrt( (ST)er->getLength() );
				if( 0==space->rankST() )
					std::cout << "restriction grad  error (" << Pimpact::toString(type) << "): " << rel_error << "\n";
				TEST_EQUALITY( rel_error<eps, true  );
				if( rel_error>eps ) {
					if( write ) {
						er->write( dir );
						fieldf->write( 10*(dir) );
						sol->write( 100*(dir) );
					}
				}
			}
		}

	}

}


//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, MGTransfersVF, CSL )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, MGTransfersVF, CSG2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MGTransfers, MGTransfersVF, CSG3D )



TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiGrid, DivGradOp, CS, MGT ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	auto x = Pimpact::create<Pimpact::ScalarField>( space );
	auto b = Pimpact::create<Pimpact::ScalarField>( space );
	auto res = x->clone( Pimpact::ECopy::Shallow );
	auto sol = x->clone( Pimpact::ECopy::Shallow );

	auto xm = Pimpact::createMultiField( x );
	auto bm = Pimpact::createMultiField( b );

	auto op   = Pimpact::create<Pimpact::DivGradOp>( space );
	auto opO2 = Pimpact::create<Pimpact::DivGradO2Op>( space );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids );

	if( space->rankST()==0 ) {
		if( print ) mgSpaces->print();
		for( int i=0; i<mgSpaces->getNGrids(); ++i ) {
			std::string fname = "coord_" + std::to_string( i ) + ".txt";
			Teuchos::RCP<std::ostream> fstream = Pimpact::createOstream( fname, 0 );
			mgSpaces->get(i)->getCoordinatesGlobal()->print( *fstream );
		}
	}

	auto mgPL = Teuchos::parameterList();
	// MG
	mgPL->set<int>( "numCycles", 1 );
	mgPL->set<bool>( "defect correction", false );
	//mgPL->set<bool>( "defect correction", true );
	//mgPL->set<bool>( "init zero", false );
	mgPL->set<bool>( "init zero", true );

	// Smoother: Line
	//mgPL->sublist("Smoother").set<int>( "numIters", 4 );
	//mgPL->sublist("Smoother").set<bool>( "X", true );
	//mgPL->sublist("Smoother").set<bool>( "Y", false );
	//mgPL->sublist("Smoother").set<bool>( "Z", false );

	// Smoother: JT
	//mgPL->sublist("Smoother").set<int>( "BC smoothing", 1 );
	//mgPL->sublist("Smoother").set<OT>( "depth", 2 );


	// Smoother: GS
	//mgPL->sublist("Smoother").set<int>( "RBGS mode", 2 );

	// Smoother: Chebyshev
	// compute EV
	//ST evMax;
	//ST evMin;

	//Teuchos::RCP<Pimpact::TeuchosEigenvalues<Pimpact::DivGradO2Op<typename CS::SpaceT> > > ev = 
	//Teuchos::rcp( new Pimpact::TeuchosEigenvalues<Pimpact::DivGradO2Op<typename CS::SpaceT> >( opO2 ) );
	//ev->computeEV( evMax, evMin );
	//std::cout << "glob: " << evMax << "\t" <<evMin << "\n";
	//////ev->computeFullEV( evMax, evMin );
	//////std::cout << "glob: " << evMax << "\t" <<evMin << "\n";

	//mgPL->sublist("Smoother").set<int>( "numIters", 8 );
	//mgPL->sublist("Smoother").set<ST>( "min EV", evMin*1.1 );
	//mgPL->sublist("Smoother").set<ST>( "max EV", evMin*1.1/30. );

	Teuchos::RCP<MGT> mg =
		Teuchos::rcp( new MGT( mgSpaces, mgPL ) );

	auto prec = Pimpact::createMultiOperatorBase( mg );

	//std::string solvName = "Pseudoblock GMRES";
	//std::string solvName = "Block GMRES";
	//std::string solvName = "Block CG";
	//std::string solvName = "Pseudoblock CG";
	//std::string solvName = "Pseudoblock Stochastic CG";
	//std::string solvName = "GCRODR";
	//std::string solvName = "RCG";
	//std::string solvName = "MINRES";
	//std::string solvName = "LSQR";
	std::string solvName = "TFQMR";
	//std::string solvName = "Pseudoblock TFQMR";
	//std::string solvName = "Hybrid Block GMRES";
	//std::string solvName = "PCPG";
	//std::string solvName = "Fixed Point";
	//std::string solvName = "BiCGStab";

	Teuchos::RCP< Teuchos::ParameterList > param = Teuchos::parameterList();// = Pimpact::createLinSolverParameter( solvName, 1.e-6 );
	param->set( "Output Style", Belos::Brief );
	param->set( "Output Stream", Teuchos::rcp(&std::cout,false) );
	param->set( "Output Frequency", 1 );
	param->set( "Verbosity",			        
			Belos::Errors +
			Belos::Warnings +
			Belos::IterationDetails +
			Belos::OrthoDetails +
			Belos::FinalSummary +
			Belos::TimingDetails +
			Belos::StatusTestDetails +
			Belos::Debug );
	param->set( "Maximum Iterations", nMax );
	//param->set( "Flexible Gmres", false );

	auto bop = Pimpact::createMultiOperatorBase( op );

	auto linprob = Pimpact::createLinearProblem<Pimpact::MultiField<Pimpact::ScalarField<typename CS::SpaceT> > >( bop, xm, bm, param, solvName );
	if( solvName!="Fixed Point" ) 
		linprob->setRightPrec(prec);
	//if( solvName!="BiCGStab" ) 
	if( solvName=="Fixed Point" ) 
		linprob->setLeftPrec(prec);

	// --- zero rhs test ---
	//b->init(0.);
	//x->random();

	//ST error0 = x->norm();
	//ST errorp = error0;

	////opO2->apply( *x, *res );
	//op->apply( *x, *res );
	//ST res0 = res->norm();
	//ST resP = res0;

	//if( space()->rankST()==0 ) {
	//std::cout << "\n\n\t\t\t--- zero rhs test ---\n";
	//std::cout << "\tresidual:\trate:\t\t\terror:\t\trate: \n";
	//std::cout <<  std::scientific;
	//std::cout << "\t" << 1. << "\t\t\t\t"  << 1.  << "\n";
	//}

	//for( int i=0; i<nMax; ++i ) {
	//mg->apply( *b, *x );
	//x->level();

	//ST error = x->norm();

	////opO2->apply( *x, *res );
	//op->apply( *x, *res );
	//ST residual = res->norm();

	//if( space()->rankST()==0 )
	//std::cout << "\t" << residual/res0 << "\t" << residual/resP << "\t\t" << error/error0 << "\t" <<  error/errorp << "\n";

	////if( error>= errorp )
	////break;
	////else
	//errorp = error;
	//resP = residual;
	//}

	////x->print();
	//x->write();
	//TEST_EQUALITY( x->norm()/std::sqrt( static_cast<ST>(x->getLength()) )<1.e-3, true );

	bm->init( 0. );
	xm->random();
	linprob->solve( xm, bm );

	// --- grad test ---
	auto e = x->clone();
	for( int dir=1; dir<=1; ++dir ) {

		Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

		x->init( 0 );
		sol->initField( type );
		sol->level();

		// construct RHS
		opO2->apply( *sol, *b );

		// residual
		opO2->apply( *x, *res );
		res->add( -1, *b, 1., *res );
		ST res0 = res->norm();
		ST resP = res0;

		// error
		res->add( 1., *sol, -1., *x );
		ST err0 = res->norm();
		ST errP = err0;

		// output
		if( space()->rankST()==0 ) {
			std::cout << "\n\n\t\t\t--- " << Pimpact::toString(type) << " test ---\n";
			std::cout << "\tresidual:\trate:\t\t\terror:\t\trate:\n";
			std::cout <<  std::scientific;
			std::cout << "\t"  << 1.<< "\t\t\t\t" << 1.  << "\n";
		}
		for( int i=0; i<nMax; ++i ) {
			// mg cycle
			mg->apply( *b, *x );
			x->level();

			// residual
			opO2->apply( *x, *res );
			res->add( -1, *b, 1., *res );
			ST residual = res->norm();
			res->add( 1., *sol, -1., *x );
			ST err = res->norm();

			if( space()->rankST()==0 )
				std::cout << "\t" << residual/res0 << "\t" <<  residual/resP << "\t\t" << err/err0 << "\t" <<  err/errP  << "\n";
			resP = residual;
			errP = err;
		}
		TEST_EQUALITY( res->norm()/std::sqrt( static_cast<ST>( res->getLength() ) )<1.e-3, true );

		op->apply( *sol, *b );
		xm->init( 0. );
		linprob->solve( xm, bm );
		if( write ) xm->write();
	}


}


//TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiGrid, DivGradOp, CSG2D, DGJMGT2D )
//TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiGrid, DivGradOp, CSG2D, DGLMGT2D )
//TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiGrid, DivGradOp, CSG2D, DGCMGT2D )

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiGrid, DivGradOp, CSG3D, DGJMGT3D )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiGrid, DivGradOp, CSG3D, DGLMGT3D )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiGrid, DivGradOp, CSG3D, DGCMGT3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, ConvDiffSOR, CS ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids );

	auto mgPL = Teuchos::parameterList();
	mgPL->set<int>("numCycles", 1);
	mgPL->sublist("Smoother").set<ST>( "omega", 1. );
	mgPL->sublist("Smoother").set( "numIters", 2 );
	mgPL->sublist("Smoother").set( "Ordering", 1 );
	mgPL->sublist("Smoother").set<short int>( "dir X", -1 );
	mgPL->sublist("Smoother").set<short int>( "dir Y", -1 );
	mgPL->sublist("Smoother").set<short int>( "dir Z", -1 );

	mgPL->sublist("Coarse Grid Solver").set<std::string>("Solver name", "GMRES" );

	//	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<Teuchos::RCP<std::ostream> >( "Output Stream", Teuchos::rcp( &std::cout, false ) );
	//	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set("Verbosity",
	//			Belos::Errors + Belos::Warnings + Belos::IterationDetails +
	//			Belos::OrthoDetails + Belos::FinalSummary + Belos::TimingDetails +
	//			Belos::StatusTestDetails + Belos::Debug );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<std::string>("Timer Label", "Coarse Grid Solver" );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<ST>("Convergence Tolerance"
			, 0.1 );

	auto mg =
		Pimpact::createMultiGrid<
		Pimpact::VectorField,
		TransVF,
		RestrVF,
		InterVF,
		ConvDiffOpT,
		ConvDiffOpT,
		ConvDiffSORT,
		MOP
			>( mgSpaces, mgPL );

	auto op = Pimpact::create< ConvDiffOpT >( space );

	auto x = Pimpact::create<Pimpact::VectorField>( space );
	auto b = Pimpact::create<Pimpact::VectorField>( space );
	auto temp = x->clone();

	{
		auto wind = x->clone();
		wind->initField();
		//		wind->initField( Pimpact::ConstFlow, 1., 1., 1. );
		op->assignField( *wind );
		mg->assignField( *wind );
	}

	std::ofstream ofs;
	if( space()->rankST()==0 )
		ofs.open("MG2.txt", std::ofstream::out);

	// 
	(*x)(Pimpact::F::U).initField( Pimpact::Grad2D_inX );
	(*x)(Pimpact::F::V).initField( Pimpact::Grad2D_inY );
	(*x)(Pimpact::F::W).initField( Pimpact::Grad2D_inZ );
	auto sol = x->clone( Pimpact::ECopy::Deep );

	op->apply(*x,*b);
	{
		x->init(0);
		auto bc = x->clone( Pimpact::ECopy::Shallow );
		op->apply( *x, *bc );
		b->add( 1., *b, -1., *bc );
	}
	if( write ) b->write(1);

	x->initField();
	//x->random();

	temp->add( -1, *x, 1., *sol );
	ST res = temp->norm();
	ST res_0 = res;
	ST res_p = res;

	if( space()->rankST()==0 ) {
		std::cout << "\n\n\tresidual:\trate: \n";
		std::cout << "\t"  << 1.  << "\n";
		ofs << res << "\n";
	}

	for( int i=0; i<20; ++i ) {
		mg->apply( *b, *x );
		if( write ) x->write(i+10);

		temp->add( -1, *x, 1., *sol );
		ST res = temp->norm();

		if( space()->rankST()==0 ) {
			std::cout << "\t" << res/res_0 << "\t" <<  res/res_p << "\n";
			ofs << res << "\n";
		}
		res_p = res;
	}

	if( space()->rankST()==0 ) {
		std::cout << "\n";
	}

	TEST_EQUALITY( temp->norm()<0.5, true );

	if( write ) x->write(2);

	if( space()->rankST()==0 )
		ofs.close();

}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, ConvDiffSOR, CSG2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, ConvDiffSOR, CSG3D )



//TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiGrid, ConvDiffJ, CS, SmootherT ) {
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiGrid, ConvDiffJ, CS ) {

	setParameter( CS::SpaceT::sdim );

	auto space = Pimpact::create<typename CS::SpaceT>( pl );

	auto mgSpaces = Pimpact::createMGSpaces<CS>( space, maxGrids );

	auto mgPL = Teuchos::parameterList();
	mgPL->set<int>("numCycles", 1 );

	mgPL->sublist("Coarse Grid Solver").set<std::string>("Solver name", "GMRES" );

	//	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<Teuchos::RCP<std::ostream> >( "Output Stream", Teuchos::rcp( &std::cout, false ) );
	//	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set("Verbosity", Belos::Errors + Belos::Warnings +
	//			Belos::IterationDetails + Belos::OrthoDetails +
	//			Belos::FinalSummary + Belos::TimingDetails +
	//			Belos::StatusTestDetails + Belos::Debug );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<std::string>("Timer Label", "Coarse Grid Solver" );
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<ST>("Convergence Tolerance" , 1.e-12 );

	auto mg =
		Pimpact::createMultiGrid<
		Pimpact::VectorField,
		TransVF,
		RestrVF,
		InterVF,
		ConvDiffOpT,
		ConvDiffOpT,
		ConvDiffJT,
		ConvDiffJT
		//		ConvDiffSORT,
		//MOP
			>( mgSpaces, mgPL );

	if( print ) mg->print();

	auto op = Pimpact::create< ConvDiffOpT >( space );

	Pimpact::VectorField<typename CS::SpaceT> x( space );
	Pimpact::VectorField<typename CS::SpaceT> b( space );
	Pimpact::VectorField<typename CS::SpaceT> temp( space );

	{
		Pimpact::VectorField<typename CS::SpaceT> wind( space );
		//wind.initField();
		////		wind->initField();
		////		wind->random();
		//wind(Pimpact::F::U).initField( Pimpact::Poiseuille2D_inX );
		//wind(Pimpact::F::V).random();
		//wind(Pimpact::F::V).scale(0.1);
		//wind(Pimpact::F::W).random();
		//wind(Pimpact::F::W).scale(0.1);
		wind.init( 1. );
		op->assignField( wind );
		mg->assignField( wind );
	}

	std::ofstream ofs;
	if( space()->rankST()==0 )
		ofs.open("MG2.txt", std::ofstream::out);

	// 
	//	x(Pimpact::F::V).initField( Pimpact::Poiseuille2D_inY );
	//	x(Pimpact::F::U).initField( Pimpact::Grad2D_inX );
	//	x(Pimpact::F::V).initField( Pimpact::Grad2D_inY );
	x(Pimpact::F::W).initField( Pimpact::Grad2D_inY );
	auto sol = x.clone( Pimpact::ECopy::Deep );

	//	op->apply( x, b );
	//	{
	//		x.init(0);
	//		auto bc = x.clone( Pimpact::ECopy::Shallow );
	//		op->apply( x, *bc );
	//		b.add( 1., b, -1., *bc );
	//	}
	//	b.write(1);

	x.initField();
	sol->initField();
	x.random();

	temp.add( -1, x, 1., *sol );
	ST res = temp.norm();
	ST res_0 = res;
	ST res_p = res;

	if( space()->rankST()==0 ) {
		std::cout <<  std::scientific;
		std::cout << "\n\n\tresidual:\trate: \n";
		std::cout << "\t"  << 1.  << "\n";
		ofs << res << "\n";
	}

	for( int i=0; i<20; ++i ) {
		mg->apply( b, x );
		if( write ) x.write(i+10);

		temp.add( -1, x, 1., *sol );
		//ST res = temp.norm(Belos::TwoNorm, Pimpact::B::N );
		ST res = temp.norm();

		if( space()->rankST()==0 ) {
			std::cout << "\t" << res/res_0 << "\t" <<  res/res_p << "\n";
			ofs << res << "\n";
		}
		res_p = res;
	}

	if( space()->rankST()==0 ) {
		std::cout << "\n";
	}

	TEST_EQUALITY( temp.norm()<0.5, true );

	if( write ) x.write(2);

	if( space()->rankST()==0 )
		ofs.close();
}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, ConvDiffJ, CSG2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiGrid, ConvDiffJ, CSG3D )


} // end of namespace
