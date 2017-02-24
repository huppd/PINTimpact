#include <iostream>

#include "mpi.h"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "BelosOutputManager.hpp"
#include "BelosMVOPTester.hpp"

#include "BelosPimpactAdapter.hpp"

#include "Pimpact_Fields.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_TimeOpWrap.hpp"
#include "Pimpact_DtTimeOp.hpp"
#include "Pimpact_TimeNSOp.hpp"

#include "Pimpact_TimeStokesBSmoother.hpp"

#include "Pimpact_Test.hpp"




namespace {



const int sd = 3;

using SpaceT = Pimpact::Space<ST,OT,sd,d,dNC>;

using SF = Pimpact::ScalarField<SpaceT>;
using VF = Pimpact::VectorField<SpaceT>;

using TSF = Pimpact::TimeField<SF>;
using TVF = Pimpact::TimeField<VF>;

using MTSF = Pimpact::MultiField<TSF>;
using MTVF = Pimpact::MultiField<TVF>;

using SOpBase = Pimpact::OperatorBase<MTSF>;
using VOpBase = Pimpact::OperatorBase<MTVF>;




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, create_init_print, FType ) {

	setParameter( SpaceT::sdim );

	auto space = Pimpact::create<SpaceT>( pl );

	space->print();

	FType p(space);

	p.init( space->rankST() );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, create_init_print, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, create_init_print, TVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, InfNorm_and_init, FType ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

  FType p( space );

  double norm;

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 10.1; ++i ) {
    p.init(i/2.);
    norm = p.norm(Belos::InfNorm);
    TEST_FLOATING_EQUALITY( i/2., norm, eps );

  }

  // one test with infty-norm
  int rank;
  int size;
  double init;
  MPI_Comm_rank(space->comm(),&rank);
  MPI_Comm_size(space->comm(),&size);
  for( double i = 0.; i<10.1; ++i) {
    init = (size-1)*i-1.;
    init = (init<0)?-init:init;
    p.init(rank*i-1.);
    norm = p.norm(Belos::InfNorm);
    TEST_FLOATING_EQUALITY( init, norm, eps );
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNorm_and_init, TSF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNorm_and_init, TVF ) 



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, OneNorm_and_init, FType ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 10.1; ++i ) {
    p->init(i/2.);
//    TEST_EQUALITY( (i/2.)*p->getLength(), p->norm(Belos::OneNorm) );
    TEST_FLOATING_EQUALITY( (i/2.)*p->getLength(), p->norm(Belos::OneNorm), eps );

  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNorm_and_init, TSF ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNorm_and_init, TVF ) 



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, TwoNorm_and_init, FType ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  // test different float values, assures that initial and norm work smoothly
  for( double i=0.; i< 10.1; ++i ) {
    p->init(i/2.);
    TEST_FLOATING_EQUALITY( std::sqrt( std::pow(i/2.,2)*p->getLength() ), p->norm(Belos::TwoNorm), eps );
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNorm_and_init, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNorm_and_init, TVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, dot, FType ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

  auto vel1 = Pimpact::create<FType>(space);
  auto vel2 = Pimpact::create<FType>(space);

  int Np = vel1->getLength();
  int Nq = vel2->getLength();
  double dot;

  TEST_EQUALITY( Np , Nq );
  int N = Np;

  vel1->init(0.);
  vel2->init(1.);
  dot = vel1->dot(*vel2);
  TEST_EQUALITY( dot<eps, true );

  vel1->init(1.);
  vel2->init(1.);
  dot = vel2->dot(*vel1);
  TEST_FLOATING_EQUALITY( static_cast<ST>(N), dot, eps );

  vel1->init(2.);
  vel2->init(1.);
  dot = vel1->dot(*vel2);
  TEST_FLOATING_EQUALITY( 2.*N, dot, eps );

  vel1->init(1.);
  vel2->init(2.);
  dot = vel1->dot(*vel2);
  TEST_FLOATING_EQUALITY( 2.*N, dot, eps );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, TVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, scale, FType ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  double norm;
  int N = p->getLength();

  p->init(1.);
  p->scale(2.);
  norm = p->norm(Belos::TwoNorm);
  TEST_FLOATING_EQUALITY( std::sqrt(4*N), norm, eps )

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, TVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, random, FType ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  double norm;
  int N = p->getLength();

  p->init(1.);
  p->random();
  norm = p->norm(Belos::TwoNorm);
  TEST_INEQUALITY( N, norm)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, TVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TemplateField, add, FType ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

  auto vel1 = Pimpact::create<FType>(space);
  auto vel2 = Pimpact::create<FType>(space);
  auto vel3 = Pimpact::create<FType>(space);

  TEST_EQUALITY( vel1->getLength(), vel2->getLength() )
  TEST_EQUALITY( vel2->getLength(), vel3->getLength() )
  TEST_EQUALITY( vel1->getLength(), vel3->getLength() )

  double norm;
  int N = vel1->getLength();

  vel1->init(0.);
  vel2->init(1./2.);
  vel3->init(1./3.);

  vel1->add( 2., *vel2, 0., *vel3);
  norm = vel1->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(N), norm )

  vel1->init(0.);
  vel2->init(1./2.);
  vel3->init(1./3.);

  vel1->add( 0., *vel2, 3., *vel3);
  norm = vel1->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(N), norm )

  vel1->init(0.);
  vel2->init(1.);
  vel3->init(1.);

  vel1->add( 0.5, *vel2, 0.5, *vel3);
  norm = vel1->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(N), norm )
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, TVF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, write, FType ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

  auto p = Pimpact::create<FType>( space );

  p->init(1.);
  p->write();

  p->random();
  p->write(1);

  TEST_EQUALITY( 0, 0)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, TVF )



// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TimeField, all, FType ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

	space->print();
  //---------------------------------------------------
  auto field1 = Pimpact::createTimeField<FType>( space );

  auto field2 = field1->clone();

  std::cout << "field1: length: " << field1->getLength() << "\n";
  std::cout << "field2: length: " << field2->getLength() << "\n";

  for( ST i=0.; i< 10.1; ++i ) {
    field1->init(i/2.);
    ST norm_ = field1->norm(Belos::InfNorm);
    TEST_EQUALITY( i/2., norm_ );
    field2->init(i/2.);
    norm_ = field2->norm(Belos::InfNorm);
    TEST_EQUALITY( i/2., norm_ );
  }
  field2->init(1.);
  for( double i=0.; i< 10.1; ++i ) {
    field1->init(i/2.);
    TEST_EQUALITY( (i/2.)*field1->getLength(), field1->norm(Belos::OneNorm) );
    TEST_EQUALITY( (i/2.)*field1->getLength(), field1->dot(*field2) );
  }

	// some test
	auto fields = Pimpact::createTimeField<FType>( space );

	auto msca = Pimpact::createMultiField(*fields,5);

	auto field2s = fields->clone();
	field2s->random();
	MPI_Barrier( MPI_COMM_WORLD );
	field2s->norm();
	field2s->init( space->getProcGrid()->getIB(3) );
	field2s->exchange();
	field2s->write();
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TimeField, all, SF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TimeField, all, VF )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TimeField, BelosMVTest, FType ) {

	setParameter( SpaceT::sdim );

  auto space = Pimpact::create<SpaceT>( pl );

  auto vel = Pimpact::create<FType>( space );

  auto mvec = Pimpact::createMultiField( *vel, 5 );

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<ST> > MyOM =
		Teuchos::rcp( new Belos::OutputManager<ST>(Belos::Warnings,rcp(&out,false)) );

  bool res = Belos::TestMultiVecTraits<ST,Pimpact::MultiField<FType> >(MyOM,mvec);
  TEST_EQUALITY_CONST(res,true);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TimeField, BelosMVTest, TSF )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TimeField, BelosMVTest, TVF )



TEUCHOS_UNIT_TEST( TimeOpearotr, TimeOpWrap ) {

	setParameter( SpaceT::sdim );

	auto space = Pimpact::create<SpaceT>( pl );

	auto field = Pimpact::create<TVF>( space );
	auto field1 = field->clone();
	auto field2 = field->clone();

	auto mv = Pimpact::createMultiField( *field, 10 );

	// op test
	auto op = Pimpact::createTimeOpWrap<Pimpact::HelmholtzOp<SpaceT>,true>(
			Pimpact::create<Pimpact::HelmholtzOp>(space ) );

	field2->random();
	Pimpact::initVectorTimeField( field1, Pimpact::Poiseuille_inX );
	op->apply( *field1, *field2 );
	field2->write();
	Pimpact::initVectorTimeField( field2, Pimpact::Zero2DFlow );
	Pimpact::initVectorTimeField( field2, Pimpact::Poiseuille_inX );
	field2->write(10);
	Pimpact::initVectorTimeField( field2, Pimpact::Poiseuille_inY );
	field2->write(20);
	Pimpact::initVectorTimeField( field2, Pimpact::Streaming2DFlow );
	field2->write(30);
	Pimpact::initVectorTimeField( field2, Pimpact::OscilatingDisc2D );
	field2->write(40);

	auto bop =
		Pimpact::createMultiOperatorBase( op );


	Teuchos::RCP<Belos::OutputManager<ST> > myOM = Teuchos::rcp(
			new Belos::OutputManager<ST>(Belos::Warnings+Belos::TimingDetails,rcp(&out,false)) );

	bool res = Belos::TestOperatorTraits< ST, MTVF, VOpBase >( myOM, mv, bop );

  TEST_EQUALITY( res, true );
}



TEUCHOS_UNIT_TEST( TimeOperator, DtTimeOp ) {

	setParameter( SpaceT::sdim );
	// processor grid size
    
	double pi = 4.*std::atan(1.);

	int l = 4;
	std::vector<double> error(l);
	std::vector<double> idt(l);

	for( int q = 0; q < l; q++ ) {

		pl->set<OT>( "nf", 8*std::pow(2,q) );

		OT nt = 8*std::pow(2,q);

		auto space = Pimpact::create<SpaceT>( pl );
		//if( 1==space->rankST() ) {
			//space->print();
			//space->getCoordinatesLocal()->print();
		//}

		ST a2 = space->getDomainSize()->getAlpha2()/space->getDomainSize()->getRe();

		TVF x( space );
		TVF y( space );
		TVF sol( space );
		TVF er( space );

		// op test
		auto dt = Pimpact::create<Pimpact::DtTimeOp>( space );

		// zero test
		x.init(1.);
		y.init(1.);

		dt->apply( x, y );

		TEST_EQUALITY( y.norm()<eps, true );

		// cos test
		for( OT i=space->begin(Pimpact::F::S,3); i<=space->end(Pimpact::F::S,3); ++i ) {
			ST alpha = std::sin( space->getCoordinatesLocal()->getX( Pimpact::F::S, 3, i ) );
			ST beta = a2*std::cos( space->getCoordinatesLocal()->getX( Pimpact::F::S, 3, i ) );
			//std::cout << "\ti: " << i << "\tt: " << space->getCoordinatesLocal()->getX( Pimpact::F::S, 3, i ) << 
				//"\talpha: " << alpha << "\tbeta: " << beta <<"\n";
			x.getField(i).init( alpha );
			sol.getField(i).init( beta );
			x.changed();
			sol.changed();
		}

		dt->apply( x, y );

		////y.write(0);
		////sol.write(1000);
		
		//ST offset = space->getShift(3);
		//for( OT i=space->begin(Pimpact::F::S,3); i<space->end(Pimpact::F::S,3); ++i ) {
			////std::cout << "\ti: " << i << "\tt: " << 2.*pi*(i+offset)/nt << "\tt: " << space->getCoordinatesLocal()->getX(  Pimpact::F::S, Pimpact::ECoord::T, i ) << "\t" << x.getField(i).norm( Belos::InfNorm) << "\n";
		//}
		er.add( 1., y, -1., sol, Pimpact::B::N );

		//for( OT i=space->begin(Pimpact::F::S,3); i<space->end(Pimpact::F::S,3); ++i ) {
			//std::cout << "i: " << i << "er: " << er.getField(i).norm(Belos::InfNorm) << "\n";
		//}
		ST bla = er.norm(Belos::InfNorm);

		error[q] = std::log10( bla ); 
		idt[q]   = std::log10( nt );

		if( 0==space->rankST() )
				std::cout << std::pow(10.,idt[q]) << "\t" << bla << "\n";
	}

	double sl = order(idt,error);
	std::cout << "order: " << sl << "\n";

	TEST_EQUALITY( std::abs(sl + 1.) <  0.05, true );

	pl->set("nf", 8 );

}



TEUCHOS_UNIT_TEST( TimeOperator, TimeDtConvectionDiffusionOp ) {

	setParameter( SpaceT::sdim );

	auto space = Pimpact::create<SpaceT>( pl );

	auto y = Pimpact::create<TVF>( space );
	auto x = y->clone();
	auto sol = y->clone();
	auto wind = y->clone();


	auto op = Pimpact::create<Pimpact::TimeDtConvectionDiffusionOp<SpaceT,false> >( space );

	// test against flow dir
	for( OT i=space->begin(Pimpact::F::U,3); i<=space->end(Pimpact::F::U,3); ++i ) {
		wind->getField(i).init( 2. );
		x->getField(i)(Pimpact::F::U).initField( Pimpact::Grad2D_inX );
		x->getField(i)(Pimpact::F::V).initField( Pimpact::Grad2D_inY );
		x->getField(i)(Pimpact::F::W).initField( Pimpact::Grad2D_inZ );
		sol->getField(i).init( Teuchos::tuple(2.,2.,2.) );
		wind->changed();
		x->changed();
		sol->changed();
	}


	op->assignField( *wind );
	op->apply( *x, *y );

	y->add( 1., *y, -1., *sol );

	std::cout << "error: " << y->norm() << "\n";
}



TEUCHOS_UNIT_TEST( TimeOperator, TimeStokesOp ) {

	setParameter( SpaceT::sdim );

	using OpT = Pimpact::TimeStokesOp<SpaceT>;

	auto space = Pimpact::create<SpaceT>( pl );


	auto op = Pimpact::create<OpT>( space );

	auto x = Pimpact::create<typename OpT::DomainFieldT>( space );
	auto y = Pimpact::create<typename OpT::RangeFieldT>( space );

	// Const diffusion test in X
	//Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Poiseuille_inX );

	//	x->write(1000);
	y->random();

	op->apply( *x, *y );
	//	y->write(0);

	//initVectorTimeField( x->getVFieldPtr(), Pimpact::Const2DFlow, 8./std::pow(space->getDomainSize()->getSize(Pimpact::Y),2)/re, 0., 0. );
	x->add( 1., *x, -1., *y );
	//	x->write(100);

	std::cout << "error in x: " << x->norm() << "\n";
	TEST_EQUALITY( x->norm()<eps, true );

	// Const diffusion test in Y
	//Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Poiseuille_inY );

	//	x->write(1000);
	y->random();

	op->apply( *x, *y );
	//	y->write(0);

	//Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Const2DFlow, 0., 8./std::pow(space->getDomainSize()->getSize(Pimpact::X),2)/re, 0. );
	x->add( 1., *x, -1., *y );
	//	x->write(100);

	std::cout << "error in Y: " << x->norm() << "\n";
	TEST_EQUALITY( x->norm()<eps, true );

	// Gradient test in X
	//Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Zero2DFlow );
	for( OT i=space->begin(Pimpact::F::U,3); i<space->end(Pimpact::F::U,3); ++i ) {
		x->getSFieldPtr()->getField(i).initField( Pimpact::Grad2D_inX );
	}
//	x->write();
	y->random();
	
	op->apply( *x, *y );
//	y->write(100);
	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Const2DFlow, 1., 0., 0. );
	x->getSFieldPtr()->init( 0. );
	x->add( 1., *x, -1., *y );
//	x->write(200);

	std::cout << "error grad in X: " << x->norm() << "\n";
	TEST_EQUALITY( x->norm()<eps, true );

	// Gradient test in Y
	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Zero2DFlow );
	for( OT i=space->begin(Pimpact::F::U,3); i<space->end(Pimpact::F::U,3); ++i ) {
		x->getSFieldPtr()->getField(i).initField( Pimpact::Grad2D_inY );
	}
//	x->write();
	y->random();
	
	op->apply( *x, *y );
//	y->write(100);
	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Const2DFlow, 0., 1., 0. );
	x->getSFieldPtr()->init( 0. );
	x->add( 1., *x, -1., *y );
//	x->write(200);

	std::cout << "error grad in Y: " << x->norm() << "\n";
	TEST_EQUALITY( x->norm()<eps, true );

	// Gradient test in Z
	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Zero2DFlow );
	for( OT i=space->begin(Pimpact::F::U,3); i<space->end(Pimpact::F::U,3); ++i ) {
		x->getSFieldPtr()->getField(i).initField( Pimpact::Grad2D_inZ );
	}
//	x->write();
	y->random();
	
	op->apply( *x, *y );
//	y->write(100);
	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Const2DFlow, 0., 0., 1. );
	x->getSFieldPtr()->init( 0. );
	x->add( 1., *x, -1., *y );
//	x->write(200);

	std::cout << "error grad in Z: " << x->norm() << "\n";
	TEST_EQUALITY( x->norm()<eps, true );

}



TEUCHOS_UNIT_TEST( TimeOperator, TimeStokesOpDT ) {

	using OpT = Pimpact::TimeStokesOp<SpaceT>;

	setParameter( SpaceT::sdim );


	double pi = 4.*std::atan(1.);

	int l = 7;
	std::vector<double> error(l);
	std::vector<double> idt(l);

	for (int q = 3; q < 3+l; q++ ) {

		pl->set<OT>("nf", int(std::pow(2,q)) );

		auto space = Pimpact::create<SpaceT>( pl );


		auto field = Pimpact::create<typename OpT::DomainFieldT>( space );
		auto field1 = Pimpact::create<typename OpT::RangeFieldT>( space );
		auto field2 = Pimpact::create<typename OpT::RangeFieldT>( space );


		// op test
		auto dt = Pimpact::create<OpT>( space );

		// zero test
		Pimpact::initVectorTimeField( field1->getVFieldPtr(), Pimpact::Const2DFlow, 1., 1., 1. );
		field2->random();

		dt->apply( *field1, *field2 );

		TEST_EQUALITY( field2()->norm()<eps, true );

		// cos test
		Pimpact::initVectorTimeField( field->getVFieldPtr(), Pimpact::OscilatingDisc2DVel,0.,0.,0.,1. );
		//field->write();
		field1->init(0);
		field2->init(0);

		dt->apply( *field, *field1 );
		dt->apply( *field1, *field2 );

		ST a2 = space->getDomainSize()->getAlpha2()/space->getDomainSize()->getRe();

		ST bla = a2*a2;

		field->add( bla, *field, 1, *field2 );

		std::cout << "a^4: " << bla << "\n";
		std::cout << "error: " << field()->norm() << "\n";

		error[q-3] = std::log(field()->norm()); 
		idt[q-3]   = std::log(((double)space->nGlo()[3])/2./pi);

	}

	double sl = std::abs(order(idt,error));
	std::cout << "order: " << sl << "\n";

	TEST_EQUALITY( std::abs(sl - 0.5) <  0.05, true );

	pl->set<OT>("nf", 8 );


}


TEUCHOS_UNIT_TEST( TimeOperator, TimeNSOp ) {

	setParameter( SpaceT::sdim );

	using OpT = Pimpact::TimeNSOp<SpaceT>;

	auto space = Pimpact::create<SpaceT>( pl );


	auto op = Pimpact::create<OpT>( space );

	auto x = Pimpact::create<typename OpT::DomainFieldT>( space );
	auto y = Pimpact::create<typename OpT::RangeFieldT>( space );

	// Const diffusion test in X
	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Poiseuille_inX );

	//	x->write(1000);
	y->random();

	op->apply( *x, *y );
	//	y->write(0);

	initVectorTimeField( x->getVFieldPtr(), Pimpact::Const2DFlow, 8./std::pow(space->getDomainSize()->getSize(Pimpact::Y),2), 0., 0. );
	x->add( 1., *x, -1., *y );
	//	x->write(100);

	std::cout << "error in x: " << x->norm() << "\n";
	TEST_EQUALITY( x->norm()<eps, true );

	// Const diffusion test in Y
	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Poiseuille_inY );

	//	x->write(1000);
	y->random();

	op->apply( *x, *y );
	//	y->write(0);

	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Const2DFlow, 0., 8./std::pow(space->getDomainSize()->getSize(Pimpact::X),2)/re, 0. );
	x->add( 1., *x, -1., *y );
	//	x->write(100);

	std::cout << "error in Y: " << x->norm() << "\n";
	TEST_EQUALITY( x->norm()<eps, true );

	// Gradient test in X
	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Zero2DFlow );
	for( OT i=space->begin(Pimpact::F::U,3); i<space->end(Pimpact::F::U,3); ++i ) {
		x->getSFieldPtr()->getField(i).initField( Pimpact::Grad2D_inX );
	}
//	x->write();
	y->random();
	
	op->apply( *x, *y );
//	y->write(100);
	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Const2DFlow, 1., 0., 0. );
	x->getSFieldPtr()->init( 0. );
	x->add( 1., *x, -1., *y );
//	x->write(200);

	std::cout << "error grad in X: " << x->norm() << "\n";
	TEST_EQUALITY( x->norm()<eps, true );

	// Gradient test in Y
	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Zero2DFlow );
	for( OT i=space->begin(Pimpact::F::U,3); i<space->end(Pimpact::F::U,3); ++i ) {
		x->getSFieldPtr()->getField(i).initField( Pimpact::Grad2D_inY );
	}
//	x->write();
	y->random();
	
	op->apply( *x, *y );
//	y->write(100);
	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Const2DFlow, 0., 1., 0. );
	x->getSFieldPtr()->init( 0. );
	x->add( 1., *x, -1., *y );
//	x->write(200);

	std::cout << "error grad in Y: " << x->norm() << "\n";
	TEST_EQUALITY( x->norm()<eps, true );

	// Gradient test in Z
	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Zero2DFlow );
	for( OT i=space->begin(Pimpact::F::U,3); i<space->end(Pimpact::F::U,3); ++i ) {
		x->getSFieldPtr()->getField(i).initField( Pimpact::Grad2D_inZ );
	}
//	x->write();
	y->random();
	
	op->apply( *x, *y );
//	y->write(100);
	Pimpact::initVectorTimeField( x->getVFieldPtr(), Pimpact::Const2DFlow, 0., 0., 1. );
	x->getSFieldPtr()->init( 0. );
	x->add( 1., *x, -1., *y );
//	x->write(200);

	std::cout << "error grad in Z: " << x->norm() << "\n";
	TEST_EQUALITY( x->norm()<eps, true );


	// test against flow dir
	auto wind   = Pimpact::create<typename OpT::DomainFieldT>( space );
	auto field  = Pimpact::create<typename OpT::DomainFieldT>( space );
	auto field1 = Pimpact::create<typename OpT::DomainFieldT>( space );
	auto field2 = Pimpact::create<typename OpT::DomainFieldT>( space );

	for( OT i=space->begin(Pimpact::F::U,3); i<space->end(Pimpact::F::U,3); ++i ) {
		wind->getVFieldPtr()->getField(i)( Pimpact::F::U ).initField( Pimpact::ConstField, 2. );
		wind->getVFieldPtr()->getField(i)( Pimpact::F::V ).initField( Pimpact::ConstField, 2. );
		wind->getVFieldPtr()->getField(i)( Pimpact::F::W ).initField( Pimpact::ConstField, 2. );
		field1->getVFieldPtr()->getField(i)(Pimpact::F::U).initField( Pimpact::Grad2D_inY );
		field1->getVFieldPtr()->getField(i)(Pimpact::F::V).initField( Pimpact::Grad2D_inZ );
		field1->getVFieldPtr()->getField(i)(Pimpact::F::W).initField( Pimpact::Grad2D_inX );
		field2->getVFieldPtr()->getField(i).init( Teuchos::tuple(2.,2.,2.) );
		wind->getVFieldPtr()->changed();
		field1->getVFieldPtr()->changed();
		field2->getVFieldPtr()->changed();
	}

//	wind->write();
//	field1->write(10);
//	field2->write(20);

 op->assignField( *wind );
 op->apply( *field1, *field );

//	field->write(30);

	field->add( 1., *field, -1., *field2 );

//	field->write(40);

	std::cout << "conv errorV: " << field()->getVFieldPtr()->norm() << "\n";
	std::cout << "conv errorS: " << field()->getSFieldPtr()->norm() << "\n";
}



TEUCHOS_UNIT_TEST( TimeOperator, TimeNSOpDT ) {

	using OpT = Pimpact::TimeNSOp<SpaceT>;

	setParameter( SpaceT::sdim );

	double pi = 4.*std::atan(1.);

	int l = 7;
	std::vector<double> error(l);
	std::vector<double> idt(l);

	for (int q = 3; q < 3+l; q++ ) {

		pl->set("nf", int(std::pow(2,q)) );

		auto space = Pimpact::create<SpaceT>( pl );


		auto field = Pimpact::create<typename OpT::DomainFieldT>( space );
		auto field1 = Pimpact::create<typename OpT::RangeFieldT>( space );
		auto field2 = Pimpact::create<typename OpT::RangeFieldT>( space );


		// op test
		auto dt = Pimpact::create<OpT>( space );

		// zero test
		Pimpact::initVectorTimeField( field1->getVFieldPtr(), Pimpact::Const2DFlow, 1., 1., 1. );
		field2->random();

		dt->apply( *field1, *field2 );

		TEST_EQUALITY( field2()->norm()<eps, true );

		// cos test
		Pimpact::initVectorTimeField( field->getVFieldPtr(), Pimpact::OscilatingDisc2DVel,0.,0.,0.,1. );
		//field->write();
		field1->init(0);
		field2->init(0);

		dt->apply( *field, *field1 );
		dt->apply( *field1, *field2 );

		ST a2 = space->getDomainSize()->getAlpha2()/space->getDomainSize()->getRe();

		ST bla = a2*a2;


		field->add( bla, *field, 1, *field2 );


		std::cout << "a^4: " << bla << "\n";
		std::cout << "error: " << field()->norm() << "\n";

		error[q-3] = std::log(field()->norm()); 
		idt[q-3]   = std::log(((double)space->nGlo()[3])/2./pi);

	}

	double sl = std::abs(order(idt,error));
	std::cout << "order: " << sl << "\n";

	TEST_EQUALITY( std::abs(sl - 0.5) <  0.05, true );

	pl->set<OT>("nf", 8 );

}

} // end of namespace
