/// Pimpact 
/// \author huppd
/// \date 2018


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

using GridT = Pimpact::Grid<ST, OT, sd, d, dNC>;

using SF = Pimpact::ScalarField<GridT>;
using VF = Pimpact::VectorField<GridT>;

using TSF = Pimpact::TimeField<SF>;
using TVF = Pimpact::TimeField<VF>;

using MTSF = Pimpact::MultiField<TSF>;
using MTVF = Pimpact::MultiField<TVF>;

using SOpBase = Pimpact::OperatorBase<MTSF>;
using VOpBase = Pimpact::OperatorBase<MTVF>;




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TempField, create_init_print, FType) {

	setParameter(GridT::sdim);

	auto grid = Pimpact::create<GridT>(pl);

	grid->print();

	FType p(grid);

	p.init(grid->rankST());
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, create_init_print, TSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, create_init_print, TVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TempField, InfNorm_and_init, FType) {

	setParameter(GridT::sdim);

  auto grid = Pimpact::create<GridT>(pl);

  FType p(grid);

  double norm;

  // test different float values, assures that initial and norm work smoothly
  for(double i=0.; i<10.1; ++i) {
    p.init(i/2.);
    norm = p.norm(Pimpact::ENorm::Inf);
    TEST_FLOATING_EQUALITY(i/2., norm, eps);

  }

  // one test with infty-norm
  int rank;
  int size;
  double init;
  MPI_Comm_rank(grid->comm(), &rank);
  MPI_Comm_size(grid->comm(), &size);
  for(double i = 0.; i<10.1; ++i) {
    init = (size-1)*i-1.;
    init = (init<0)?-init:init;
    p.init(rank*i-1.);
    norm = p.norm(Pimpact::ENorm::Inf);
    TEST_FLOATING_EQUALITY(init, norm, eps);
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, InfNorm_and_init, TSF) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, InfNorm_and_init, TVF) 



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TempField, OneNorm_and_init, FType) {

	setParameter(GridT::sdim);

  auto grid = Pimpact::create<GridT>(pl);

  auto p = Pimpact::create<FType>(grid);

  // test different float values, assures that initial and norm work smoothly
  for(double i=0.; i<10.1; ++i) {
    p->init(i/2.);
//    TEST_EQUALITY((i/2.)*p->getLength(), p->norm(Pimpact::ENorm::One));
    TEST_FLOATING_EQUALITY((i/2.)*p->getLength(), p->norm(Pimpact::ENorm::One), eps);

  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, OneNorm_and_init, TSF) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, OneNorm_and_init, TVF) 



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TempField, TwoNorm_and_init, FType) {

	setParameter(GridT::sdim);

  auto grid = Pimpact::create<GridT>(pl);

  auto p = Pimpact::create<FType>(grid);

  // test different float values, assures that initial and norm work smoothly
  for(double i=0.; i<10.1; ++i) {
    p->init(i/2.);
    TEST_FLOATING_EQUALITY(std::sqrt(std::pow(i/2., 2)*p->getLength()), p->norm(Pimpact::ENorm::Two), eps);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, TwoNorm_and_init, TSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, TwoNorm_and_init, TVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TempField, dot, FType) {

	setParameter(GridT::sdim);

  auto grid = Pimpact::create<GridT>(pl);

  auto vel1 = Pimpact::create<FType>(grid);
  auto vel2 = Pimpact::create<FType>(grid);

  int Np = vel1->getLength();
  int Nq = vel2->getLength();
  double dot;

  TEST_EQUALITY(Np , Nq);
  int N = Np;

  vel1->init(0.);
  vel2->init(1.);
  dot = vel1->dot(*vel2);
  TEST_EQUALITY(dot<eps, true);

  vel1->init(1.);
  vel2->init(1.);
  dot = vel2->dot(*vel1);
  TEST_FLOATING_EQUALITY(static_cast<ST>(N), dot, eps);

  vel1->init(2.);
  vel2->init(1.);
  dot = vel1->dot(*vel2);
  TEST_FLOATING_EQUALITY(2.*N, dot, eps);

  vel1->init(1.);
  vel2->init(2.);
  dot = vel1->dot(*vel2);
  TEST_FLOATING_EQUALITY(2.*N, dot, eps);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, dot, TSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, dot, TVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TempField, scale, FType) {

	setParameter(GridT::sdim);

  auto grid = Pimpact::create<GridT>(pl);

  auto p = Pimpact::create<FType>(grid);

  double norm;
  int N = p->getLength();

  p->init(1.);
  p->scale(2.);
  norm = p->norm(Pimpact::ENorm::Two);
  TEST_FLOATING_EQUALITY(std::sqrt(4*N), norm, eps)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, scale, TSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, scale, TVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TempField, random, FType) {

	setParameter(GridT::sdim);

  auto grid = Pimpact::create<GridT>(pl);

  auto p = Pimpact::create<FType>(grid);

  double norm;
  int N = p->getLength();

  p->init(1.);
  p->random();
  norm = p->norm(Pimpact::ENorm::Two);
  TEST_INEQUALITY(N, norm)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, random, TSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, random, TVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TemplateField, add, FType) {

	setParameter(GridT::sdim);

  auto grid = Pimpact::create<GridT>(pl);

  auto vel1 = Pimpact::create<FType>(grid);
  auto vel2 = Pimpact::create<FType>(grid);
  auto vel3 = Pimpact::create<FType>(grid);

  TEST_EQUALITY(vel1->getLength(), vel2->getLength())
  TEST_EQUALITY(vel2->getLength(), vel3->getLength())
  TEST_EQUALITY(vel1->getLength(), vel3->getLength())

  double norm;
  int N = vel1->getLength();

  vel1->init(0.);
  vel2->init(1./2.);
  vel3->init(1./3.);

  vel1->add(2., *vel2, 0., *vel3);
  norm = vel1->norm(Pimpact::ENorm::Two);
  TEST_EQUALITY(std::sqrt(N), norm)

  vel1->init(0.);
  vel2->init(1./2.);
  vel3->init(1./3.);

  vel1->add(0., *vel2, 3., *vel3);
  norm = vel1->norm(Pimpact::ENorm::Two);
  TEST_EQUALITY(std::sqrt(N), norm)

  vel1->init(0.);
  vel2->init(1.);
  vel3->init(1.);

  vel1->add(0.5, *vel2, 0.5, *vel3);
  norm = vel1->norm(Pimpact::ENorm::Two);
  TEST_EQUALITY(std::sqrt(N), norm)
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TemplateField, add, TSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TemplateField, add, TVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TempField, write, FType) {

	setParameter(GridT::sdim);

  auto grid = Pimpact::create<GridT>(pl);

  auto p = Pimpact::create<FType>(grid);

  p->init(1.);
  p->write();

  p->random();
  p->write(1);

  TEST_EQUALITY(0, 0)

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, write, TSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TempField, write, TVF)



// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TimeField, all, FType) {

	setParameter(GridT::sdim);

  auto grid = Pimpact::create<GridT>(pl);

	grid->print();
  //---------------------------------------------------
  auto field1 = Pimpact::createTimeField<FType>(grid);

  auto field2 = field1->clone();

  std::cout <<"field1: length: " <<field1->getLength() <<"\n";
  std::cout <<"field2: length: " <<field2->getLength() <<"\n";

  for(ST i=0.; i<10.1; ++i) {
    field1->init(i/2.);
    ST norm_ = field1->norm(Pimpact::ENorm::Inf);
    TEST_EQUALITY(i/2., norm_);
    field2->init(i/2.);
    norm_ = field2->norm(Pimpact::ENorm::Inf);
    TEST_EQUALITY(i/2., norm_);
  }
  field2->init(1.);
  for(double i=0.; i<10.1; ++i) {
    field1->init(i/2.);
    TEST_EQUALITY((i/2.)*field1->getLength(), field1->norm(Pimpact::ENorm::One));
    TEST_EQUALITY((i/2.)*field1->getLength(), field1->dot(*field2));
  }

	// some test
	auto fields = Pimpact::createTimeField<FType>(grid);

	auto field2s = fields->clone();
	field2s->random();
	MPI_Barrier(MPI_COMM_WORLD);
	field2s->norm();
	field2s->init(grid->getProcGrid()->getIB(3));
	field2s->exchange();
	field2s->write();
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TimeField, all, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TimeField, all, VF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TimeField, BelosMVTest, FType) {

	setParameter(GridT::sdim);

  auto grid = Pimpact::create<GridT>(pl);

  auto mvec = Teuchos::rcp(new Pimpact::MultiField<FType>(grid, 5));

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<ST> > MyOM =
		Teuchos::rcp(new Belos::OutputManager<ST>(Belos::Warnings, rcp(&out, false)));

  bool res = Belos::TestMultiVecTraits<ST, Pimpact::MultiField<FType> >(MyOM, mvec);
  TEST_EQUALITY_CONST(res, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TimeField, BelosMVTest, TSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TimeField, BelosMVTest, TVF)



TEUCHOS_UNIT_TEST(TimeOpearotr, TimeOpWrap) {

	setParameter(GridT::sdim);

	auto grid = Pimpact::create<GridT>(pl);

	auto mv = Teuchos::rcp(new Pimpact::MultiField<TVF>(grid, 10));

	// op test
	auto op = Pimpact::createTimeOpWrap<Pimpact::DiffusionOp<GridT>, true>(
			Pimpact::create<Pimpact::DiffusionOp>(grid));

	auto bop = Pimpact::createMultiOperatorBase(op);

	Teuchos::RCP<Belos::OutputManager<ST> > myOM = Teuchos::rcp(
			new Belos::OutputManager<ST>(Belos::Warnings+Belos::TimingDetails, rcp(&out, false)));

	bool res = Belos::TestOperatorTraits<ST, MTVF, VOpBase >(myOM, mv, bop);

  TEST_EQUALITY(res, true);
}



TEUCHOS_UNIT_TEST(TimeOperator, DtTimeOp) {

	setParameter(GridT::sdim);
	// processor grid size
    
	//double pi = 4.*std::atan(1.);

	std::vector<ST> error2(ns);
	std::vector<ST> errorInf(ns);
	std::vector<ST> dofs(ns);

	for(OT n=0; n<ns; ++n) {

		pl->set<OT>("nf", n0*std::pow(2, n));

		auto grid = Pimpact::create<GridT>(pl);

		ST a2 = grid->getDomainSize()->getAlpha2()/grid->getDomainSize()->getRe();

		TVF x(grid);
		TVF y(grid);
		TVF sol(grid);
		TVF er(grid);

		// op test
		auto dt = Pimpact::create<Pimpact::DtTimeOp>(grid);

		// zero test
		x.init(1.);
		y.init(1.);

		dt->apply(x, y);

		TEST_EQUALITY(y.norm()<eps, true);

		// cos test
		for(OT i=grid->si(Pimpact::F::S, 3); i<=grid->ei(Pimpact::F::S, 3); ++i) {
			ST alpha = std::sin(grid->getCoordinatesLocal()->getX(Pimpact::F::S, 3, i));
			ST beta = a2*std::cos(grid->getCoordinatesLocal()->getX(Pimpact::F::S, 3, i));
			//std::cout <<"\ti: " <<i <<"\tt: " <<grid->getCoordinatesLocal()->getX(Pimpact::F::S, 3, i) <<
				//"\talpha: " <<alpha <<"\tbeta: " <<beta <<"\n";
			x(i).init(alpha);
			sol(i).init(beta);
			x.changed();
			sol.changed();
		}

		dt->apply(x, y);

		er.add(1., y, -1., sol);

		error2[n] = std::log10(er.norm(Pimpact::ENorm::Two) / sol.norm(Pimpact::ENorm::Two));
		errorInf[n] = std::log10(er.norm(Pimpact::ENorm::Inf) / sol.norm(Pimpact::ENorm::Inf));
		dofs[n] = std::log10(n0*std::pow(2., n));

		if(0==grid->rankST())
			std::cout <<std::pow(10., dofs[n]) <<"\t" <<std::pow(10., error2[n]) <<"\t" <<std::pow(10., errorInf[n]) <<"\n";
	}

	// compute order
	ST order2 = order<ST>(dofs, error2);
	if(0==rank)	
		std::cout <<"DtOp: order two norm in: " <<order2 <<"\n";

	ST orderInf = order<ST>(dofs, errorInf);
	if(0==rank)	
		std::cout <<"DtOp: order inf norm: " <<orderInf <<"\n";
	// test
	TEST_EQUALITY(-order2>1., true);
	TEST_EQUALITY(-orderInf>1., true);
	pl->set("nf", 8);
}



TEUCHOS_UNIT_TEST(TimeOperator, TimeDtConvectionDiffusionOp) {

	setParameter(GridT::sdim);

	auto grid = Pimpact::create<GridT>(pl);

	TVF y(grid);
	TVF x(grid);
	TVF sol(grid);
	TVF wind(grid);

	auto op = Pimpact::create<Pimpact::TimeDtConvectionDiffusionOp<GridT, true> >(grid);

	// test against flow dir
	for(OT i=grid->si(Pimpact::F::U, 3); i<=grid->ei(Pimpact::F::U, 3); ++i) {
		wind(i).init(2.);
		x(i)(Pimpact::F::U).initField(Pimpact::Grad2D_inX);
		x(i)(Pimpact::F::V).initField(Pimpact::Grad2D_inY);
		x(i)(Pimpact::F::W).initField(Pimpact::Grad2D_inZ);
		sol(i).init(2.);
		wind.changed();
		x.changed();
		sol.changed();
	}


	op->assignField(wind);
	op->apply(x, y);

	y.add(1., y, -1., sol, Pimpact::B::N);

	if(write) y.write();

	ST error = y.norm(Pimpact::ENorm::Two, Pimpact::B::N);
	std::cout <<"error: " << error <<"\n";
	TEST_EQUALITY(error<eps, true);
}



TEUCHOS_UNIT_TEST(TimeOperator, TimeStokesOp) {

	setParameter(GridT::sdim);

	using OpT = Pimpact::TimeStokesOp<GridT>;

	auto grid = Pimpact::create<GridT>(pl);


	auto op = Pimpact::create<OpT>(grid);

	auto x = Pimpact::create<typename OpT::DomainFieldT>(grid);
	auto y = Pimpact::create<typename OpT::RangeFieldT>(grid);

	// Const diffusion test in X
	//Pimpact::initVectorTimeField(x->getVField(), Pimpact::Poiseuille_inX);

	//	x->write(1000);
	y->random();

	op->apply(*x, *y);
	//	y->write(0);

	//initVectorTimeField(x->getVField(), Pimpact::Const2DFlow, 8./std::pow(grid->getDomainSize()->getSize(Pimpact::Y), 2)/re, 0., 0.);
	x->add(1., *x, -1., *y);
	//	x->write(100);

	std::cout <<"error in x: " <<x->norm() <<"\n";
	TEST_EQUALITY(x->norm()<eps, true);

	// Const diffusion test in Y
	//Pimpact::initVectorTimeField(x->getVField(), Pimpact::Poiseuille_inY);

	//	x->write(1000);
	y->random();

	op->apply(*x, *y);
	//	y->write(0);

	//Pimpact::initVectorTimeField(x->getVField(), Pimpact::Const2DFlow, 0., 8./std::pow(grid->getDomainSize()->getSize(Pimpact::X), 2)/re, 0.);
	x->add(1., *x, -1., *y);
	//	x->write(100);

	std::cout <<"error in Y: " <<x->norm() <<"\n";
	TEST_EQUALITY(x->norm()<eps, true);

	// Gradient test in X
	//Pimpact::initVectorTimeField(x->getVField(), Pimpact::Zero2DFlow);
	for(OT i=grid->si(Pimpact::F::U, 3); i<grid->ei(Pimpact::F::U, 3); ++i) {
		x->getSField()(i).initField(Pimpact::Grad2D_inX);
	}
//	x->write();
	y->random();
	
	op->apply(*x, *y);
//	y->write(100);
	Pimpact::initVectorTimeField(x->getVField(), Pimpact::Const2DFlow, 1., 0., 0.);
	x->getSField().init(0.);
	x->add(1., *x, -1., *y);
//	x->write(200);

	std::cout <<"error grad in X: " <<x->norm() <<"\n";
	TEST_EQUALITY(x->norm()<eps, true);

	// Gradient test in Y
	Pimpact::initVectorTimeField(x->getVField(), Pimpact::Zero2DFlow);
	for(OT i=grid->si(Pimpact::F::U, 3); i<grid->ei(Pimpact::F::U, 3); ++i) {
		x->getSField()(i).initField(Pimpact::Grad2D_inY);
	}
//	x->write();
	y->random();
	
	op->apply(*x, *y);
//	y->write(100);
	Pimpact::initVectorTimeField(x->getVField(), Pimpact::Const2DFlow, 0., 1., 0.);
	x->getSField().init(0.);
	x->add(1., *x, -1., *y);
//	x->write(200);

	std::cout <<"error grad in Y: " <<x->norm() <<"\n";
	TEST_EQUALITY(x->norm()<eps, true);

	// Gradient test in Z
	Pimpact::initVectorTimeField(x->getVField(), Pimpact::Zero2DFlow);
	for(OT i=grid->si(Pimpact::F::U, 3); i<grid->ei(Pimpact::F::U, 3); ++i) {
		x->getSField()(i).initField(Pimpact::Grad2D_inZ);
	}
//	x->write();
	y->random();
	
	op->apply(*x, *y);
//	y->write(100);
	Pimpact::initVectorTimeField(x->getVField(), Pimpact::Const2DFlow, 0., 0., 1.);
	x->getSField().init(0.);
	x->add(1., *x, -1., *y);
//	x->write(200);

	std::cout <<"error grad in Z: " <<x->norm() <<"\n";
	TEST_EQUALITY(x->norm()<eps, true);

}



TEUCHOS_UNIT_TEST(TimeOperator, TimeStokesOpDT) {

	using OpT = Pimpact::TimeStokesOp<GridT>;

	setParameter(GridT::sdim);


	double pi = 4.*std::atan(1.);

	int l = 7;
	std::vector<double> error(l);
	std::vector<double> idt(l);

	for (int q = 3; q <3+l; q++) {

		pl->set<OT>("nf", int(std::pow(2, q)));

		auto grid = Pimpact::create<GridT>(pl);


		auto field = Pimpact::create<typename OpT::DomainFieldT>(grid);
		auto field1 = Pimpact::create<typename OpT::RangeFieldT>(grid);
		auto field2 = Pimpact::create<typename OpT::RangeFieldT>(grid);


		// op test
		auto dt = Pimpact::create<OpT>(grid);

		// zero test
		Pimpact::initVectorTimeField(field1->getVField(), Pimpact::Const2DFlow, 1., 1., 1.);
		field2->random();

		dt->apply(*field1, *field2);

		TEST_EQUALITY(field2()->norm()<eps, true);

		// cos test
		Pimpact::initVectorTimeField(field->getVField(), Pimpact::OscilatingDisc2DVel, 0., 0., 0., 1.);
		//field->write();
		field1->init(0);
		field2->init(0);

		dt->apply(*field, *field1);
		dt->apply(*field1, *field2);

		ST a2 = grid->getDomainSize()->getAlpha2()/grid->getDomainSize()->getRe();

		ST bla = a2*a2;

		field->add(bla, *field, 1, *field2);

		std::cout <<"a^4: " <<bla <<"\n";
		std::cout <<"error: " <<field()->norm() <<"\n";

		error[q-3] = std::log(field()->norm()); 
		idt[q-3]   = std::log(((double)grid->nGlo()[3])/2./pi);

	}

	double sl = std::abs(order(idt, error));
	std::cout <<"order: " <<sl <<"\n";

	TEST_EQUALITY(std::abs(sl - 0.5) < 0.05, true);

	pl->set<OT>("nf", 8);


}


TEUCHOS_UNIT_TEST(TimeOperator, TimeNSOp) {

	setParameter(GridT::sdim);

	using OpT = Pimpact::TimeNSOp<GridT>;

	auto grid = Pimpact::create<GridT>(pl);


	auto op = Pimpact::create<OpT>(grid);

	typename OpT::DomainFieldT e(grid);
	typename OpT::RangeFieldT  y(grid);
	typename OpT::DomainFieldT  sol(grid);
	typename OpT::RangeFieldT  fu(grid);

	ST pi2 = 2.*std::acos(-1.);
	ST A =   0.5;
	ST B =  -0.5;
	ST C =   0.;
	ST a =   1.;
	ST b =   1.;
	ST c =   1.;

	for(OT i=grid->si(Pimpact::F::U, 3); i<=grid->ei(Pimpact::F::U, 3); ++i) {
		{
			// init RHS
			ST timei = grid->getCoordinatesLocal()->getX(Pimpact::F::U, 3, i);
			ST timeim =grid->getCoordinatesLocal()->getX(Pimpact::F::U, 3, i-1);
			//std::cout <<time <<"\n";
			//ST ctime = (std::cos(timei)+std::cos(timeim))/2.;
			ST ctime = std::cos((timei+timeim)/2.);
			ST stime = (std::sin(timei)+std::sin(timeim))/2.;
			fu.getVField()(i)(Pimpact::F::U).initFromFunction(
					[=](ST x, ST y, ST z) ->ST {
					return(
						alpha2/re*A*std::cos(a*pi2*x)*std::sin(b*pi2*y)*(ctime)									// \alpha^2 dt u
						+2.*A/re*std::cos(a*pi2*x)*std::sin(b*pi2*y)*(1.+stime)); });	// -\lap u

			fu.getVField()(i)(Pimpact::F::V).initFromFunction(
					[=](ST x, ST y, ST z) ->ST {
					return(
						alpha2/re*B*std::sin(a*pi2*x)*std::cos(b*pi2*y)*ctime									// \alpha^2 dt v
						+2.*B/re*std::sin(a*pi2*x)*std::cos(b*pi2*y)*(1.+stime)); });	// -\lap u


		}
		// init sol
		{
			ST time = grid->getCoordinatesLocal()->getX(Pimpact::F::U, 3, i);
			//std::cout<<time <<"\n";
			ST stime = std::sin(time);
			sol.getVField()(i)(Pimpact::F::U).initFromFunction(
					[=](ST x, ST y, ST z) ->ST {
					return(A*std::cos(a*pi2*x)*std::sin(b*pi2*y)*(1.+stime)); });
			sol.getVField()(i)(Pimpact::F::V).initFromFunction(
					[=](ST x, ST y, ST z) ->ST {
					return(B*std::sin(a*pi2*x)*std::cos(b*pi2*y)*(1.+stime)); });
		}

	}
	sol.getVField().changed();
	fu.getVField().changed();

	op->apply(sol, y);

	e.add(1., fu, -1., y);

	ST error = e.norm()/sol.norm();
	std::cout <<"error in x: " <<error <<"\n";
	TEST_EQUALITY(error<eps, true);

	//// Const diffusion test in Y
	//Pimpact::initVectorTimeField(x->getVField(), Pimpact::Poiseuille_inY);

	////	x->write(1000);
	//y->random();

	//op->apply(*x, *y);
	////	y->write(0);

	////Pimpact::initVectorTimeField(x->getVField(), Pimpact::Const2DFlow, 0., 8./std::pow(grid->getDomainSize()->getSize(Pimpact::X), 2)/re, 0.);
	//x->add(1., *x, -1., *y);
	////	x->write(100);

	//std::cout <<"error in Y: " <<x->norm() <<"\n";
	//TEST_EQUALITY(x->norm()<eps, true);

	//// Gradient test in X
	//Pimpact::initVectorTimeField(x->getVField(), Pimpact::Zero2DFlow);
	//for(OT i=grid->si(Pimpact::F::U, 3); i<grid->ei(Pimpact::F::U, 3); ++i) {
		//x->getSField()(i).initField(Pimpact::Grad2D_inX);
	//}
////	x->write();
	//y->random();
	
	//op->apply(*x, *y);
////	y->write(100);
	//Pimpact::initVectorTimeField(x->getVField(), Pimpact::Const2DFlow, 1., 0., 0.);
	//x->getSField().init(0.);
	//x->add(1., *x, -1., *y);
////	x->write(200);

	//std::cout <<"error grad in X: " <<x->norm() <<"\n";
	//TEST_EQUALITY(x->norm()<eps, true);

	//// Gradient test in Y
	//Pimpact::initVectorTimeField(x->getVField(), Pimpact::Zero2DFlow);
	//for(OT i=grid->si(Pimpact::F::U, 3); i<grid->ei(Pimpact::F::U, 3); ++i) {
		//x->getSField()(i).initField(Pimpact::Grad2D_inY);
	//}
////	x->write();
	//y->random();
	
	//op->apply(*x, *y);
	////	y->write(100);
	//Pimpact::initVectorTimeField(x->getVField(), Pimpact::Const2DFlow, 0., 1., 0.);
	//x->getSField().init(0.);
	//x->add(1., *x, -1., *y);
	////	x->write(200);

	//std::cout <<"error grad in Y: " <<x->norm() <<"\n";
	//TEST_EQUALITY(x->norm()<eps, true);

	//// Gradient test in Z
	//Pimpact::initVectorTimeField(x->getVField(), Pimpact::Zero2DFlow);
	//for(OT i=grid->si(Pimpact::F::U, 3); i<grid->ei(Pimpact::F::U, 3); ++i) {
		//x->getSField()(i).initField(Pimpact::Grad2D_inZ);
	//}
	////	x->write();
	//y->random();
	
	//op->apply(*x, *y);
////	y->write(100);
	//Pimpact::initVectorTimeField(x->getVField(), Pimpact::Const2DFlow, 0., 0., 1.);
	//x->getSField().init(0.);
	//x->add(1., *x, -1., *y);
////	x->write(200);

	//std::cout <<"error grad in Z: " <<x->norm() <<"\n";
	//TEST_EQUALITY(x->norm()<eps, true);


	//// test against flow dir
	//auto wind   = Pimpact::create<typename OpT::DomainFieldT>(grid);
	//auto field  = Pimpact::create<typename OpT::DomainFieldT>(grid);
	//auto field1 = Pimpact::create<typename OpT::DomainFieldT>(grid);
	//auto field2 = Pimpact::create<typename OpT::DomainFieldT>(grid);

	//for(OT i=grid->si(Pimpact::F::U, 3); i<grid->ei(Pimpact::F::U, 3); ++i) {
		//wind->getVField()(i)(Pimpact::F::U).initField(Pimpact::ConstField, 2.);
		//wind->getVField()(i)(Pimpact::F::V).initField(Pimpact::ConstField, 2.);
		//wind->getVField()(i)(Pimpact::F::W).initField(Pimpact::ConstField, 2.);
		//field1->getVField()(i)(Pimpact::F::U).initField(Pimpact::Grad2D_inY);
		//field1->getVField()(i)(Pimpact::F::V).initField(Pimpact::Grad2D_inZ);
		//field1->getVField()(i)(Pimpact::F::W).initField(Pimpact::Grad2D_inX);
		//field2->getVField()(i).init(2.);
		//wind->getVField().changed();
		//field1->getVField().changed();
		//field2->getVField().changed();
	//}

////	wind->write();
////	field1->write(10);
////	field2->write(20);

 //op->assignField(*wind);
 //op->apply(*field1, *field);

////	field->write(30);

	//field->add(1., *field, -1., *field2);

////	field->write(40);

	//std::cout <<"conv errorV: " <<field()->getVField().norm() <<"\n";
	//std::cout <<"conv errorS: " <<field()->getSField().norm() <<"\n";
}



TEUCHOS_UNIT_TEST(TimeOperator, TimeNSOpDT) {

	using OpT = Pimpact::TimeNSOp<GridT>;

	setParameter(GridT::sdim);

	double pi = 4.*std::atan(1.);

	int l = 7;
	std::vector<double> error(l);
	std::vector<double> idt(l);

	for (int q = 3; q <3+l; q++) {

		pl->set("nf", int(std::pow(2, q)));

		auto grid = Pimpact::create<GridT>(pl);


		auto field = Pimpact::create<typename OpT::DomainFieldT>(grid);
		auto field1 = Pimpact::create<typename OpT::RangeFieldT>(grid);
		auto field2 = Pimpact::create<typename OpT::RangeFieldT>(grid);


		// op test
		auto dt = Pimpact::create<OpT>(grid);

		// zero test
		Pimpact::initVectorTimeField(field1->getVField(), Pimpact::Const2DFlow, 1., 1., 1.);
		field2->random();

		dt->apply(*field1, *field2);

		TEST_EQUALITY(field2()->norm()<eps, true);

		// cos test
		Pimpact::initVectorTimeField(field->getVField(), Pimpact::OscilatingDisc2DVel, 0., 0., 0., 1.);
		//field->write();
		field1->init(0);
		field2->init(0);

		dt->apply(*field, *field1);
		dt->apply(*field1, *field2);

		ST a2 = grid->getDomainSize()->getAlpha2()/grid->getDomainSize()->getRe();

		ST bla = a2*a2;


		field->add(bla, *field, 1, *field2);


		std::cout <<"a^4: " <<bla <<"\n";
		std::cout <<"error: " <<field()->norm() <<"\n";

		error[q-3] = std::log(field()->norm()); 
		idt[q-3]   = std::log(((double)grid->nGlo()[3])/2./pi);

	}

	double sl = std::abs(order(idt, error));
	std::cout <<"order: " <<sl <<"\n";

	TEST_EQUALITY(std::abs(sl - 0.5) < 0.05, true);

	pl->set<OT>("nf", 8);

}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TimeOperator, TimeDtConvectionDiffusionOp2, OpT) {

	using SpT = typename OpT::GridT;

	pl->set<bool>("spectral in time", false);

	const ST& pi2 = 2.*std::acos(-1.);

	const ST& a =  1.;
	const ST& A =  1.;
	const ST& b =  1.;
	const ST& B = -1.;

	lx = pi2 ;
	ly = pi2 ;
	lz = pi2 ;

	setParameter(SpT::sdim);

	Teuchos::SerialDenseMatrix<OT, ST> error2(ns, ns);
	Teuchos::SerialDenseMatrix<OT, ST> errorInf(ns, ns);
	std::vector<ST> dofS(ns);
	std::vector<ST> dofT(ns);

	//for(OT nX=0; nX<ns; ++nX) {
	{ OT nX = 2;
		for(OT nT=0; nT<ns; ++nT) {

			pl->set<OT>("nx", n0*std::pow(2, nX)+1);
			pl->set<OT>("ny", n0*std::pow(2, nX)+1);
			pl->set<OT>("nz", 5);
			pl->set<OT>("nf", 8*n0*std::pow(2, nT));
			//pl->set<OT>("nf", 8);
			dofS[nX] = pl->get<OT>("nx");
			dofT[nT] = pl->get<OT>("nf");

			// grid stretching
			setStretching();

			Teuchos::RCP<const SpT> grid = Pimpact::create<SpT>(pl);
			//std::cout <<"re: " <<grid->getDomainSize()->getRe() <<"\n";
			//std::cout <<"re: " <<re <<"\n";
			//std::cout <<"alpha: " <<grid->getDomainSize()->getAlpha2() <<"\n";
			//std::cout <<"alpha: " <<alpha2 <<"\n";

			Pimpact::TimeField<Pimpact::VectorField<SpT> > x(grid);
			Pimpact::TimeField<Pimpact::VectorField<SpT> > y(grid);
			Pimpact::TimeField<Pimpact::VectorField<SpT> > sol(grid);
			Pimpact::TimeField<Pimpact::VectorField<SpT> > err(grid);

			auto op = Pimpact::create<OpT>(grid);

			// initializtion
			//std::cout <<"\n";
			for(OT i=grid->si(Pimpact::F::U, 3); i<=grid->ei(Pimpact::F::U, 3); ++i) {
				ST time = grid->getCoordinatesLocal()->getX(Pimpact::F::U, 3, i);
				ST stime = std::sin(time);
				//std::cout <<i <<"\t" <<time <<"\t" <<stime <<"\n";
				x(i)(Pimpact::F::U).initFromFunction(
						[&](ST x, ST y, ST z) ->ST { return(A*std::cos(a*pi2*x)*std::sin(b*pi2*y)*stime); });
				x(i)(Pimpact::F::V).initFromFunction(
						[&](ST x, ST y, ST z) ->ST { return(B*std::sin(a*pi2*x)*std::cos(b*pi2*y)*stime); });
			}
			x.changed();
			//if(write) x.write();


			// solution init
			for(OT i=grid->si(Pimpact::F::U, 3); i<=grid->ei(Pimpact::F::U, 3); ++i) {
				ST time = (grid->getCoordinatesLocal()->getX(Pimpact::F::U, 3, i-1) + grid->getCoordinatesLocal()->getX(Pimpact::F::U, 3, i))/2./*-grid->getCoordinatesLocal()->getX(Pimpact::F::U, 3, 1)/2.*/;
				ST stime = std::sin(time);
				ST s2time = std::pow(stime, 2);
				ST ctime = std::cos(time);
				sol(i)(Pimpact::F::U).initFromFunction(
						[=](ST x, ST y, ST z) ->ST {
						if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(0)>0) ||
							( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(0)>0) ||
							( (y  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(1)>0) ||
							( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(1)>0) ||
							( (z  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(2)>0) ||
							( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(2)>0))
							return(A*std::cos(a*pi2*std::min(std::max(x, 0.), 1.))*std::sin(b*pi2*y)*stime);
						else
							return(
								alpha2/re*A*std::cos(a*pi2*x)*std::sin(b*pi2*y)*ctime									// \alpha^2 dt u
								-a*A*A/2.*std::sin(2.*a*pi2*x)*s2time 																// (\u * \na) u
								+A*(a*a + b*b)/re*std::cos(a*pi2*x)*std::sin(b*pi2*y)*stime); });	// -\lap u

				sol(i)(Pimpact::F::V).initFromFunction(
						[=](ST x, ST y, ST z) ->ST {
						if(((x  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(0)>0) ||
							( (x-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(0)>0) ||
							( (y  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(1)>0) ||
							( (y-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(1)>0) ||
							( (z  )<= Teuchos::ScalarTraits<ST>::eps() && grid->bcl(2)>0) ||
							( (z-1.)>=-Teuchos::ScalarTraits<ST>::eps() && grid->bcu(2)>0))
							return(B*std::sin(a*pi2*x)*std::cos(b*pi2*std::max(std::min(y, 1.), 0.))*stime);
						else
							return(
								alpha2/re*B*std::sin(a*pi2*x)*std::cos(b*pi2*y)*ctime									// \alpha^2 dt v
								-b*B*B/2.*std::sin(2.*b*pi2*y)*s2time 																// (\u * \na) v
								+B*(a*a + b*b)/re*std::sin(a*pi2*x)*std::cos(b*pi2*y)*stime); });	// -\lap u
			}
			y.changed();

			//if(write) sol.write(4*dofT[nT]);
			//if(write) sol.write();

			//if(print) {
				//std::cout <<"\n--- x ---\n";
				//std::cout <<"\n--- x ---\n";
				//std::cout <<"\n--- x ---\n";
				//x.print();
			//}
			op->assignField(x);
			op->apply(x, y);

			//if(print) {
				//std::cout <<"\n--- y ---\n";
				//std::cout <<"\n--- y ---\n";
				//std::cout <<"\n--- y ---\n";
				//y.print();
			//}
			//if(write) y.write(8*dofT[nT]);
			if(write) y.write();

			err.add(1., sol, -1., y);
			//if(write) err.write();

			//if(print) {
				//std::cout <<"\n--- sol ---\n";
				//sol.print();
			//}
			if(print) {
				std::cout <<"\n--- err ---\n";
				err.print();
			}
			//for(OT i=grid->si(Pimpact::F::U, 3); i<=grid->ei(Pimpact::F::U, 3); ++i) {
				//std::cout <<"err(" <<grid->getCoordinatesLocal()->getX(Pimpact::F::U, 3, i) <<"): " <<err(i).norm() <<"\n";
			//}

			error2(nX, nT)   = err.norm()/sol.norm();
			errorInf(nX, nT) = err.norm(Pimpact::ENorm::Inf)/sol.norm(Pimpact::ENorm::Inf);

			std::cout <<"\nnx: " <<dofS[nX] <<"^2\n";
			std::cout <<"\nnt: " <<dofT[nT] <<"\n";
			//TEST_EQUALITY(error<(1./nx/ny), true);
		}
	}
	std::cout <<"error 2:\n";
	std::cout <<error2 <<"\n";
	std::cout <<"error Inf:\n";
	std::cout <<errorInf <<"\n";
}

//using Op2DFDT = Pimpact::TimeDtConvectionDiffusionOp<D2, 0>;
//using Op3DFDT = Pimpact::TimeDtConvectionDiffusionOp<D3, 0>;
using Op2DCNT = Pimpact::TimeDtConvectionDiffusionOp<D2, 1> ;
using Op3DCNT = Pimpact::TimeDtConvectionDiffusionOp<D3, 1> ;
//using Op2DFD2T = Pimpact::TimeDtConvectionDiffusionOp<D2, 2>;

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TimeOperator, TimeDtConvectionDiffusionOp2, Op2DFDT)
//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TimeOperator, TimeDtConvectionDiffusionOp2, Op2DFD2T)
//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TimeOperator, TimeDtConvectionDiffusionOp2, Op3DFDT)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TimeOperator, TimeDtConvectionDiffusionOp2, Op2DCNT)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TimeOperator, TimeDtConvectionDiffusionOp2, Op3DCNT)


} // end of namespace
