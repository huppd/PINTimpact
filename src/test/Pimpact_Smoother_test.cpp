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

#include "Pimpact_TimeStokesBSmoother.hpp"
#include "Pimpact_TimeStokesLSmoother.hpp"

#include "Pimpact_TimeNSBSmoother.hpp"
#include "Pimpact_TimeNS4DBSmoother.hpp"


namespace {


using ST = double;
using OT = int;

const int sd = 3;
const int d = 4;
const int dNC = 2;

using SpaceT = Pimpact::Space<ST, OT, sd, d, dNC>;

bool testMpi = true;
bool output = false;
double eps = 3e-1;
int domain = 0;

int npx = 1;
int npy = 1;
int npz = 1;
int npf = 1;

using SF = Pimpact::ScalarField<SpaceT>;
using VF = Pimpact::VectorField<SpaceT>;

using CF = Pimpact::CompoundField<
  Pimpact::TimeField<Pimpact::VectorField<SpaceT> >,
  Pimpact::TimeField<Pimpact::ScalarField<SpaceT> > >;

using TSF = Pimpact::TimeField<SF>;
using TVF = Pimpact::TimeField<VF>;

using MTSF = Pimpact::MultiField<TSF>;
using MTVF = Pimpact::MultiField<TVF>;

using SOpBase = Pimpact::OperatorBase<MTSF>;
using VOpBase = Pimpact::OperatorBase<MTVF>;

auto pl = Teuchos::parameterList();

TEUCHOS_STATIC_SETUP() {
	Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
	clp.addOutputSetupOptions(true);
	clp.setOption(
			"test-mpi", "test-serial", &testMpi,
			"Test MPI (if available) or force test of serial.  In a serial build, "
			" this option is ignored and a serial comm is always used.");
	clp.setOption(
			"output", "noutput", &output,
			"Test MPI (if available) or force test of serial.  In a serial build, "
			" this option is ignored and a serial comm is always used.");
	clp.setOption(
			"error-tol-slack", &eps,
			"Slack off of machine epsilon used to check test results");
	clp.setOption("domain", &domain, "domain");
	clp.setOption("npx", &npx, "");
	clp.setOption("npy", &npy, "");
	clp.setOption("npz", &npz, "");
	clp.setOption("npf", &npf, "");

	pl->set("Re", 1.);
	pl->set("alpha2", 1.);
	Pimpact::setBoundaryConditions(pl, domain);

	pl->set("lx", 2.);
	pl->set("ly", 2.);
	pl->set("lz", 2.);


	pl->set("nx", 17) ;
	pl->set("ny", 17) ;
	pl->set("nz", 17) ;
	pl->set("nf", 16) ;
}



TEUCHOS_UNIT_TEST(TimeOperator, TimeStokesOperator) {

	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	using OpT = Pimpact::TimeStokesOp<SpaceT>;

	auto space = Pimpact::create<SpaceT>(pl);

	auto op = Pimpact::create<OpT>(space);

	auto x = Pimpact::create<CF>(space);

	auto y = x->clone();
	auto y_cc = x->clone();
	auto error = x->clone();
	auto pulsatile = x->clone();

	double p = 1;
	double alpha = std::sqrt(pl->get<double>("alpha2"));

	//RHS
	Pimpact::initVectorTimeField(y->getVField(), Pimpact::ConstVel_inX, p);

	// true solution //pl->get<int> ("Block Size");
	Pimpact::initVectorTimeField(pulsatile->getVField(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha);

	pulsatile->getSField().init(0);

	// consistecy check
	op->apply(*pulsatile, *y_cc);
	error->add(1., *y_cc, -1., *y);

	if(output && error()->norm()/std::sqrt(error->getLength()) > 0.001) {
		error->write();
		//y->write();
		//y_cc->write(100);
	}

	std::cout <<"|| RHS - RHS_true || = " <<error()->norm()/std::sqrt(error->getLength()) <<std::endl;

	TEST_EQUALITY(error()->norm()/std::sqrt(error->getLength())<0.03, true);


	// constistency check fos Bsmoother
	auto pulsatile2 = pulsatile->clone();
	auto pulsatile3 = pulsatile->clone();

	auto bSmoother = Teuchos::rcp(new Pimpact::TimeStokesBSmoother<OpT>(op));


	bSmoother->apply(*y_cc, *pulsatile2);
	bSmoother->apply(*y, *pulsatile3);

	error->add(1., *pulsatile, -1., *pulsatile2);
	std::cout <<"Consistency error = " <<error()->norm()/std::sqrt(error->getLength())  <<"\n";
	TEST_EQUALITY(error()->norm()<eps, true);

	TEST_EQUALITY(y_cc->getSField().norm()<eps, true);
	std::cout <<"RHS pressure norm= " <<y_cc->getSField().norm()  <<"\n";

	if (output && error()->norm()>eps)
		error->write();

	error->add(1., *pulsatile, -1., *pulsatile3);
	std::cout <<"Consistency error wrt analytical RHS = " <<error()->norm()/std::sqrt(error->getLength())  <<"\n";

}



TEUCHOS_UNIT_TEST(TimeOperator, TimeStokesBSmooth) {

	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	using OpT = Pimpact::TimeStokesOp<SpaceT>;
	auto space = Pimpact::create<SpaceT>(pl);

	auto op = Pimpact::create<OpT>(space);

	auto bSmoother = Teuchos::rcp(new Pimpact::TimeStokesBSmoother<OpT>(op));

	// test smoothing properties
	auto x = Pimpact::create<CF>(space);
	auto y = x->clone();
	auto error = x->clone();
	auto true_sol = x->clone();	

	double p = 1;
	double alpha = std::sqrt(pl->get<double>("alpha2"));

	// RHS
	Pimpact::initVectorTimeField(y->getVField(), Pimpact::ConstVel_inX, p);	

	// true solution //pl->get<int> ("Block Size");
	Pimpact::initVectorTimeField(true_sol->getVField(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha);

	// test smoothing
	x->random();
	//x->scale(10);

	// initial error
	error->add(1., *x, -1., *true_sol);
	if(output)
		error->write();
	std::cout <<"err: " <<error->norm()/std::sqrt(error->getLength()) <<"\n";	

	// aprrox. solution
	bSmoother->apply(*y, *x);		

	// final error
	error->add(1., *x, -1., *true_sol);
	if(output)
		error->write(100);
	std::cout <<"err: " <<error->norm()/std::sqrt(error->getLength()) <<"\n";
	// aprrox. solution2
	bSmoother->apply(*y, *x);

	// final error2
	error->add(1., *x, -1., *true_sol);
	if(output)
		error->write(200);	
	std::cout <<"err: " <<error->norm()/std::sqrt(error->getLength()) <<"\n";
}

TEUCHOS_UNIT_TEST(TimeOperator, TimeStokesLSmooth) {

	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	using OpT = Pimpact::TimeStokesOp<SpaceT>;
	auto space = Pimpact::create<SpaceT>(pl);

	auto op = Pimpact::create<OpT>(space);

	auto bSmoother = Teuchos::rcp(new Pimpact::TimeStokesBSmoother<OpT>(op));      
	auto lSmoother = Teuchos::rcp(new Pimpact::TimeStokesLSmoother<OpT>(op));


	auto x = Pimpact::create<CF>(space);

	auto error = x->clone(); 
	x->random();
	auto x2 = x->clone();
	auto y = x->clone();
	y->random();

	bSmoother->apply(*y, *x);
	lSmoother->apply(*y, *x2, 2);                                        

	error->add(1., *x, -1., *x2);
	//error->write();

	std::cout <<"|| error || = " <<error()->norm()/std::sqrt(error->getLength()) <<std::endl;

	TEST_EQUALITY(error()->norm()<eps, true);
}



TEUCHOS_UNIT_TEST(TimeOperator, TimeStokesBSmooth_conv) {

	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	using OpT = Pimpact::TimeStokesOp<SpaceT>;
	auto space = Pimpact::create<SpaceT>(pl);

	auto op = Pimpact::create<OpT>(space);

	auto bSmoother = Teuchos::rcp(new Pimpact::TimeStokesBSmoother<OpT>(op));

	auto x = Pimpact::create<CF>(space);

	auto y = x->clone();
	auto Ax = x->clone();
	auto error = x->clone();
	auto true_sol = x->clone();

	double p = 1;
	double alpha = std::sqrt(pl->get<double>("alpha2"));

	//Pimpact::initVectorTimeField(y->getVField(), Pimpact::ConstVel_inX, p);

	Pimpact::initVectorTimeField(true_sol->getVField(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha);

	op->apply(*true_sol, *y);

	x->random();
	x->scale(10);

	for (int i = 0; i <20; i++){
		//error->add(1., *x, -1., *true_sol);
		//std::cout  <<"error = " <<error->norm()/std::sqrt(error->getLength());	


		//if (output)
		//      error->write(i*100);

		op->apply(*x, *Ax);
		error->add(1., *Ax, -1., *y);
		std::cout  <<error->norm()/std::sqrt(error->getLength()) <<"\n";

		bSmoother->apply(*y, *x);

		x->level();
	}
}



TEUCHOS_UNIT_TEST(TimeOperator, TimeStokesLSmooth_conv) {

	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	using OpT = Pimpact::TimeStokesOp<SpaceT>;
	auto space = Pimpact::create<SpaceT>(pl);

	auto op = Pimpact::create<OpT>(space);

	auto lSmoother = Teuchos::rcp(new Pimpact::TimeStokesLSmoother<OpT>(op));

	auto x = Pimpact::create<CF>(space);
	auto y = x->clone();
	auto error = x->clone();
	auto true_sol = x->clone();

	//double p = 1;
	//double alpha = std::sqrt(pl->get<double>("alpha2"));

	//Pimpact::initVectorTimeField(y->getVField(), Pimpact::ConstVel_inX, p);

	//Pimpact::initVectorTimeField(true_sol->getVField(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha);

	x->random();
	x->scale(10);

	for (int i = 1; i <50; i++){
		error->add(1., *x, -1., *true_sol);
		std::cout  <<error->norm()/std::sqrt(error->getLength()) <<"\n";
		lSmoother->apply(*y, *x, 2);

		if (i%5==0 &&  output)
			error->write(300+i*100);
	}
}



TEUCHOS_UNIT_TEST(TimeOperator, TimeNSBSmooth) {

	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	using OpT = Pimpact::TimeNSOp<SpaceT>;

	auto space = Pimpact::create<SpaceT>(pl);

	auto op = Pimpact::create<OpT>(space);

	auto bSmoother = Teuchos::rcp(new Pimpact::TimeNSBSmoother<OpT>(op));

	auto x = Pimpact::create<CF>(space);
	auto rhs = x->clone();
	auto error = x->clone();
	auto zero_wind = x->clone();

	//x->random();
	Pimpact::initVectorTimeField(x->getVField(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), 1., 1.);
	op->assignField(*x);
	//op->assignField(*zero_wind);

	auto x2 = x->clone();

	op->apply(*x, *rhs);
	bSmoother->apply(*rhs, *x);

	error->add(1., *x, -1., *x2);
	std::cout <<"Consistency error = " <<error()->norm()/std::sqrt(error->getLength())  <<"\n";
	TEST_EQUALITY(error()->norm()<eps, true);

	//TEST_EQUALITY(rhs->getSField().norm()<eps, true);
	//std::cout <<"RHS pressure norm = " <<y_cc->getSField().norm()  <<"\n";

	if (output && error()->norm()>eps)
		error->write();	

}




TEUCHOS_UNIT_TEST(TimeOperator, TimeNSBSmooth_conv) {

	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	using OpT = Pimpact::TimeNSOp<SpaceT>;
	auto space = Pimpact::create<SpaceT>(pl);

	auto op = Pimpact::create<OpT>(space);

	auto bSmoother = Teuchos::rcp(new Pimpact::TimeNSBSmoother<OpT>(op));

	auto x = Pimpact::create<CF>(space);
	auto y = x->clone();
	auto error = x->clone();
	auto true_sol = x->clone();
	auto Ax = x->clone();
	auto wind = x->clone();

	double p = 1;
	double alpha = std::sqrt(pl->get<double>("alpha2"));

	//Pimpact::initVectorTimeField(y->getVField(), Pimpact::ConstVel_inX, p);

	Pimpact::initVectorTimeField(true_sol->getVField(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha);

	//	op->assignField(*true_sol);
	op->assignField(*wind);
	op->apply(*true_sol, *y);

	// print out convergence
	x->random();
	x->scale(10);

	for (int i = 1; i <20; i++){

		op->apply(*x, *Ax);
		error->add(1., *Ax, -1., *y);
		std::cout  <<error->norm()/std::sqrt(error->getLength()) <<"\n";

		if ((i==1 || i%5==0) && output)
			error->write(300+i*100);

		bSmoother->apply(*y, *x);
	}
}



TEUCHOS_UNIT_TEST(TimeOperator, TimeNS4DBSmoothasfd) {

	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	using OpT = Pimpact::TimeNSOp<SpaceT>;

	auto space = Pimpact::create<SpaceT>(pl);

	auto op = Pimpact::create<OpT>(space);

	auto bSmoother = Teuchos::rcp(new Pimpact::TimeNS4DBSmoother<OpT>(op));

	auto x = Pimpact::create<CF>(space);
	auto rhs = x->clone();
	auto error = x->clone();

	//x->random();
	Pimpact::initVectorTimeField(x->getVField(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), 1, 1);
	op->assignField(*x);

	auto x2 = x->clone();

	op->apply(*x, *rhs);
	bSmoother->apply(*rhs, *x);

	error->add(1., *x, -1., *x2);
	std::cout <<"Consistency error = " <<error()->norm()/std::sqrt(error->getLength())  <<"\n";
	TEST_EQUALITY(error()->norm()<eps, true);

	//TEST_EQUALITY(rhs->getSField().norm()<eps, true);
	//std::cout <<"RHS pressure norm = " <<y_cc->getSField().norm()  <<"\n";

	if (output && error()->norm()>eps)
		error->write();	

}



TEUCHOS_UNIT_TEST(TimeOperator, TimeNS4DBSmooth_conv) {

	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	using OpT = Pimpact::TimeNSOp<SpaceT>;
	auto space = Pimpact::create<SpaceT>(pl);

	auto op = Pimpact::create<OpT>(space);

	auto bSmoother = Teuchos::rcp(new Pimpact::TimeNS4DBSmoother<OpT>(op));

	auto x = Pimpact::create<CF>(space);
	auto y = x->clone();
	auto error = x->clone();
	auto true_sol = x->clone();
	auto Ax = x->clone();
	auto zero_wind = x->clone();

	double p = 1;
	double alpha = std::sqrt(pl->get<double>("alpha2"));

	//Pimpact::initVectorTimeField(y->getVField(), Pimpact::ConstVel_inX, p);

	Pimpact::initVectorTimeField(true_sol->getVField(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), p, alpha);

	op->assignField(*true_sol);
	//op->assignField(*zero_wind);
	op->apply(*true_sol, *y);

	// print out convergence
	x->random();
	x->scale(10);

	for (int i = 1; i <20; i++){

		op->apply(*x, *Ax);
		error->add(1., *Ax, -1., *y);
		std::cout  <<error->norm()/std::sqrt(error->getLength()) <<"\n";

		if ((i==1 || i%5==0) && output)
			error->write(300+i*100);

		bSmoother->apply(*y, *x);
	}
}



TEUCHOS_UNIT_TEST(TimeOperator, TimeNS4DBSmooth) {

	pl->set("npx", npx) ;
	pl->set("npy", npy) ;
	pl->set("npz", npz) ;
	pl->set("npf", npf) ;

	using OpT = Pimpact::TimeNSOp<SpaceT>;

	auto space = Pimpact::create<SpaceT>(pl);

	auto op = Pimpact::create<OpT>(space);

	auto bSmoother = Teuchos::rcp(new Pimpact::TimeNS4DBSmoother<OpT>(op));

	auto x = Pimpact::create<CF>(space);
	auto rhs = x->clone();
	auto error = x->clone();

	//x->random();
	Pimpact::initVectorTimeField(x->getVField(), Pimpact::Pulsatile_inX, pl->get<double>("Re"), 1, 1);
	op->assignField(*x);

	auto x2 = x->clone();

	op->apply(*x, *rhs);
	bSmoother->apply(*rhs, *x);

	error->add(1., *x, -1., *x2);
	std::cout <<"Consistency error = " <<error()->norm()/std::sqrt(error->getLength())  <<"\n";
	TEST_EQUALITY(error()->norm()<eps, true);

	//TEST_EQUALITY(rhs->getSField().norm()<eps, true);
	//std::cout <<"RHS pressure norm = " <<y_cc->getSField().norm()  <<"\n";

	if (output && error()->norm()>eps)
		error->write();	

}



} // end of namespace
