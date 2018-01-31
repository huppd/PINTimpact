#include <functional>
#include <iostream>
#include <cmath>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Fields.hpp"
#include "Pimpact_LinSolverParameter.hpp"
#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"
#include "Pimpact_OperatorFactory.hpp"
#include "Pimpact_TransferOp.hpp"
#include "Pimpact_SimpleVectorIteration.hpp"
#include "Pimpact_VectorFieldOpWrap.hpp"



template<typename ScalarT>
ScalarT order(const std::vector<ScalarT>& x, const std::vector<ScalarT>& y) {

	TEUCHOS_TEST_FOR_EXCEPT(x.size()!=y.size());

	const ScalarT n    = x.size();
	if(n<2) return(0.);

	const ScalarT s_x  = std::accumulate(x.begin(), x.end(), 0.0);
	const ScalarT s_y  = std::accumulate(y.begin(), y.end(), 0.0);

	const ScalarT s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
	const ScalarT s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);

	const ScalarT a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);

	return(a);
}


namespace {


using ST = double;
using OT= int;
const int sd = 3;
const int d = 3;
const int dNC = 4;
//const int dNC = 2;

using SpaceT = Pimpact::Space<ST, OT, sd, d, dNC>;

using SF = typename Pimpact::ScalarField<SpaceT>;
using VF = typename Pimpact::VectorField<SpaceT>;
using MSF = typename Pimpact::ModeField<SF>;
using MVF = typename Pimpact::ModeField<VF>;



bool testMpi = true;
ST eps = 1e-8;

int domain = 0;

ST lx = 1.;
ST ly = 1.;
ST lz = 1.;

ST omega = 0.8;
ST winds = 1;
int sweeps = 12;
int nIter = 1;

OT nx = 33;
OT ny = 33;
OT nz = 33;
OT nf = 1;

int sx = 0;
int sy = 0;
int sz = 0;

int npx = 1;
int npy = 1;
int npz = 1;
int npf = 1;

Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

//Teuchos::RCP<std::ostream> outp;
//std::ostream& out;

int rank=0;


TEUCHOS_STATIC_SETUP() {

  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build, "
      " this option is ignored and a serial comm is always used.");
  clp.setOption(
      "eps", &eps,
      "Slack off of machine epsilon used to check test results");
	clp.setOption("domain", &domain, "domain");
	clp.setOption("omega", &omega,
      "Slack off of machine epsilon used to check test results");
	clp.setOption("wind", &winds,
      "Slack off of machine epsilon used to check test results");
	clp.setOption("sweeps", &sweeps, "");
	clp.setOption("nIter", &nIter, "");

	clp.setOption("lx", &lx, "");
	clp.setOption("ly", &ly, "");
	clp.setOption("lz", &lz, "");

	clp.setOption("nx", &nx, "");
	clp.setOption("ny", &ny, "");
	clp.setOption("nz", &nz, "");
	clp.setOption("nf", &nf, "");

	clp.setOption("sx", &sx, "");
	clp.setOption("sy", &sy, "");
	clp.setOption("sz", &sz, "");

	clp.setOption("npx", &npx, "");
	clp.setOption("npy", &npy, "");
	clp.setOption("npz", &npz, "");
	clp.setOption("npf", &npf, "");

	//int rank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//if(0==rank)
		//outp =  Teuchos::rcpFromRef(std::cout);
	//else
		//outp = Teuchos::rcp(new Teuchos::oblackholestream());
	//out = *outp;

}


TEUCHOS_UNIT_TEST(BasicOperator, HelmholtzOp) {

	Pimpact::setBoundaryConditions(pl, domain);

	pl->set("lx", lx);
	pl->set("ly", ly);
	pl->set("lz", lz);

	//  grid size
	pl->set("nx", nx);
	pl->set("ny", ny);
	pl->set("nz", nz);
	pl->set("nf", nf);

	// grid stretching
	if(sx!=0) {
		pl->sublist("Stretching in X").set<std::string>("Stretch Type", "cos");
		pl->sublist("Stretching in X").set<ST>("N metr L", static_cast<ST>(nx)/2.);
		pl->sublist("Stretching in X").set<ST>("N metr U", static_cast<ST>(nx)/2.);
	}
	if(sy!=0) {
		pl->sublist("Stretching in Y").set<std::string>("Stretch Type", "cos");
		pl->sublist("Stretching in Y").set<ST>("N metr L", static_cast<ST>(ny)/2.);
		pl->sublist("Stretching in Y").set<ST>("N metr U", static_cast<ST>(ny)/2.);
	}
	if(sz!=0) {
		pl->sublist("Stretching in Z").set<std::string>("Stretch Type", "cos");
		pl->sublist("Stretching in Z").set<ST>("N metr L", static_cast<ST>(nz)/2.);
		pl->sublist("Stretching in Z").set<ST>("N metr U", static_cast<ST>(nz)/2.);
	}

  // processor grid size
  pl->set("npx", npx);
  pl->set("npy", npy);
  pl->set("npz", npz);
  pl->set("npf", npf);

  ST mulI = 5.;
  ST mulL = 3.;

  pl->set<ST>("alpha2", mulI);
  pl->set<ST>("Re", 1./mulL);

  Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>(pl);


  auto op = Pimpact::create<Pimpact::HelmholtzOp>(space);

	Teuchos::RCP<Pimpact::SimpleVectorIteration<Pimpact::HelmholtzOp<SpaceT> > >
		svi = Teuchos::rcp(new
				Pimpact::SimpleVectorIteration<Pimpact::HelmholtzOp<SpaceT> >(op));

	if(0==space->rankST())
		std::cout <<"\n" <<op->getLabel() <<":\tmax EV:\t" <<svi->getMaxEV() <<"\n";
}



TEUCHOS_UNIT_TEST(BasicOperator, DivGradO2Inv) {

	Pimpact::setBoundaryConditions(pl, domain);

	pl->set("lx", lx);
	pl->set("ly", ly);
	pl->set("lz", lz);

	//  grid size
	pl->set("nx", 33);
	pl->set("ny", 17);
	pl->set("nz", 9);
	pl->set("nf", nf);

	// grid stretching
	if(sx!=0) {
		pl->sublist("Stretching in X").set<std::string>("Stretch Type", "cos");
		pl->sublist("Stretching in X").set<ST>("N metr L", 33./2.);
		pl->sublist("Stretching in X").set<ST>("N metr U", 33./2.);
	}
	if(sy!=0) {
		pl->sublist("Stretching in Y").set<std::string>("Stretch Type", "cos");
		pl->sublist("Stretching in Y").set<ST>("N metr L", 17./2.);
		pl->sublist("Stretching in Y").set<ST>("N metr U", 17./2.);
	}
	if(sz!=0) {
		pl->sublist("Stretching in Z").set<std::string>("Stretch Type", "cos");
		pl->sublist("Stretching in Z").set<ST>("N metr L", 9./2.);
		pl->sublist("Stretching in Z").set<ST>("N metr U", 9./2.);
	}

  // processor grid size
  pl->set("npx", npx);
  pl->set("npy", npy);
  pl->set("npz", npz);
  pl->set("npf", npf);

	//const int dNC=2;

  Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>(pl);

  auto op = Pimpact::create<Pimpact::DivGradO2Op>(space);

	op->print2Mat();

  auto ppl = Teuchos::parameterList();

  auto solver = Pimpact::create<Pimpact::DivGradO2Inv>(op, ppl);

	{
		auto svi = Teuchos::rcp(new
				Pimpact::SimpleVectorIteration<Pimpact::DivGradO2Op<SpaceT> >(op));

		if(0==space->rankST())
			std::cout <<"\n" <<op->getLabel() <<":\tmax EV:\t" <<svi->getMaxEV() <<"\n";
	}
	{
		auto svi = Teuchos::rcp(new
				Pimpact::SimpleVectorIteration<Pimpact::DivGradO2Inv<Pimpact::DivGradO2Op<SpaceT> > >(solver));

		if(0==space->rankST())
			std::cout <<"\n" <<op->getLabel() <<":\tmax EV:\t" <<svi->getMaxEV() <<"\n";
	}

}


} // namespace
