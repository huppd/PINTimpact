/// Pimpact 
/// \author huppd
/// \date 2018


#include <iostream>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Pimpact_Fields.hpp"




namespace {



using ST = double;
using OT = int;

const int sd = 3;
const int d = 3;
const int dNC=4;

using GridT = Pimpact::Grid<ST, OT, sd, d, dNC>;

using SF  = Pimpact::MultiField<Pimpact::ScalarField<GridT> >;
using VF  = Pimpact::MultiField<Pimpact::VectorField<GridT> >;
using MSF = Pimpact::MultiField<Pimpact::ModeField<Pimpact::ScalarField<GridT> > >;
using MVF = Pimpact::MultiField<Pimpact::ModeField<Pimpact::VectorField<GridT> > >;


bool testMpi = true;
double eps = 1e-6;
int domain = 0;

auto pl = Teuchos::parameterList();


TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
	clp.setOption("test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build, "
      " this option is ignored and a serial comm is always used.");
	clp.setOption("error-tol-slack", &eps,
      "Slack off of machine epsilon used to check test results");
	clp.setOption("domain", &domain, "domain");

  pl->set("nx", 33);
  pl->set("ny", 17);
  pl->set("nz", 9);
}



// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, constructor, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv = Teuchos::rcp(new FType(grid, 10));

  const int m = mv->getNumberVecs();

  TEST_EQUALITY(10, m);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, constructor, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, constructor, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, constructor, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, constructor, MVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, TwoNorm, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv = Teuchos::rcp(new FType(grid, 10));

  const OT m = mv->getNumberVecs();
  const OT n = mv->getLength();
  std::vector<ST> normval(m);

  // test different float values, assures that initial and norm work smoothly
  for(ST i=0.; i<200.1; ++i) {
    mv->init(i/2.);
    mv->norm(normval, Pimpact::ENorm::Two);
    for(int j=0; j<m; ++j)
      TEST_FLOATING_EQUALITY(std::sqrt(std::pow(i/2., 2)*n), normval[j], eps);
  }

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, TwoNorm, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, TwoNorm, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, TwoNorm, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, TwoNorm, MVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, clone, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv = Teuchos::rcp(new FType(grid, 1));

  auto mv2 = mv->clone(10);

  int n1(mv->getNumberVecs());
  int n2(mv2->getNumberVecs());

  TEST_EQUALITY(1, n1);
  TEST_EQUALITY(10, n2);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, clone, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, clone, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, clone, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, clone, MVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, CloneCopy, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv = Teuchos::rcp(new FType(grid, 10));

  mv->random();
  auto mv2 = mv->CloneCopy();

  int n1(mv->getNumberVecs());
  int n2(mv2->getNumberVecs());

  TEST_EQUALITY(10, n1);
  TEST_EQUALITY(n1, n2);

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv->norm(norm1);
  mv2->norm(norm2);
  for(int i=0; i<n1; ++i)
    TEST_FLOATING_EQUALITY(norm1[i], norm2[i], eps);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneCopy, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneCopy, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneCopy, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneCopy, MVF)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, CloneCopy2, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv = Teuchos::rcp(new FType(grid, 10));

  mv->random();

  std::vector<int> index(5);
  for(int i=0; i<5; ++i)
    index[i] = 2*i;

  auto mv2 = mv->CloneCopy(index);

  unsigned int n1 = (mv->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());

  TEST_EQUALITY(10, n1);
  TEST_EQUALITY(5, n2);
  TEST_EQUALITY(index.size(), n2);

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv->norm(norm1);
  mv2->norm(norm2);

  for(unsigned int i=0; i<index.size(); ++i)
    TEST_FLOATING_EQUALITY(norm1[index[i]], norm2[i], eps);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneCopy2, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneCopy2, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneCopy2, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneCopy2, MVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, CloneCopy3, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv1 = Teuchos::rcp(new FType(grid, 10));

  mv1->random();

  Teuchos::Range1D range(2, 7);

  auto mv2 = mv1->CloneCopy(range);

  unsigned int n1(mv1->getNumberVecs());
  unsigned int n2(mv2->getNumberVecs());

  TEST_EQUALITY(10, n1);
  TEST_EQUALITY(range.size(), n2);

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv1->norm(norm1);
  mv2->norm(norm2);
  for(int i=0; i<range.size(); ++i)
    TEST_FLOATING_EQUALITY(norm1[i+range.lbound()], norm2[i], eps);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneCopy3, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneCopy3, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneCopy3, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneCopy3, MVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, CloneViewNonConst1, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv1 = Teuchos::rcp(new FType(grid, 10));

  mv1->init(0.);

  std::vector<int> index(5);
  for(int i=0; i<5; ++i)
    index[i] = 2*i;

  auto mv2 = mv1->CloneViewNonConst(index);

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());

  TEST_EQUALITY(10, n1);
  TEST_EQUALITY(5, n2);
  TEST_EQUALITY(index.size(), n2);

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv2->random();

  mv1->norm(norm1);
  mv2->norm(norm2);

  for(unsigned int i=0; i<index.size(); ++i)
    TEST_FLOATING_EQUALITY(norm1[index[i]], norm2[i], eps);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneViewNonConst1, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneViewNonConst1, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneViewNonConst1, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneViewNonConst1, MVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, CloneViewNonConst2, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv1 = Teuchos::rcp(new FType(grid, 10));

  mv1->init(0.);

  Teuchos::Range1D range(2, 7);

  auto mv2 = mv1->CloneViewNonConst(range);

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());

  TEST_EQUALITY(10, n1);
  TEST_EQUALITY(range.size(), n2);

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv2->random();

  mv1->norm(norm1);
  mv2->norm(norm2);

  for(unsigned int i=0; i<range.size(); ++i)
    TEST_EQUALITY(norm1[i+range.lbound()], norm2[i]);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneViewNonConst2, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneViewNonConst2, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneViewNonConst2, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneViewNonConst2, MVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, CloneView1, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv1 = Teuchos::rcp(new FType(grid, 10));

  mv1->init(0.);

  std::vector<int> index(5);
  for(int i=0; i<5; ++i)
    index[i] = 2*i;

  auto mv2 = mv1->CloneView(index);

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());

  TEST_EQUALITY(10, n1);
  TEST_EQUALITY(5, n2);
  TEST_EQUALITY(index.size(), n2);

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  //			mv2->Random(); //<this should give compile error
  mv1->random();

  mv1->norm(norm1);
  mv2->norm(norm2);

  for(unsigned int i=0; i<index.size(); ++i)
    TEST_EQUALITY(norm1[index[i]], norm2[i]);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneView1, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneView1, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneView1, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneView1, MVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, CloneView2, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv1 = Teuchos::rcp(new FType(grid, 10));

  mv1->init(0.);

  Teuchos::Range1D range(2, 7);

  auto mv2 = mv1->CloneView(range);

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());

  TEST_EQUALITY(10, n1);
  TEST_EQUALITY(range.size(), n2);

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  //mv2->random(); // has to give compilation error
  mv1->random();

  mv1->norm(norm1);
  mv2->norm(norm2);

  for(unsigned int i=0; i<range.size(); ++i)
    TEST_EQUALITY(norm1[i+range.lbound()], norm2[i]);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneView2, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneView2, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneView2, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, CloneView2, MVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, TimesMatAdd, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv1 = Teuchos::rcp(new FType(grid, 10));

  mv1->init(0.);

  Teuchos::Range1D index1(0, 9);
  std::vector<int> index2(10);
  for(int i=0; i<10; ++i)
    index2[i] = i;

  auto mv2 = mv1->CloneView(index1);
  auto mv3 = mv1->clone(mv1->getNumberVecs());

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());
  unsigned int n3 = (mv3->getNumberVecs());

  TEST_EQUALITY(n1, n2);
  TEST_EQUALITY(n3, n2);
  TEST_EQUALITY(n3, n1);

  //  mv2->Random(); // has to give compilation error
  mv1->init(1.);
  //  mv2->Random();
  //	mv3->Assign(2.);

  Teuchos::SerialDenseMatrix<int, double> B(n1, n2);

  for(unsigned int j=0; j<n1; ++j)
    for(unsigned int i=0; i<n1; ++i)
      B(j, i) = 1./n1;

  mv1->TimesMatAdd(0.5, *mv2, B, 0.5);

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv1->norm(norm1);
  mv2->norm(norm2);

  for(unsigned int i=0; i<n1; ++i) {
    TEST_FLOATING_EQUALITY(norm1[i], norm2[i], eps);
    TEST_FLOATING_EQUALITY(std::sqrt((double)mv1->getLength()), norm2[i], 3.e-1);
  }

  std::vector<double> scales(n1);
  for(unsigned int j=0; j<n1; ++j){
    scales[j] = (j+1);
    for(unsigned int i=0; i<n1; ++i)
      B(j, i) = 1./n1/(j+1);
  }
  mv1->init(1.);
  mv3->init(1.);

  mv3->scale(scales);

  mv1->TimesMatAdd(1., *mv3, B, 0.);

  mv1->norm(norm1);
  mv2->norm(norm2);

  for(unsigned int i=0; i<n1; ++i) {
    TEST_FLOATING_EQUALITY(std::sqrt((double)mv1->getLength()), norm2[i], eps);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, TimesMatAdd, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, TimesMatAdd, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, TimesMatAdd, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, TimesMatAdd, MVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, add, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv1 = Teuchos::rcp(new FType(grid, 10));

  mv1->init(0.);

  Teuchos::Range1D range(0, 9);
  std::vector<int> index(10);
  for(int i=0; i<10; ++i)
    index[i] = i;

  auto mv2 = mv1->CloneViewNonConst(range);
  auto mv3 = mv1->CloneView(range);

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());
  unsigned int n3 = (mv3->getNumberVecs());

  TEST_EQUALITY(n1, n2);
  TEST_EQUALITY(n3, n2);
  TEST_EQUALITY(n3, n1);

  mv1->init(1.);
  mv2->init(1.);

  mv1->add(0.5, *mv2, 0.5, *mv3);

  std::vector<double> norm1(n1);
  std::vector<double> norm2(n2);

  mv1->norm(norm1);
  mv2->norm(norm2);

  for(unsigned int i=0; i<n1; ++i)
    TEST_FLOATING_EQUALITY(std::sqrt((double)mv1->getLength()), norm2[i], eps);

  mv1->init(1.);
  mv2->init(1.);

  mv2->scale(0.5);

  mv1->add(1., *mv2, 1., *mv3);

  mv1->norm(norm1);
  *mv2 = *mv1;
  mv2->norm(norm2);

  for(unsigned int i=0; i<n1; ++i) {
    TEST_FLOATING_EQUALITY(norm1[i], norm2[i], eps);
    TEST_FLOATING_EQUALITY(std::sqrt((double)mv1->getLength()), norm1[i], eps);
  }

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, add, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, add, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, add, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, add, MVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, dot, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv1 = Teuchos::rcp(new FType(grid, 10));

  mv1->init(0.);

  Teuchos::Range1D index1(0, 9);
  std::vector<int> index2(10);
  for(int i=0; i<10; ++i)
    index2[i] = i;

  auto mv2 = mv1->CloneViewNonConst(index1);
  auto mv3 = mv1->CloneView(index1);

  unsigned int n1 = (mv1->getNumberVecs());
  unsigned int n2 = (mv2->getNumberVecs());
  unsigned int n3 = (mv3->getNumberVecs());

  TEST_EQUALITY(n1, n2);
  TEST_EQUALITY(n3, n2);
  TEST_EQUALITY(n3, n1);

  mv1->init(1.);
  mv2->init(1.);

  std::vector<double> dots(n1);

  mv1->dot(*mv2, dots);


  for(unsigned int i=0; i<n1; ++i) {
    TEST_EQUALITY(mv1->getLength(), dots[i]);
  }

  mv2->init(2.);

  mv3->dot(*mv2, dots);

  for(unsigned int i=0; i<n1; ++i) {
    TEST_EQUALITY(4*mv1->getLength(), dots[i]);
  }

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, dot, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, dot, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, dot, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, dot, MVF)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MultiField, Trans, FType) {

  Pimpact::setBoundaryConditions(pl, domain);

  // processor grid size
  pl->set("npx", (2==FType::GridT::sdim)?4:2);
  pl->set("npy",            2);
  pl->set("npz", (2==FType::GridT::sdim)?1:2);

  auto grid = Pimpact::create<GridT>(pl);

  auto mv1 = Teuchos::rcp(new FType(grid, 10));

  mv1->init(0.);

  Teuchos::Range1D range(0, 9);
  std::vector<int> index(10);
  for(int i=0; i<10; ++i)
    index[i] = i;

  auto mv2 = mv1->CloneView(range);
  auto mv3 = mv1->CloneView(range);

  unsigned int n1 = mv1->getNumberVecs();
  unsigned int n2 = mv2->getNumberVecs();
  unsigned int n3 = mv3->getNumberVecs();

  TEST_EQUALITY(n1, n2);
  TEST_EQUALITY(n3, n2);
  TEST_EQUALITY(n3, n1);

  mv1->init(1.);

  Teuchos::SerialDenseMatrix<int, double> B(n1, n2);

  mv1->Trans(1., *mv2, B);

  for(unsigned int j=0; j<n1; ++j) {
    for(unsigned int i=0; i<n1; ++i)
      TEST_EQUALITY(mv1->getLength(), B(j, i));
  }

  std::vector<double> scales(n1);

  for(unsigned int i=0; i<scales.size(); ++i)
    scales[i] = i*2;

  mv1->scale(scales);

  mv2->Trans(1., *mv3, B);

  for(unsigned int j=0; j<n1; ++j) {
    for(unsigned int i=0; i<n1; ++i)
      TEST_EQUALITY(scales[i]*scales[j]*mv1->getLength(), B(j, i));
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, Trans, SF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, Trans, VF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, Trans, MSF)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MultiField, Trans, MVF)


} // end of namespace
