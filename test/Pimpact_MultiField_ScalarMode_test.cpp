// Pimpact_SalarVectorSpace_test.cpp

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include "Teuchos_Range1D.hpp"
#include <Teuchos_CommHelpers.hpp>

#include "pimpact.hpp"
#include "Pimpact_FieldSpace.hpp"
#include "Pimpact_IndexSpace.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiField.hpp"

#include <iostream>

namespace {

bool testMpi = true;
double errorTolSlack = 1e+1;

TEUCHOS_STATIC_SETUP() {
	Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
	clp.addOutputSetupOptions(true);
	clp.setOption(
			"test-mpi", "test-serial", &testMpi,
			"Test MPI (if available) or force test of serial.  In a serial build,"
			" this option is ignored and a serial comm is always used." );
	clp.setOption(
			"error-tol-slack", &errorTolSlack,
			"Slack off of machine epsilon used to check test results" );
}


// test shows that nLoc is not consistent with start and end indexes
TEUCHOS_UNIT_TEST( MultiField, constructor ) {
	// init impact
	init_impact(0,0);

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>(*vel,10);

	const int m = mv->GetNumberVecs();

	TEST_EQUALITY( 10, m );
}


TEUCHOS_UNIT_TEST( MultiField, TwoNorm_and_init ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>(*vel,10);

	const int m = mv->GetNumberVecs();
	const int n = mv->GetVecLength();
	std::vector<double> normval(m);

	// test different float values, assures that initial and norm work smoothly
	for( double i=0.; i< 200.1; ++i ) {
		mv->Init(i/2.);
		mv->Norm(normval,Belos::TwoNorm);
//			mv->Assign(*mv);
//			mv->Dot(*mv,normval);
		for( int j=0; j<m; ++j )
			TEST_EQUALITY( std::pow(i/2.,2)*n, normval[j] );
	}
}


TEUCHOS_UNIT_TEST( MultiField, clone ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>( *vel, 1 );

	auto mv2 = mv->Clone(10);

	int m1(mv->GetNumberVecs());
	int m2(mv2->GetNumberVecs());

	int n1(mv->GetVecLength());
	int n2(mv2->GetVecLength());

	TEST_EQUALITY( 1, m1 );
	TEST_EQUALITY( 10, m2 );

	TEST_EQUALITY( n1, n2);

}


TEUCHOS_UNIT_TEST( MultiField, CloneCopy ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>(*vel,10);

	mv->Random();
	auto mv2 = mv->CloneCopy();

	int n1(mv->GetNumberVecs());
	int n2(mv2->GetNumberVecs());

	TEST_EQUALITY( 10, n1 );
	TEST_EQUALITY( n1, n2 );

	int m1(mv->GetVecLength());
	int m2(mv2->GetVecLength());
	TEST_EQUALITY( m1, m2);


	std::vector<double> norm1(n1);
	std::vector<double> norm2(n2);

	mv->Norm(norm1);
	mv2->Norm(norm2);
	for( int i=0; i<n1; ++i)
		TEST_EQUALITY( norm1[i], norm2[i] );
}


TEUCHOS_UNIT_TEST( MultiField, CloneCopy2 ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>(*vel,10);

	mv->Random();

	std::vector<int> index(5);
	for(int i=0; i<5; ++i)
		index[i] = 2*i;

	auto mv2 = mv->CloneCopy(index);

	unsigned int n1 = (mv->GetNumberVecs());
	unsigned int n2 = (mv2->GetNumberVecs());

	TEST_EQUALITY( 10, n1 );
	TEST_EQUALITY( 5, n2 );
	TEST_EQUALITY( index.size(), n2 );

	int m1(mv->GetVecLength());
	int m2(mv2->GetVecLength());
	TEST_EQUALITY( m1, m2);

	std::vector<double> norm1(n1);
	std::vector<double> norm2(n2);

	mv->Norm(norm1);
	mv2->Norm(norm2);

	for( unsigned int i=0; i<index.size(); ++i)
		TEST_EQUALITY( norm1[index[i]], norm2[i] );
}


TEUCHOS_UNIT_TEST( MultiField, CloneCopy3 ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>(*vel,10);

	mv1->Random();

	Teuchos::Range1D index(2,7);

	auto mv2 = mv1->CloneCopy(index);

	unsigned int n1(mv1->GetNumberVecs());
	unsigned int n2(mv2->GetNumberVecs());

	TEST_EQUALITY( 10, n1 );
	TEST_EQUALITY( index.size(), n2 );

	int m1(mv1->GetVecLength());
	int m2(mv2->GetVecLength());
	TEST_EQUALITY( m1, m2);


	std::vector<double> norm1(n1);
	std::vector<double> norm2(n2);

	mv1->Norm(norm1);
	mv2->Norm(norm2);
	for( int i=0; i<index.size(); ++i)
		TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );

}


TEUCHOS_UNIT_TEST( MultiField, CloneViewNonConst1 ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>(*vel,10);

	mv1->Init(0.);

	std::vector<int> index(5);
	for(int i=0; i<5; ++i)
		index[i] = 2*i;

	auto mv2 = mv1->CloneViewNonConst(index);

	unsigned int n1 = (mv1->GetNumberVecs());
	unsigned int n2 = (mv2->GetNumberVecs());
//
	TEST_EQUALITY( 10, n1 );
	TEST_EQUALITY( 5, n2 );
	TEST_EQUALITY( index.size(), n2 );

	int m1(mv1->GetVecLength());
	int m2(mv2->GetVecLength());
	TEST_EQUALITY( m1, m2);


	std::vector<double> norm1(n1);
	std::vector<double> norm2(n2);

	mv2->Random();

	mv1->Norm(norm1);
	mv2->Norm(norm2);

	for( unsigned int i=0; i<index.size(); ++i)
		TEST_EQUALITY( norm1[index[i]], norm2[i] );
}


TEUCHOS_UNIT_TEST( MultiField, CloneViewNonConst2 ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>(*vel,10);

	mv1->Init(0.);

	Teuchos::Range1D index(2,7);

	auto mv2 = mv1->CloneViewNonConst(index);

	unsigned int n1 = (mv1->GetNumberVecs());
	unsigned int n2 = (mv2->GetNumberVecs());

	TEST_EQUALITY( 10, n1 );
	TEST_EQUALITY( index.size(), n2 );

	int m1(mv1->GetVecLength());
	int m2(mv2->GetVecLength());
	TEST_EQUALITY( m1, m2);


	std::vector<double> norm1(n1);
	std::vector<double> norm2(n2);

	mv2->Random();

	mv1->Norm(norm1);
	mv2->Norm(norm2);

	for( unsigned int i=0; i<index.size(); ++i)
		TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
}


TEUCHOS_UNIT_TEST( MultiField, CloneView1 ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>(*vel,10);

	mv1->Init(0.);

	std::vector<int> index(5);
	for(int i=0; i<5; ++i)
		index[i] = 2*i;

	auto mv2 = mv1->CloneView(index);

	unsigned int n1 = (mv1->GetNumberVecs());
	unsigned int n2 = (mv2->GetNumberVecs());

	TEST_EQUALITY( 10, n1 );
	TEST_EQUALITY( 5, n2 );
	TEST_EQUALITY( index.size(), n2 );


	int m1(mv1->GetVecLength());
	int m2(mv2->GetVecLength());
	TEST_EQUALITY( m1, m2);

	std::vector<double> norm1(n1);
	std::vector<double> norm2(n2);

	mv1->Random();

	mv1->Norm(norm1);
	mv2->Norm(norm2);

	for( unsigned int i=0; i<index.size(); ++i)
		TEST_EQUALITY( norm1[index[i]], norm2[i] );
}


TEUCHOS_UNIT_TEST( MultiField, CloneViewt2 ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>(*vel,10);

	mv1->Init(0.);

	Teuchos::Range1D index(2,7);

	auto mv2 = mv1->CloneView(index);

	unsigned int n1 = (mv1->GetNumberVecs());
	unsigned int n2 = (mv2->GetNumberVecs());

	TEST_EQUALITY( 10, n1 );
	TEST_EQUALITY( index.size(), n2 );


	int m1(mv1->GetVecLength());
	int m2(mv2->GetVecLength());
	TEST_EQUALITY( m1, m2);

	std::vector<double> norm1(n1);
	std::vector<double> norm2(n2);

//			mv2->Random(); // has to give compilation error
	mv1->Random();

	mv1->Norm(norm1);
	mv2->Norm(norm2);

	for( unsigned int i=0; i<index.size(); ++i)
		TEST_EQUALITY( norm1[i+index.lbound()], norm2[i] );
}


TEUCHOS_UNIT_TEST( MultiField, TimesMatAdd ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>(*vel,10);

	mv1->Init(0.);

	Teuchos::Range1D index1(0,9);
	std::vector<int> index2(10);
	for(int i=0; i<10; ++i)
		index2[i] = i;

	auto mv2 = mv1->CloneView(index1);
	auto mv3 = mv1->Clone(mv1->GetNumberVecs() );


	unsigned int n1 = (mv1->GetNumberVecs());
	unsigned int n2 = (mv2->GetNumberVecs());
	unsigned int n3 = (mv3->GetNumberVecs());

	TEST_EQUALITY( n1, n2 );
	TEST_EQUALITY( n3, n2 );
	TEST_EQUALITY( n3, n1 );


	int m1(mv1->GetVecLength());
	int m2(mv2->GetVecLength());
	int m3(mv3->GetVecLength());
	TEST_EQUALITY( m1, m2);
	TEST_EQUALITY( m2, m3);


	//			mv2->Random(); // has to give compilation error
	mv1->Init(1.);
	//				mv2->Random();
	//				mv3->Assign(2.);


	Teuchos::SerialDenseMatrix<int,double> B(n1,n2);
	for( unsigned int j=0; j<n1; ++j)
		for( unsigned int i=0; i<n1; ++i)
			B(j,i) = 1./n1;

	mv1->TimesMatAdd( 0.5, *mv2, B, 0.5 );


	std::vector<double> norm1(n1);
	std::vector<double> norm2(n2);

	mv1->Norm(norm1);
	mv2->Norm(norm2);

	for( unsigned int i=0; i<n1; ++i) {
		TEST_EQUALITY( norm1[i], norm2[i] );
		TEST_EQUALITY( mv1->GetVecLength(), norm2[i] );
	}

	std::vector<double> scales(n1);
	for( unsigned int j=0; j<n1; ++j){
		scales[j] = (j+1);
	//					scales[j] = 0;
		for( unsigned int i=0; i<n1; ++i)
			B(j,i) = 1./n1/(j+1);
	//						B(j,i) = 1./5.;
	}
	//				std::cout << B;
	//				scales[5] = 5.;
	mv1->Init(1.);
	mv3->Init(1.);

	mv3->Scale(scales);

	mv1->TimesMatAdd( 1., *mv3, B, 0. );


	mv1->Norm(norm1);
	mv2->Norm(norm2);

	for( unsigned int i=0; i<n1; ++i) {
	//					TEST_EQUALITY( norm1[i], norm2[i] );
		TEST_FLOATING_EQUALITY( (double)mv1->GetVecLength(), norm2[i], 0.01 );
	}
}


TEUCHOS_UNIT_TEST( MultiField, Add ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>(*vel,10);

	mv1->Init(0.);

	Teuchos::Range1D index1(0,9);
	std::vector<int> index2(10);
	for(int i=0; i<10; ++i)
		index2[i] = i;

//					auto mv2 = mv1->CloneCopy(index1);
	auto mv2 = mv1->CloneViewNonConst(index1);
	auto mv3 = mv1->CloneView(index1);

	unsigned int n1 = (mv1->GetNumberVecs());
	unsigned int n2 = (mv2->GetNumberVecs());
	unsigned int n3 = (mv3->GetNumberVecs());

	TEST_EQUALITY( n1, n2 );
	TEST_EQUALITY( n3, n2 );
	TEST_EQUALITY( n3, n1 );

	int m1(mv1->GetVecLength());
	int m2(mv2->GetVecLength());
	int m3(mv3->GetVecLength());
	TEST_EQUALITY( m1, m2);
	TEST_EQUALITY( m2, m3);

	mv1->Init(1.);
	mv2->Init(1.);

	mv1->Add( 0.5, *mv2, 0.5, *mv3);

	std::vector<double> norm1(n1);
	std::vector<double> norm2(n2);

	mv1->Norm(norm1);
	mv2->Norm(norm2);

	for( unsigned int i=0; i<n1; ++i) {
		TEST_EQUALITY( mv1->GetVecLength(), norm2[i] );
	}

	mv1->Init(1.);
	mv2->Init(1.);

	mv2->Scale(0.5);

	mv1->Add( 1., *mv2, 1., *mv3 );


	mv1->Norm(norm1);
	mv2->Assign(*mv1);
	mv2->Norm(norm2);

	for( unsigned int i=0; i<n1; ++i) {
		TEST_EQUALITY( norm1[i], norm2[i] );
		TEST_FLOATING_EQUALITY( (double)mv1->GetVecLength(), norm1[i], 0.01 );
	}
}


TEUCHOS_UNIT_TEST( MultiField, Dot ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>(*vel,10);

	mv1->Init(0.);

	Teuchos::Range1D index1(0,9);
	std::vector<int> index2(10);
	for(int i=0; i<10; ++i)
	index2[i] = i;

	//					auto mv2 = mv1->CloneCopy(index1);
	auto mv2 = mv1->CloneViewNonConst(index1);
	auto mv3 = mv1->CloneView(index1);

	unsigned int n1 = (mv1->GetNumberVecs());
	unsigned int n2 = (mv2->GetNumberVecs());
	unsigned int n3 = (mv3->GetNumberVecs());

	TEST_EQUALITY( n1, n2 );
	TEST_EQUALITY( n3, n2 );
	TEST_EQUALITY( n3, n1 );

	int m1(mv1->GetVecLength());
	int m2(mv2->GetVecLength());
	int m3(mv3->GetVecLength());

	TEST_EQUALITY( m1, m2);
	TEST_EQUALITY( m2, m3);

	mv1->Init(1.);
	mv2->Init(1.);

	std::vector<double> dots(n1);

	mv1->Dot( *mv2, dots );


	for( unsigned int i=0; i<n1; ++i) {
		TEST_EQUALITY( mv1->GetVecLength(), dots[i] );
	}

	mv2->Init(2.);

	mv3->Dot( *mv2, dots  );
	for( unsigned int i=0; i<n1; ++i) {
		TEST_EQUALITY( 4*mv1->GetVecLength(), dots[i] );
	}
}


TEUCHOS_UNIT_TEST( MultiField, Trans ) {

	auto sVS = Pimpact::createFieldSpace<int>();

	auto pc = Pimpact::createScalarField<double,int>(sVS);
	auto ps = Pimpact::createScalarField<double,int>(sVS);

	auto vel = Pimpact::createModeField( pc, ps );

	auto mv1 = Pimpact::createMultiField<Pimpact::ModeField<Pimpact::ScalarField<double,int> >,double,int>(*vel,10);

	mv1->Init(0.);

	Teuchos::Range1D index1(0,9);
	std::vector<int> index2(10);
	for(int i=0; i<10; ++i)
		index2[i] = i;

	auto mv2 = mv1->CloneView(index1);
	auto mv3 = mv1->CloneView(index1);

	unsigned int n1 = (mv1->GetNumberVecs());
	unsigned int n2 = (mv2->GetNumberVecs());
	unsigned int n3 = (mv3->GetNumberVecs());

	TEST_EQUALITY( n1, n2 );
	TEST_EQUALITY( n3, n2 );
	TEST_EQUALITY( n3, n1 );


	mv1->Init(1.);

	Teuchos::SerialDenseMatrix<int,double> B(n1,n2);

	mv1->Trans( 1., *mv2, B );

	for( unsigned int j=0; j<n1; ++j){
		for( unsigned int i=0; i<n1; ++i)
			TEST_EQUALITY( mv1->GetVecLength(), B(j,i) );
	}

	std::vector<double> scales(n1);
	for( unsigned int i=0; i<scales.size(); ++i)
		scales[i] = i*2;
	mv1->Scale(scales);

	mv2->Trans( 1., *mv3, B );

	for( unsigned int j=0; j<n1; ++j){
		for( unsigned int i=0; i<n1; ++i)
			TEST_EQUALITY( scales[i]*scales[j]*mv1->GetVecLength(), B(j,i) );
	}
}


} // namespace

