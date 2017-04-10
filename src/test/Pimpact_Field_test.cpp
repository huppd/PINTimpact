#include <cmath>
#include <iostream>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Pimpact_Fields.hpp"

#include "Pimpact_Test.hpp"



namespace {



using SF2D  = typename Pimpact::ScalarField<D2>;
using VF2D  = typename Pimpact::VectorField<D2>;
using MSF2D = typename Pimpact::ModeField<SF2D>;
using MVF2D = typename Pimpact::ModeField<VF2D>;
using CF2D  = typename Pimpact::CompoundField<VF2D,SF2D>;
using CMF2D = typename Pimpact::CompoundField<MVF2D,MSF2D>;

using SF3D  = typename Pimpact::ScalarField<D3>;
using VF3D  = typename Pimpact::VectorField<D3>;
using MSF3D = typename Pimpact::ModeField<SF3D>;
using MVF3D = typename Pimpact::ModeField<VF3D>;
using CF3D  = typename Pimpact::CompoundField<VF3D,SF3D>;
using CMF3D = typename Pimpact::CompoundField<MVF3D,MSF3D>;




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, print, FType ) {

	setParameter( FType::SpaceT::sdim );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

	p->init( space->rankST() );
	p->print();
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, print, CMF3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, InfNormAndInit, FType ) {

	setParameter( FType::SpaceT::sdim );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  ST norm;

  // test different float values, assures that initial and norm work smoothly
  for( ST i=0.; i< 10.1; ++i ) {
    p->init(i/2.);
    norm = p->norm(Belos::InfNorm);
    TEST_FLOATING_EQUALITY( i/2., norm, eps );

  }

  // one test with infty-norm
  int rank;
  int size;
  ST init;
  MPI_Comm_rank(space->commST(),&rank);
  MPI_Comm_size(space->commST(),&size);
  for( ST i = 0.; i<10.1; ++i) {
    init = ( size-1 )*i-1.;
    init = std::abs( init );
    p->init( rank*i-1. );
    norm = p->norm( Belos::InfNorm );
    TEST_FLOATING_EQUALITY( std::max(init,1.), norm, eps );

  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, SF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, VF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, MSF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, MVF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, CF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, CMF2D ) 

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, SF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, VF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, MSF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, MVF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, CF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, InfNormAndInit, CMF3D ) 



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, OneNormAndInit, FType ) {

	setParameter( FType::SpaceT::sdim );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  // test different float values, assures that initial and norm work smoothly
  for( ST i=0.; i< 10.1; ++i ) {
    p->init(i/2.);
//    TEST_EQUALITY( (i/2.)*p->getLength(), p->norm(Belos::OneNorm) );
    TEST_FLOATING_EQUALITY( (i/2.)*p->getLength(), p->norm(Belos::OneNorm), eps );

  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, SF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, MSF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, CF2D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, SF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, MSF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, CF3D ) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, OneNormAndInit, CMF3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, TwoNormAndInit, FType ) {

	setParameter( FType::SpaceT::sdim );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  // test different float values, assures that initial and norm work smoothly
  for( ST i=0.; i< 10.1; ++i ) {
    p->init(i/2.);
    TEST_FLOATING_EQUALITY( std::sqrt( std::pow(i/2.,2)*p->getLength() ), p->norm(Belos::TwoNorm), eps );
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, TwoNormAndInit, CMF3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, dot, FType ) {

	setParameter( FType::SpaceT::sdim );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto vel1 = Pimpact::create<FType>(space);
  auto vel2 = Pimpact::create<FType>(space);

  int Np = vel1->getLength();
  int Nq = vel2->getLength();
  ST dot;

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

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, dot, CMF3D )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, scale, FType ) {

	setParameter( FType::SpaceT::sdim );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  ST norm;
  int N = p->getLength();

  p->init(1.);
  p->scale(2.);
  norm = p->norm(Belos::TwoNorm);
  TEST_EQUALITY( std::sqrt(4*N), norm)

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, scale, CMF3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, random, FType ) {

	setParameter( FType::SpaceT::sdim );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  ST norm;
  int N = p->getLength();

  p->init(1.);
  p->random();
  norm = p->norm(Belos::TwoNorm);
  TEST_INEQUALITY( N, norm)

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, random, CMF3D )
	


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TemplateField, add, FType ) {

	setParameter( FType::SpaceT::sdim );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto vel1 = Pimpact::create<FType>(space);
  auto vel2 = Pimpact::create<FType>(space);
  auto vel3 = Pimpact::create<FType>(space);

  TEST_EQUALITY( vel1->getLength(), vel2->getLength() )
  TEST_EQUALITY( vel2->getLength(), vel3->getLength() )
  TEST_EQUALITY( vel1->getLength(), vel3->getLength() )

  ST norm;
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

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TemplateField, add, CMF3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, write, FType ) {

	setParameter( FType::SpaceT::sdim );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>( space );

  p->init(1.);
  p->write();

  p->random();
  p->write(1);

  TEST_EQUALITY( 0, 0 )

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, VF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, MSF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, MVF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, CF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, CMF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, VF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, MSF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, MVF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, CF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, write, CMF3D )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TempField, writeRestart, FType ) {

	setParameter( FType::SpaceT::sdim );

  Teuchos::RCP<const typename FType::SpaceT> space = Pimpact::create<typename FType::SpaceT>( pl );

  auto p = Pimpact::create<FType>(space);

  p->init(1.);
  p->write();

  p->random();
  p->write( 99, true );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, writeRestart, SF2D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, writeRestart, VF2D )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, writeRestart, SF3D )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TempField, writeRestart, VF3D )
	

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ScalarField, ReadWrite, SpaceT ) {

	setParameter( SpaceT::sdim );

  Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl );

	Pimpact::ScalarField<SpaceT> write( space );
	Pimpact::ScalarField<SpaceT> read( space );
	Pimpact::ScalarField<SpaceT> err( space );

  for( int i=1; i<=6; ++i ) {
    write.initField( static_cast<Pimpact::EScalarField>(i) );
    write.write( i, true );
		read.read(i);

		err.add( 1., read, -1., write );
		ST error = err.norm();
		if( 0==space->rankST() ) std::cout << "\nerror: " << error << "\n";
		TEST_EQUALITY( std::abs(error)<eps, true );
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ScalarField, ReadWrite, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ScalarField, ReadWrite, D3 )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ScalarField, level, SpaceT ) {

	setParameter( SpaceT::sdim );

  Teuchos::RCP<const SpaceT> space = Pimpact::create<SpaceT>( pl );

	Pimpact::ScalarField<SpaceT> x( space );

	x.init( 1. );
	x.level();

	ST level = x.norm();
	if( 0==space()->rankST() )
		std::cout << "\nlevel: " << level << "\n";
  TEST_EQUALITY( level<eps , true );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ScalarField, level, D2 )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ScalarField, level, D3 )

} // namespace

