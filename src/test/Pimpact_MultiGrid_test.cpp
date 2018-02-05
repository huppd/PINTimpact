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


const ST pi2 = 2.*std::acos(-1.);

using CGrid2DT = Pimpact::Grid<ST, OT, 2, d, 2>;
using CGrid3DT = Pimpact::Grid<ST, OT, 3, d, 2>;

using CSL2D = Pimpact::CoarsenStrategy<D2, CGrid2DT>;
using CSG2D = Pimpact::CoarsenStrategyGlobal<D2, CGrid2DT>;

using CSL3D = Pimpact::CoarsenStrategy<D3, CGrid3DT>;
using CSG3D = Pimpact::CoarsenStrategyGlobal<D3, CGrid3DT>;

using CCSL3D = Pimpact::CoarsenStrategy<D3, D3>;
using CCSG3D = Pimpact::CoarsenStrategyGlobal<D3, D3>;

using Grids2D = Pimpact::MGGrids<D2, CGrid2DT>;
using Grids3D = Pimpact::MGGrids<D3, CGrid3DT>;

using CGrids2D = Pimpact::MGGrids<D2, D2>;
using CGrids3D = Pimpact::MGGrids<D3, D3>;


template<class T> using ConvDiffOpT =
	Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<T> >;

template<class T> using ConvDiffSORT =
	Pimpact::NonlinearSmoother<T, Pimpact::ConvectionDiffusionSORSmoother >;

template<class T> using ConvDiffJT =
	Pimpact::NonlinearSmoother<T, Pimpact::ConvectionDiffusionJSmoother >;

template<class T1, class T2> using TransVF = Pimpact::VectorFieldOpWrap<Pimpact::TransferOp<T1, T2> >;
template<class T> using RestrVF = Pimpact::VectorFieldOpWrap<Pimpact::RestrictionVFOp<T> >;
template<class T> using InterVF = Pimpact::VectorFieldOpWrap<Pimpact::InterpolationOp<T> >;

template<class T> using MOP = Pimpact::InverseOp<T >;
template<class T> using POP = Pimpact::PrecInverseOp<T, Pimpact::DivGradO2JSmoother >;


using DGJMGT2D = Pimpact::MultiGrid<
	Grids2D,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::DivGradO2JSmoother,
	Pimpact::DivGradO2Inv >;

using DGJMGT3D = Pimpact::MultiGrid<
	Grids3D,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::DivGradO2JSmoother,
	Pimpact::DivGradO2Inv >;



using DGLMGT2D = Pimpact::MultiGrid<
	Grids2D,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::DivGradO2LSmoother,
	Pimpact::DivGradO2Inv >;

using DGLMGT3D = Pimpact::MultiGrid<
	Grids3D,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::DivGradO2LSmoother,
	Pimpact::DivGradO2Inv >;

using DGCMGT2D = Pimpact::MultiGrid<
	Grids2D,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::Chebyshev,
	Pimpact::DivGradO2Inv >;

using DGCMGT3D = Pimpact::MultiGrid<
	Grids3D,
	Pimpact::ScalarField,
	Pimpact::TransferOp,
	Pimpact::RestrictionSFOp,
	Pimpact::InterpolationOp,
	Pimpact::DivGradOp,
	Pimpact::DivGradO2Op,
	Pimpact::Chebyshev,
	Pimpact::DivGradO2Inv >;


using CDJ2D = Pimpact::MultiGrid<
	Grids2D,
	Pimpact::VectorField,
	TransVF,
	RestrVF,
	InterVF,
	ConvDiffOpT,
	ConvDiffOpT,
	ConvDiffJT,
	MOP >;

using CDJ3D = Pimpact::MultiGrid<
	Grids3D,
	Pimpact::VectorField,
	TransVF,
	RestrVF,
	InterVF,
	ConvDiffOpT,
	ConvDiffOpT,
	ConvDiffJT,
	MOP >;

using CDSOR2D = Pimpact::MultiGrid<
	Grids2D,
	Pimpact::VectorField,
	TransVF,
	RestrVF,
	InterVF,
	ConvDiffOpT,
	ConvDiffOpT,
	ConvDiffSORT,
	MOP >;

using CDSOR3D = Pimpact::MultiGrid<
	Grids3D,
	Pimpact::VectorField,
	TransVF,
	RestrVF,
	InterVF,
	ConvDiffOpT,
	ConvDiffOpT,
	ConvDiffSORT,
	MOP >;

using CCDSOR3D = Pimpact::MultiGrid<
	CGrids3D,
	Pimpact::VectorField,
	TransVF,
	RestrVF,
	InterVF,
	ConvDiffOpT,
	ConvDiffOpT,
	ConvDiffSORT,
	MOP >;


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MGGrids, constructor, CS) {

	setParameter(CS::GridT::sdim);

	auto grid = Pimpact::create<typename CS::GridT>(pl);

	int rank = grid->rankST();
	if(print) grid->print();

	{
		Teuchos::RCP<std::ostream> fstream = Pimpact::createOstream("coord_rank_"+std::to_string(rank)+".txt", 0);

		grid->getCoordinatesGlobal()->print(*fstream);
	}

	auto mgGrids = Pimpact::createMGGrids<CS>(grid, maxGrids);

	//std::cout <<"rank: " <<grid->rankST() <<"\tnGridLevels: " <<mgGrids->getNGrids() <<"\n";

	//if(grid->rankST()==0 && print)
		//mgGrids->print();

	//for(int level=0; level<mgGrids->getNGrids(); ++level) {
		//Teuchos::RCP<std::ostream> fstream = Pimpact::createOstream("coord_l"+std::to_string(level)+"rank_"+std::to_string(rank)+".txt", 0);
		//mgGrids->get(level)->getCoordinatesGlobal()->print(*fstream);
	//}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGGrids, constructor, CSG2D)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGGrids, constructor, CSL3D)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGGrids, constructor, CSG3D)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MGFields, SFconstructor, CS) {

	setParameter(CS::GridT::sdim);

	auto grid = Pimpact::create<typename CS::GridT>(pl);

	auto mgGrids = Pimpact::createMGGrids<CS>(grid, maxGrids);

	std::cout <<"nGridLevels: " <<mgGrids->getNGrids() <<"\n";
	if(grid->rankST()==0 && print)
		mgGrids->print();

	auto mgFields = Pimpact::createMGFields<Pimpact::ScalarField>(mgGrids);

	auto field = mgFields->get(-1);
	if(mgGrids->participating(-1)){

		field.init(5.);
		TEST_FLOATING_EQUALITY(std::sqrt(std::pow(5., 2)*field.getLength()), field.norm(Pimpact::ENorm::Two), eps);

		if(write) field.write(0);
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGFields, SFconstructor, CSL2D)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGFields, SFconstructor, CSL3D)

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGFields, SFconstructor, CSG2D)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGFields, SFconstructor, CSG3D)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MGFields, VFconstructor, CS) {

	setParameter(CS::GridT::sdim);

	auto grid = Pimpact::create<typename CS::GridT>(pl);

	auto mgGrids = Pimpact::createMGGrids<CS>(grid, maxGrids);

	std::cout <<"nGridLevels: " <<mgGrids->getNGrids() <<"\n";
	if(grid->rankST()==0 && print)
		mgGrids->print();

	auto mgFields = Pimpact::createMGFields<Pimpact::VectorField>(mgGrids);

	auto field = mgFields->get(-1);
	if(mgGrids->participating(-1)){

		field.init(5.);
		TEST_FLOATING_EQUALITY(std::sqrt(std::pow(5., 2)*field.getLength()), field.norm(Pimpact::ENorm::Two), eps);

		if(write) field.write(0);
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGFields, VFconstructor, CSL2D)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGFields, VFconstructor, CSL3D)

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGFields, VFconstructor, CSG2D)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGFields, VFconstructor, CSG3D)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MGOperators, SFconstructor, CS) {

	setParameter(CS::GridT::sdim);

	auto grid = Pimpact::create<typename CS::GridT>(pl);

	auto mgGrids = Pimpact::createMGGrids<CS>(grid, maxGrids);

  auto op = Pimpact::create<Pimpact::DivGradO2Op>(grid);

	auto mgOps = Pimpact::createMGOperators<Pimpact::DivGradO2Op>(mgGrids, op);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGOperators, SFconstructor, CSG2D)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGOperators, SFconstructor, CSG3D)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MGOperators, VFconstructor, CS) {

	setParameter(CS::GridT::sdim);

	auto grid = Pimpact::create<typename CS::GridT>(pl);

	auto mgGrids = Pimpact::createMGGrids<CS>(grid, maxGrids);

  auto op = Pimpact::create<ConvDiffOpT>(grid);

	auto mgOps = Pimpact::createMGOperators<ConvDiffOpT>(mgGrids, op);

	if(print && 0==grid->rankST())
		mgOps->print();
}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGOperators, VFconstructor, CSG2D)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGOperators, VFconstructor, CSG3D)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MGSmoothers, SFconstructor, CS) {

	setParameter(CS::GridT::sdim);

	auto grid = Pimpact::create<typename CS::GridT>(pl);

	auto mgGrids = Pimpact::createMGGrids<CS>(grid, maxGrids);

  auto op = Pimpact::create<Pimpact::DivGradOp>(grid);

	auto mgOps = Pimpact::createMGOperators<Pimpact::DivGradOp, Pimpact::DivGradO2Op>(mgGrids, op);

	auto mgSmoother = Pimpact::createMGSmoothers<Pimpact::DivGradO2JSmoother>(mgOps);

	auto smoother = mgSmoother->get(-2);

	if(mgGrids->participating(-2) && print)
		smoother->print();

}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGSmoothers, SFconstructor, CSG2D)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGSmoothers, SFconstructor, CSG3D)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MGSmoothers, VFconstructor, CS) {

	setParameter(CS::GridT::sdim);

	auto grid = Pimpact::create<typename CS::GridT>(pl);

	auto mgGrids = Pimpact::createMGGrids<CS>(grid, maxGrids);

  auto op = Pimpact::create<ConvDiffOpT>(grid);
	auto mgOps = Pimpact::createMGOperators<ConvDiffOpT>(mgGrids, op);

	auto mgSmoother = Pimpact::createMGSmoothers<ConvDiffSORT>(mgOps);

	auto smoother = mgSmoother->get(-2);

	if(mgGrids->participating(-2) && print)
		smoother->print();

}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGSmoothers, VFconstructor, CSG2D)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGSmoothers, VFconstructor, CSG3D)



TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(MGTransfers, Restrictor, CS, RestrictorType) {

	setParameter(CS::GridT::sdim);

	auto grid = Pimpact::create<typename CS::GridT>(pl);

	auto mgGrids = Pimpact::createMGGrids<CS>(grid, maxGrids);

	for(int level=1; level<mgGrids->getNGrids(); ++level) {
		if(0==grid->rankST()) {
			std::cout <<"\n\n\t--- level: " <<level-1 <<"---\n";
		}
		Teuchos::RCP<const RestrictorType > op = 
			Teuchos::rcp(new RestrictorType(
						mgGrids->get(level-1),
						mgGrids->get(level),
						mgGrids->get()->getProcGrid()->getNP()));

		//if(mgGrids->participating(level-1) && print) op->print();
		if(print) {
      auto bla = Pimpact::createOstream("restrictor"+std::to_string(grid->rankST())+".txt");
      op->print(*bla);
      if(grid->rankST()==rankbla)
        op->print();
		}

		std::vector<Pimpact::F> types;
		if("Restriction SF"==op->getLabel())
			types.push_back(Pimpact::F::S);
		else if("Restriction VF"==op->getLabel()) {
			types.push_back(Pimpact::F::U);
			types.push_back(Pimpact::F::V);
			types.push_back(Pimpact::F::W);
		}

		for(auto type=types.begin(); type!=types.end(); ++type) {
			if(2==CS::GridT::sdim && *type==Pimpact::F::W) break;
			if(0==grid->rankST())
				std::cout <<" --- ftype: " <<*type <<" ---\n";

			Teuchos::RCP<Pimpact::ScalarField<typename CS::CGridT> > sol;
			Teuchos::RCP<Pimpact::ScalarField<typename CS::CGridT> > er;

			Pimpact::ScalarField<typename CS::CGridT> fieldf(mgGrids->get(level-1), Pimpact::Owning::Y, *type);
			Pimpact::ScalarField<typename CS::CGridT> fieldc(mgGrids->get(level  ), Pimpact::Owning::Y, *type);

			sol = fieldc.clone();
			er = fieldc.clone();

			// the zero test
			fieldf.init(0.);
			fieldc.init(1.);

			if(mgGrids->participating(level-1)) op->apply(fieldf, fieldc);

			if(mgGrids->participating(level-1))
				TEST_EQUALITY(fieldf.norm()<eps, true);

			if(mgGrids->participating(level)) 
				TEST_EQUALITY(fieldc.norm()<eps, true);

			// the random test
      fieldf.random();

      if(mgGrids->participating(level-1)) op->apply(fieldf, fieldc);

      if(mgGrids->participating(level))
        TEST_INEQUALITY(0., fieldc.norm());

      // the const test
      fieldf.init(1.);
      fieldc.init(0.);
      sol->init(1.);

      if(mgGrids->participating(level-1)) op->apply(fieldf, fieldc);

      if(mgGrids->participating(level)) {

        er->add(1., *sol, -1., fieldc, Pimpact::B::Y);
        ST errInf = er->norm(Pimpact::ENorm::Inf, Pimpact::B::Y);
        if(0==grid->rankST())
          std::cout <<"error Const: " <<errInf <<" ("<<op->getDD() <<")\n";
        //if(i>0)
        TEST_EQUALITY(errInf<eps, true); // boundaries?
        if(print) er->print();
        if(errInf>=eps)
          if(write) er->write(0);
      }

			// the hard test
      for(int dir=1; dir<=CS::GridT::sdim; ++dir) {
      //{int dir=2;

        Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);
        fieldf.initField(type);
        fieldc.init(0.);
        sol->initField(type);

        if(mgGrids->participating(level-1)) op->apply(fieldf, fieldc);

        if(mgGrids->participating(level)) {
          er->add(1., *sol, -1., fieldc, Pimpact::B::Y);
          double errInf = er->norm(Pimpact::ENorm::Inf, Pimpact::B::Y);
          if(0==grid->rankST())
            std::cout <<"error ("<<type <<"): " <<errInf <<" ("<<op->getDD() <<")\n";
          TEST_EQUALITY(errInf<eps, true);
          if(errInf>=eps) {
            //if(print) er->print();
            if(print) fieldc.print();
            if(write) {
              er->write(dir);
              fieldc.write(10*(dir));
              sol->write(100*(dir));
            }
          }
        }
      }
		}
	} // end of for()

} // end of TEUCHOS_UNIT_TEST_TEMPLATE


using ResSF2D = Pimpact::RestrictionSFOp<CGrid2DT>;
using ResVF2D = Pimpact::RestrictionVFOp<CGrid2DT>;
using ResSF3D = Pimpact::RestrictionSFOp<CGrid3DT>;
using ResVF3D = Pimpact::RestrictionVFOp<CGrid3DT>;

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MGTransfers, Restrictor, CSG2D, ResSF2D)
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MGTransfers, Restrictor, CSG2D, ResVF2D)

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MGTransfers, Restrictor, CSG3D, ResSF3D)
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MGTransfers, Restrictor, CSG3D, ResVF3D)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MGTransfers, Interpolator, CS) {

	setParameter(CS::GridT::sdim);

	auto grid = Pimpact::create<typename CS::GridT>(pl);

	auto mgGrids = Pimpact::createMGGrids<CS>(grid, maxGrids);

	std::cout <<"\nrank: " <<grid->rankST() <<"\tnGridLevels: " <<mgGrids->getNGrids() <<"\n";

	auto mgTransfers = Pimpact::createMGTransfers<Pimpact::TransferOp, Pimpact::RestrictionSFOp, Pimpact::InterpolationOp>(mgGrids);

	for(int level=0; level<mgGrids->getNGrids()-1; ++level) {

		Teuchos::RCP<const Pimpact::InterpolationOp<typename CS::CGridT> > op = 
			Teuchos::rcp(new Pimpact::InterpolationOp<typename CS::CGridT>(
						mgGrids->get(level+1),
						mgGrids->get(level),
						mgGrids->get()->getProcGrid()->getNP()));
		if(0==grid->rankST()) {
			std::cout <<"\n\n\t--- level: " <<level <<"---\n";
		}
		if(grid->rankST()==rankbla && print)
			op->print();

		Teuchos::Tuple<Pimpact::F, 4> type =
			Teuchos::tuple(
					Pimpact::F::S,
					Pimpact::F::U,
					Pimpact::F::V,
					Pimpact::F::W);

		for(int i=fs; i<fe; ++i) {
			if(0==grid->rankST())
				std::cout <<"field type: " <<type[i] <<"\n";

			Pimpact::ScalarField<typename CS::CGridT> fieldf(mgGrids->get(level  ), Pimpact::Owning::Y, type[i]);
			Pimpact::ScalarField<typename CS::CGridT> fieldc(mgGrids->get(level+1), Pimpact::Owning::Y, type[i]);
			auto sol = fieldf.clone();
			auto er = fieldf.clone();


			// the zero test
			fieldf.init(1.);
			fieldc.init(0.);

			if(mgGrids->participating(level))
				op->apply(fieldc, fieldf);

			if(mgGrids->participating(level))
				TEST_EQUALITY(eps>fieldf.norm(Pimpact::ENorm::Inf), true);
			if(mgGrids->participating(level+1))
				TEST_EQUALITY(eps>fieldc.norm(Pimpact::ENorm::Inf), true);

			// the random test
			fieldc.random();
			fieldf.init(0.);

			if(mgGrids->participating(level+1))
				TEST_INEQUALITY(0., fieldc.norm());

			if(mgGrids->participating(level))
				op->apply(fieldc, fieldf);

			if(mgGrids->participating(level+1))
				TEST_INEQUALITY(0., fieldc.norm());


			// the stronger init test
			fieldc.init(1.);
			fieldf.init(0.);
			sol->initField(Pimpact::ConstField, 1.);
			er->random();

			if(mgGrids->participating(level))
				op->apply(fieldc, fieldf);


			er->add(1., *sol, -1., fieldf);
			if(mgGrids->participating(level)) {
				ST errInf = er->norm(Pimpact::ENorm::Inf);
				if(0==grid->rankST())
					std::cout <<"error Const: " <<errInf <<"\n";
				TEST_EQUALITY(errInf <eps, true );
				if(errInf>=eps)
					if(write) er->write(0);
			}
			if(er->normLoc(Pimpact::ENorm::Inf)>=eps && grid->rankST()==rankbla && print){
				//			std::cout <<"rank: " <<grid->rankST() <<"\n";
				er->print();
			}


			// --- hardcore test ---
			for(int dir=1; dir<=3; ++dir) {

				Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

				fieldc.initField(type);
				fieldf.init();
				sol->initField(type);
				er->random();

				if(mgGrids->participating(level)){
					op->apply(fieldc, fieldf);

					er->add(1., *sol, -1., fieldf);

					ST errInf = er->norm(Pimpact::ENorm::Inf);
					if(0==grid->rankST())
						std::cout <<"error (" <<Pimpact::toString(type) <<"): " <<errInf <<"\n";
					TEST_EQUALITY(errInf<eps, true );
					if(errInf>=eps) {
						int i = dir;
						if(write) {
							er->write(i);
							fieldf.write(10*i);
							sol->write(100*i);
						}
						if(print) er->print();
					}
				}
			}
		}
	}

}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGTransfers, Interpolator, CSG2D)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGTransfers, Interpolator, CSG3D)




TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MGTransfers, MGTransfersSF, CS) {

	setParameter(CS::GridT::sdim);

	auto grid = Pimpact::create<typename CS::GridT>(pl);

	auto mgGrids = Pimpact::createMGGrids<CS>(grid, maxGrids);

	auto mgTransfers = Pimpact::createMGTransfers<
		Pimpact::TransferOp, Pimpact::RestrictionSFOp, Pimpact::InterpolationOp>(
				mgGrids);

	Pimpact::MGFields<Pimpact::MGGrids<typename CS::GridT, typename CS::CGridT>, Pimpact::ScalarField > x(mgGrids);

	Pimpact::ScalarField<typename CS::GridT>& fieldf = x.get();
	Pimpact::ScalarField<typename CS::CGridT>& fieldc = x.get(-1);

	ST errInf = 0.;

	// interpolation 
	{
		auto sol = fieldf.clone(Pimpact::ECopy::Shallow);
		auto er = fieldf.clone(Pimpact::ECopy::Shallow);

		// the zero test
		fieldf.init(1.);
		fieldc.init();

		mgTransfers->interpolation(x);

		if(mgGrids->participating(-1)) {
			errInf = fieldc.norm(Pimpact::ENorm::Inf);
			TEST_EQUALITY(errInf<eps, true);
		}
		if(mgGrids->participating(0)) {
			errInf = fieldf.norm(Pimpact::ENorm::Inf);
			TEST_EQUALITY(errInf<eps, true);
			if(0==grid->rankST())
				std::cout <<"\ninterpolation error zero: " <<errInf <<"\n";
		}

		// the random test
		fieldc.random();
		fieldf.init(0.);

		if(mgGrids->participating(-1))
			TEST_INEQUALITY(0., fieldc.norm(Pimpact::ENorm::Inf));

		mgTransfers->interpolation(x);

		if(mgGrids->participating(0))
			TEST_INEQUALITY(0., fieldf.norm(Pimpact::ENorm::Inf));

		// the Const test
		fieldc.init(1.);
		fieldf.init();
		sol->init(1.);

		mgTransfers->interpolation(x);

		er->add(1., *sol, -1., fieldf);

		if(mgGrids->participating(0)) {
			errInf = er->norm(Pimpact::ENorm::Inf);
			if(0==grid->rankST())
				std::cout <<"interpolation error Const: " <<errInf <<"\n";
			TEST_EQUALITY(errInf<eps, true );
			if(errInf>=eps) {
				int i = 0;
				if(write) er->write(std::abs(i));
				std::string r = std::to_string(static_cast<long long>(grid->rankST())); // long long needed on brutus(intel)
				if(print) er->print(*Pimpact::createOstream("int_error_c_r"+r+".txt"));
			}
		}


		// the Grad test
		for(int dir=1; dir<=3; ++dir) {

			Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

			fieldc.initField(type);
			fieldf.init();
			sol->initField(type);

			mgTransfers->interpolation(x);

			er->add(1., *sol, -1., fieldf);

			if(mgGrids->participating(0)) {
				errInf = er->norm(Pimpact::ENorm::Inf);
				if(0==grid->rankST())
					std::cout <<"interpolation error (" <<Pimpact::toString(type) <<"): " <<errInf <<"\n";
				TEST_EQUALITY(errInf<eps, true );
				if(errInf>=eps) {
					int i = dir-2;
					if(write) er->write(std::abs(i));
					std::string d = std::to_string(static_cast<long long>(dir-1)); // long long needed on brutus(intel)
					std::string r = std::to_string(static_cast<long long>(grid->rankST())); // long long needed on brutus(intel)
					if(print) er->print(*Pimpact::createOstream("int_error_g"+d+"_r"+r+".txt"));
				}
			}
		}

	}

	// restriction 
	{
		auto sol = fieldc.clone(Pimpact::ECopy::Shallow);
		auto er = fieldc.clone(Pimpact::ECopy::Shallow);


		// the zero test
		fieldf.init();
		fieldc.init(1.);

		mgTransfers->restriction(x);

		if(mgGrids->participating(0))
			TEST_EQUALITY(fieldf.norm()<eps, true);
		if(mgGrids->participating(-1))
			TEST_EQUALITY(fieldc.norm()<eps, true);

		// the random test
		fieldf.random();
		fieldc.init(0.);

		if(mgGrids->participating(0))
			TEST_INEQUALITY(0., fieldf.norm());

		mgTransfers->restriction(x);

		if(mgGrids->participating(-1))
			TEST_INEQUALITY(0., fieldc.norm());

		// the const test
		fieldf.initField(Pimpact::ConstField, 1.);
		fieldc.init();
		sol->initField(Pimpact::ConstField, 1.);

		mgTransfers->restriction(x);

		er->add(1., *sol, -1., fieldc);

		if(mgGrids->participating(-1)) {
			errInf = er->norm(Pimpact::ENorm::Inf);
			TEST_EQUALITY(errInf<eps, true );
		}

		if(0==grid->rankST())
			std::cout <<"\nrestriction error Const: " <<errInf <<"\n";


		// hardcore Grad test
		for(int dir=1; dir<=3; ++ dir) {

			Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

			fieldf.initField(type);
			fieldc.init();

			sol->initField(type);

			mgTransfers->restriction(x);

			er->add(1., *sol, -1., fieldc);
			if(mgGrids->participating(-1)) {

				errInf = er->norm(Pimpact::ENorm::Inf);
				TEST_EQUALITY(errInf<eps, true );
				if(errInf>=eps)
					if(write) er->write(std::abs(1*(dir-2)));
			}
			if(0==grid->rankST())
				std::cout <<"restriction error (" <<Pimpact::toString(type) <<"): " <<errInf <<"\n";
		}
	}
}

//TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGTransfers, MGTransfersSF, CSG2D) 
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGTransfers, MGTransfersSF, CSG3D)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(MGTransfers, MGTransfersVF, CS) {

	setParameter(CS::GridT::sdim);

	auto grid = Pimpact::create<typename CS::GridT>(pl);

	auto mgGrids = Pimpact::createMGGrids<CS>(grid, maxGrids);

	auto mgTransfers = Pimpact::createMGTransfers<
		TransVF, RestrVF, InterVF>(mgGrids);

	Pimpact::MGFields<Pimpact::MGGrids<typename CS::GridT, typename CS::CGridT>, Pimpact::VectorField > x(mgGrids);

	Pimpact::VectorField<typename CS::GridT>& fieldf = x.get();
	Pimpact::VectorField<typename CS::CGridT>& fieldc = x.get(-1);

	// interpolation 
	{
		auto sol = fieldf.clone(Pimpact::ECopy::Shallow);
		auto er = fieldf.clone(Pimpact::ECopy::Shallow);


		// the zero test
		fieldf.init(1.);
		fieldc.init();

		mgTransfers->interpolation(x);

		if(mgGrids->participating(0))
			TEST_EQUALITY(fieldf.norm()<eps, true);
		if(mgGrids->participating(-1))
			TEST_EQUALITY(fieldc.norm()<eps, true);


		// the random test
		fieldc.random();
		fieldf.init(0.);

		if(mgGrids->participating(-1))
			TEST_INEQUALITY(0., fieldc.norm());

		mgTransfers->interpolation(x);

		if(mgGrids->participating(0))
			TEST_INEQUALITY(0., fieldf.norm());

		// the stronger init test
		fieldc(Pimpact::F::U).init(1.);
		fieldc(Pimpact::F::V).init(1.);
		fieldc(Pimpact::F::W).init(1.);
		fieldf.init(0.);
		(*sol)(Pimpact::F::U).init(1.);
		(*sol)(Pimpact::F::V).init(1.);
		(*sol)(Pimpact::F::W).init(1.);

		mgTransfers->interpolation(x);

		er->add(1., *sol, -1., fieldf);

		if(mgGrids->participating(0)) {
			ST rel_error = er->norm(Pimpact::ENorm::Inf);
			if(0==grid->rankST())
				std::cout <<"\nint. error Const: " <<rel_error <<"\n";
			TEST_EQUALITY(rel_error <eps, true );
			if(rel_error>=eps || isnan(rel_error)) {
				if(write) er->write(0);
			}
		}

		// hardcore Grad test
		for(int dir=1; dir<=3; ++ dir) {

			Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

			fieldc(Pimpact::F::U).initField(type);
			fieldc(Pimpact::F::V).initField(type);
			fieldc(Pimpact::F::W).initField(type);
			fieldf.init();

			(*sol)(Pimpact::F::U).initField(type);
			(*sol)(Pimpact::F::V).initField(type);
			(*sol)(Pimpact::F::W).initField(type);

			mgTransfers->interpolation(x);

			er->add(1., *sol, -1., fieldf);

			if(mgGrids->participating(0)) {
				ST rel_error = er->norm(Pimpact::ENorm::Inf);
				if(0==grid->rankST())
					std::cout <<"interpolation error (" <<Pimpact::toString(type) <<"): " <<rel_error <<"\n";
				TEST_EQUALITY(rel_error<eps, true );
				if(rel_error>=eps || isnan(rel_error)) {
					if(write) er->write(1*(dir-2));
					if(write) fieldf.write(10*(dir-2));
					if(write) sol->write(100*(dir-2));
				}
			}
		}

	}
	// restriction 
	{
		auto sol = fieldc.clone(Pimpact::ECopy::Shallow);
		auto er = fieldc.clone(Pimpact::ECopy::Shallow);

		// the zero test
		fieldf.init();
		fieldc.init(1.);

		mgTransfers->restriction(x);

		if(mgGrids->participating(0))
			TEST_EQUALITY(fieldf.norm()<eps, true);
		if(mgGrids->participating(-1))
			TEST_EQUALITY(fieldc.norm()<eps, true);


		// the random test
		fieldf.random();
		fieldc.init(0.);

		if(mgGrids->participating(0))
			TEST_INEQUALITY(0., fieldf.norm());

		mgTransfers->restriction(x);

		if(mgGrids->participating(-1))
			TEST_INEQUALITY(0., fieldc.norm());


		// the stronger init test
		fieldf(Pimpact::F::U).init(1.);
		fieldf(Pimpact::F::V).init(1.);
		fieldf(Pimpact::F::W).init(1.);
		fieldc.init(0.);
		(*sol)(Pimpact::F::U).init(1.);
		(*sol)(Pimpact::F::V).init(1.);
		(*sol)(Pimpact::F::W).init(1.);

		mgTransfers->restriction(x);

		er->add(1., *sol, -1., fieldc);
		if(mgGrids->participating(-1)) {
			ST rel_error = er->norm()/std::sqrt((ST)er->getLength());
			if(0==grid->rankST())
				std::cout <<"res. error Const: " <<rel_error <<"\n";
			TEST_EQUALITY(rel_error<eps, true );
			if(rel_error>eps) {
				if(write) er->write(0);
			}
		}

		// hardcore grad test
		for(int dir=1; dir<=3; ++ dir) {

			Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

			fieldf(Pimpact::F::U).initField(type);
			fieldf(Pimpact::F::V).initField(type);
			fieldf(Pimpact::F::W).initField(type);
			fieldc.init();

			(*sol)(Pimpact::F::U).initField(type);
			(*sol)(Pimpact::F::V).initField(type);
			(*sol)(Pimpact::F::W).initField(type);

			mgTransfers->restriction(x);

			er->add(1., *sol, -1., fieldc);
			if(mgGrids->participating(-1)) {
				ST rel_error = er->norm()/std::sqrt((ST)er->getLength());
				if(0==grid->rankST())
					std::cout <<"restriction grad  error (" <<type <<"): " <<rel_error <<"\n";
				TEST_EQUALITY(rel_error<eps, true );
				if(rel_error>eps) {
					if(write) {
						er->write(dir);
						fieldc.write(10*(dir));
						sol->write(100*(dir));
					}
					if(print) {
						//fieldc.print();
						er->print();
					}
				}
			}
		}
	}
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGTransfers, MGTransfersVF, CSL2D)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGTransfers, MGTransfersVF, CSL3D)

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGTransfers, MGTransfersVF, CSG2D)
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(MGTransfers, MGTransfersVF, CSG3D)



TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(MultiGrid, DivGradOp, CS, MGT) {

	setParameter(CS::GridT::sdim);

	auto grid = Pimpact::create<typename CS::GridT>(pl);

	Pimpact::ScalarField<typename CS::GridT> x(grid);
	Pimpact::ScalarField<typename CS::GridT> b(grid);
	Pimpact::ScalarField<typename CS::GridT> res(grid);
	Pimpact::ScalarField<typename CS::GridT> sol(grid);

	auto xm = Pimpact::wrapMultiField(Teuchos::rcpFromRef(x));
	auto bm = Pimpact::wrapMultiField(Teuchos::rcpFromRef(b));
	//Teuchos::RCP<Pimpact::MultiField<Pimpact::ScalarField<typename CS::GridT> >
		//> xm = Teuchos::rcp(new Pimpact::MultiField<Pimpact::ScalarField<typename
				//CS::GridT> >(grid));
	//Teuchos::RCP<Pimpact::MultiField<Pimpact::ScalarField<typename CS::GridT> >
		//> bm = Teuchos::rcp(new Pimpact::MultiField<Pimpact::ScalarField<typename
				//CS::GridT> >(grid));

	auto op   = Pimpact::create<Pimpact::DivGradOp>(grid);
	auto opO2 = Pimpact::create<Pimpact::DivGradO2Op>(grid);

	auto mgGrids = Pimpact::createMGGrids<CS>(grid, maxGrids);

	auto mgPL = Teuchos::parameterList();

	// MG
	mgPL->set<int>("numCycles", 1);
	mgPL->set<bool>("defect correction", false);
	//mgPL->set<bool>("defect correction", true);
	//mgPL->set<bool>("init zero", false);
	mgPL->set<bool>("init zero", true);

	// Smoother: Line
	//mgPL->sublist("Smoother").set<int>("numIters", 4);
	//mgPL->sublist("Smoother").set<bool>("X", true);
	//mgPL->sublist("Smoother").set<bool>("Y", false);
	//mgPL->sublist("Smoother").set<bool>("Z", false);

	// Smoother: JT

	// Smoother: Chebyshev
	// compute EV
	//ST evMax;
	//ST evMin;

	//Teuchos::RCP<Pimpact::TeuchosEigenvalues<Pimpact::DivGradO2Op<typename CS::GridT> > > ev = 
	//Teuchos::rcp(new Pimpact::TeuchosEigenvalues<Pimpact::DivGradO2Op<typename CS::GridT> >(opO2));
	//ev->computeEV(evMax, evMin);
	//std::cout <<"glob: " <<evMax <<"\t" <<evMin <<"\n";
	//////ev->computeFullEV(evMax, evMin);
	//////std::cout <<"glob: " <<evMax <<"\t" <<evMin <<"\n";

	//mgPL->sublist("Smoother").set<int>("numIters", 8);
	//mgPL->sublist("Smoother").set<ST>("min EV", evMin*1.1);
	//mgPL->sublist("Smoother").set<ST>("max EV", evMin*1.1/30.);

	Teuchos::RCP<MGT> mg =
		Teuchos::rcp(new MGT(mgGrids, op, mgPL));

	auto prec = Pimpact::createMultiOperatorBase(mg);

	std::string solvName = "Pseudoblock GMRES";
	//std::string solvName = "Block GMRES";
	//std::string solvName = "Block CG";
	//std::string solvName = "Pseudoblock CG";
	//std::string solvName = "Pseudoblock Stochastic CG";
	//std::string solvName = "GCRODR";
	//std::string solvName = "RCG";
	//std::string solvName = "MINRES";
	//std::string solvName = "LSQR";
	//std::string solvName = "TFQMR";
	//std::string solvName = "Pseudoblock TFQMR";
	//std::string solvName = "Hybrid Block GMRES";
	//std::string solvName = "PCPG";
	//std::string solvName = "Fixed Point";
	//std::string solvName = "BiCGStab";

	Teuchos::RCP<Teuchos::ParameterList > param = Teuchos::parameterList();// = Pimpact::createLinSolverParameter(solvName, 1.e-6);
	param->set("Output Style", Belos::Brief);
	//param->set("Output Stream", Teuchos::rcpFromRef(&std::cout));
	param->set("Output Frequency", 1);
	param->set("Verbosity",			        
			Belos::Errors +
			Belos::Warnings +
			Belos::IterationDetails +
			Belos::OrthoDetails +
			Belos::FinalSummary +
			Belos::TimingDetails +
			Belos::StatusTestDetails +
			Belos::Debug);
	param->set("Maximum Iterations", nMax);
	//param->set("Flexible Gmres", false);

	auto bop = Pimpact::createMultiOperatorBase(op);

	auto linprob = Pimpact::createLinearProblem<Pimpact::MultiField<Pimpact::ScalarField<typename CS::GridT> > >(bop, xm, bm, param, solvName);
	if(solvName!="Fixed Point") 
		linprob->setRightPrec(prec);
	//if(solvName!="BiCGStab") 
	if(solvName=="Fixed Point") 
		linprob->setLeftPrec(prec);

	// --- zero rhs test ---
	//b.init(0.);
	//x.random();

	//ST error0 = x.norm();
	//ST errorp = error0;

	////opO2->apply(x, res);
	//op->apply(x, res);
	//ST res0 = res.norm();
	//ST resP = res0;

	//if(grid()->rankST()==0) {
	//std::cout <<"\n\n\t\t\t--- zero rhs test ---\n";
	//std::cout <<"\tresidual:\trate:\t\t\terror:\t\trate: \n";
	//std::cout << std::scientific;
	//std::cout <<"\t" <<1. <<"\t\t\t\t"  <<1.  <<"\n";
	//}

	//for(int i=0; i<nMax; ++i) {
	//mg->apply(b, x);
	//x.level();

	//ST error = x.norm();

	////opO2->apply(x, res);
	//op->apply(x, res);
	//ST residual = res.norm();

	//if(grid()->rankST()==0)
	//std::cout <<"\t" <<residual/res0 <<"\t" <<residual/resP <<"\t\t" <<error/error0 <<"\t" << error/errorp <<"\n";

	////if(error>= errorp)
	////break;
	////else
	//errorp = error;
	//resP = residual;
	//}

	////x.print();
	//x.write();
	//TEST_EQUALITY(x.norm()/std::sqrt(static_cast<ST>(x.getLength()))<1.e-3, true);

	//bm->init(0.);
	//xm->random();
	//linprob->solve(xm, bm);

	// --- grad test ---
	auto e = x.clone();
	for(int dir=1; dir<=1; ++dir) {

		Pimpact::EScalarField type = static_cast<Pimpact::EScalarField>(dir);

		x.init(0);
		sol.initField(type);
		sol.level();

		// construct RHS
		opO2->apply(sol, b);
		if(print) b.print();

		// residual
		opO2->apply(x, res);
		res.add(-1., b, 1., res);
		ST res0 = res.norm();
		ST resP = res0;

		// error
		res.add(1., sol, -1., x);
		ST err0 = res.norm();
		ST errP = err0;

		// output
		if(grid()->rankST()==0) {
			std::cout <<"\n\n\t\t\t--- " <<Pimpact::toString(type) <<" test ---\n";
			std::cout <<"\tresidual:\trate:\t\t\terror:\t\trate:\n";
			std::cout << std::scientific;
			std::cout <<"\t"  <<1.<<"\t\t\t\t" <<1.  <<"\n";
		}
		for(int i=0; i<1; ++i) {
			// mg cycle
			mg->apply(b, x);
			//mg->applyVFull(b, x);
			x.level();

			// residual
			opO2->apply(x, res);
			res.add(-1, b, 1., res);
			ST residual = res.norm();
			res.add(1., sol, -1., x);
			ST err = res.norm();

			if(grid()->rankST()==0)
				std::cout <<"\t" <<residual/res0 <<"\t" << residual/resP <<"\t\t" <<err/err0 <<"\t" << err/errP  <<"\n";
			resP = residual;
			errP = err;
		}
		TEST_EQUALITY(res.norm()/std::sqrt(static_cast<ST>(res.getLength()))<1.e-3, true);

		op->apply(sol, b);
		if(print) b.print();
		xm->init(0.);
		//xm->random();
		linprob->solve(xm, bm);
		if(write) xm->write();
	}
}


//TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MultiGrid, DivGradOp, CSG2D, DGJMGT2D)
//TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MultiGrid, DivGradOp, CSG2D, DGLMGT2D)
//TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MultiGrid, DivGradOp, CSG2D, DGCMGT2D)

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MultiGrid, DivGradOp, CSG3D, DGJMGT3D)
//TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MultiGrid, DivGradOp, CSG3D, DGLMGT3D)
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MultiGrid, DivGradOp, CSG3D, DGCMGT3D)



TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(MultiGrid, ConvDiffOp, CS, MGT) {

	setParameter(CS::GridT::sdim);

	auto grid = Pimpact::create<typename CS::GridT>(pl);

	auto mgGrids = Pimpact::createMGGrids<CS>(grid, maxGrids);
	if(-1==print)  mgGrids->print();

	auto mgPL = Teuchos::parameterList();
	mgPL->set<int>("numCycles", 1);
	mgPL->set<bool>("defect correction", false);
	//mgPL->set<bool>("defect correction", true);
	mgPL->set<bool>("init zero", true);
	//mgPL->set<bool>("init zero", false);

	mgPL->sublist("Smoother").set<ST>("omega", 0.5);
	mgPL->sublist("Smoother").set("numIters", 10);
	//mgPL->sublist("Smoother").set("Ordering", 1);
	//mgPL->sublist("Smoother").set<short int>("dir X", -1);
	//mgPL->sublist("Smoother").set<short int>("dir Y", -1);
	//mgPL->sublist("Smoother").set<short int>("dir Z", -1);

	mgPL->sublist("Coarse Grid Solver").set<std::string>("Solver name", "GMRES");

	//	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<Teuchos::RCP<std::ostream> >("Output Stream", Teuchos::rcp(&std::cout, false));
		mgPL->sublist("Coarse Grid Solver").sublist("Solver").set("Verbosity",
				Belos::Errors +
				Belos::Warnings +
				////Belos::IterationDetails +
				////Belos::OrthoDetails +
				Belos::FinalSummary +
				////Belos::TimingDetails +
				////Belos::StatusTestDetails +
				Belos::Debug);
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<std::string>("Timer Label", "Coarse Grid Solver");
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set<ST>("Convergence Tolerance"
			, 1.0e-6);
	mgPL->sublist("Coarse Grid Solver").sublist("Solver").set("Maximum Iterations", 100);

	auto op = Pimpact::create<ConvDiffOpT>(grid);

	Teuchos::RCP<MGT> mg = Teuchos::rcp(new MGT(mgGrids, op, mgPL));


	Pimpact::VectorField<typename CS::GridT> x   (grid);
	Pimpact::VectorField<typename CS::GridT> b   (grid);
	Pimpact::VectorField<typename CS::GridT> temp(grid);
	Pimpact::VectorField<typename CS::GridT> sol( grid);

	//{
		//Pimpact::VectorField<typename CS::GridT> wind(grid);
		//auto windfunc = [&pi2](ST x) -> ST { return(-std::cos(pi2*x/2)); };
		//wind(Pimpact::F::U).initFromFunction(
				//[&windfunc](ST x, ST y, ST z) ->ST { return(windfunc(x)); });
		//wind(Pimpact::F::V).initFromFunction(
				//[&windfunc](ST x, ST y, ST z) ->ST { return(windfunc(y)); });
		//if(3==CS::GridT::sdim)
			//wind(Pimpact::F::W).initFromFunction(
					//[&windfunc](ST x, ST y, ST z) ->ST { return(windfunc(z)); });

		//op->assignField(wind);
		//mg->assignField(wind);
	//}

	std::ofstream ofs;
	if(grid()->rankST()==0)
		ofs.open("MG2.txt", std::ofstream::out);

	// 
	sol(Pimpact::F::U).initField(Pimpact::Grad2D_inX);
	sol(Pimpact::F::V).initField(Pimpact::Grad2D_inY);
	if(3==CS::GridT::sdim) sol(Pimpact::F::W).initField(Pimpact::Grad2D_inZ);

	op->apply(sol, b);

	if(write) b.write(1);

	x.init();

	temp.add(-1., x, 1., sol);
	ST res = temp.norm();
	ST res_0 = res;
	ST res_p = res;

	if(grid()->rankST()==0) {
		std::cout << std::scientific;
		std::cout <<"\n\n\tresidual:\trate: \n";
		std::cout <<"\t"  <<1.  <<"\n";
		ofs <<res <<"\n";
	}

	for(int i=0; i<nIter; ++i) {
		mg->apply(b, x);
		//mg->applyVFull(b, x);
		if(write) x.write(i+10);

		temp.add(-1., x, 1., sol);
		ST res = temp.norm();

		if(grid()->rankST()==0) {
			std::cout <<"\t" <<res/res_0 <<"\t" << res/res_p <<"\n";
			ofs <<res <<"\n";
		}
		res_p = res;
	}

	if(grid()->rankST()==0) {
		std::cout <<"\n";
	}

	TEST_EQUALITY(temp.norm()<0.5, true);

	if(write) x.write(2);

	if(grid()->rankST()==0)
		ofs.close();


	// MG as preconditioner
	std::string solvName = "Pseudoblock GMRES";
	//std::string solvName = "Block GMRES";
	//std::string solvName = "Block CG";
	//std::string solvName = "Pseudoblock CG";
	//std::string solvName = "Pseudoblock Stochastic CG";
	//std::string solvName = "GCRODR";
	//std::string solvName = "RCG";
	//std::string solvName = "MINRES";
	//std::string solvName = "LSQR";
	//std::string solvName = "TFQMR";
	//std::string solvName = "Pseudoblock TFQMR";
	//std::string solvName = "Hybrid Block GMRES";
	//std::string solvName = "PCPG";
	//std::string solvName = "Fixed Point";
	//std::string solvName = "BiCGStab";

	Teuchos::RCP<Teuchos::ParameterList > param = Teuchos::parameterList();// = Pimpact::createLinSolverParameter(solvName, 1.e-6);
	param->set("Output Style", Belos::Brief);
	//param->set("Output Stream", Teuchos::rcpFromRef(&std::cout));
	param->set("Output Frequency", 1);
	param->set("Verbosity",			        
			Belos::Errors +
			Belos::Warnings +
			Belos::IterationDetails +
			Belos::OrthoDetails +
			Belos::FinalSummary +
			Belos::TimingDetails +
			Belos::StatusTestDetails +
			Belos::Debug);
	param->set("Maximum Iterations", nMax);
	//param->set("Flexible Gmres", false);
	//
	auto xm = Pimpact::wrapMultiField(Teuchos::rcpFromRef(x));
	auto bm = Pimpact::wrapMultiField(Teuchos::rcpFromRef(b));

	auto bop = Pimpact::createMultiOperatorBase(op);
	auto prec = Pimpact::createMultiOperatorBase(mg);

	auto linprob = Pimpact::createLinearProblem<Pimpact::MultiField<Pimpact::VectorField<typename CS::GridT> > >(bop, xm, bm, param, solvName);
	if(solvName!="Fixed Point") 
		linprob->setRightPrec(prec);
	//if(solvName!="BiCGStab") 
	if(solvName=="Fixed Point") 
		linprob->setLeftPrec(prec);

	bm->init(0.);
	xm->random();
	linprob->solve(xm, bm);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MultiGrid, ConvDiffOp, CSG2D, CDJ2D)
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MultiGrid, ConvDiffOp, CSG3D, CDJ3D)

TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MultiGrid, ConvDiffOp, CSG2D, CDSOR2D)
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MultiGrid, ConvDiffOp, CSG3D, CDSOR3D)
//TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(MultiGrid, ConvDiffOp, CCSG3D, CCDSOR3D)



} // end of namespace
