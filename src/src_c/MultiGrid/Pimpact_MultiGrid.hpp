#pragma once
#ifndef PIMPACT_MULTIGRID_HPP
#define PIMPACT_MULTIGRID_HPP


#include "Pimpact_MGFields.hpp"
#include "Pimpact_MGOperators.hpp"
#include "Pimpact_MGTransfers.hpp"
#include "Pimpact_MGSmoothers.hpp"
#include "Pimpact_MGSpaces.hpp"




/// \defgroup MG MultiGrid
///
/// Multi Grid



namespace Pimpact {




/// \brief basic multi grid calls
/// 
/// \tparam MGSpacesT grid hierarchy
/// \tparam FieldT field type
/// \tparam TransT transfer operator type
/// \tparam RestrT restriction operator type
/// \tparam InterT interpolation operator type
/// \tparam FOperatorT high order operator type
/// \tparam COperatorT low order operator type
/// \tparam SmootherT smoother type
/// \tparam CGST coarse grid solver type
/// 
/// \ingroup MG
template<
	class MGSpacesT,
	template<class> class FieldT,
	template<class,class> class TransT,
	template<class> class RestrT,
	template<class> class InterT,
	template<class> class FOperatorT,
	template<class> class COperatorT,
	template<class> class SmootherT,
	template<class> class CGST >
class MultiGrid {

	using FSpaceT = typename MGSpacesT::FSpaceT;
	using CSpaceT = typename MGSpacesT::CSpaceT;


	using MGTransfersT = MGTransfers<MGSpacesT,TransT,RestrT,InterT>;

	using MGFieldsT = MGFields<MGSpacesT,FieldT>;

	using MGOperatorsT = MGOperators<MGSpacesT,FOperatorT,COperatorT>;

	using MGSmoothersT = MGSmoothers<MGOperatorsT, SmootherT>;

public:

	using SpaceT =  FSpaceT;

	using DomainFieldT = FieldT<FSpaceT>;
	using RangeFieldT = FieldT<FSpaceT>;

	using CGridSolverT = CGST< COperatorT<CSpaceT> >;

protected:

	bool defectCorrection_;
	bool initZero_;
	int numCycles_;

	const Teuchos::RCP<const MGSpacesT> mgSpaces_;

	Teuchos::RCP<const MGTransfersT> mgTrans_;

	Teuchos::RCP<const MGOperatorsT> mgOps_;

	Teuchos::RCP<const MGSmoothersT> mgSms_;

	Teuchos::RCP<CGridSolverT> cGridSolver_;

public:

	/// \brief constructor
	///
	/// \param mgSpaces
	/// \param pl  Parameter list of options for the multi grid solver.
	///   These are the options accepted by the solver manager:
	///   - "defect correction" - a \c bool specifying if defect correction is
	///     applied before cycling. Default: true  /
	///   - "init zero" - a \c bool specifying if defect correction is
	///     applied before cylcing. Default: false  /
	///   - "numCycles" - a \c int number of cycles. Default:
	///     FSpaceT::dimNC - CSpaceT::dimNC+1)  /
	///   - "Smoother" - a \c sublist for smoothers
	///   - "Coarse Grid Solver" - a \c sublist for coarse grid solver
	MultiGrid(
			const Teuchos::RCP<const MGSpacesT>& mgSpaces,
			const Teuchos::RCP<Teuchos::ParameterList>& pl ):
		defectCorrection_( pl->get<bool>("defect correction", true ) ),
		initZero_( pl->get<bool>("init zero", false ) ),
		numCycles_( pl->get<int>("numCycles",FSpaceT::dimNC-CSpaceT::dimNC+1) ),
		mgSpaces_(mgSpaces),
		mgTrans_( createMGTransfers<TransT,RestrT,InterT>(mgSpaces) ),
		mgOps_(   createMGOperators<FOperatorT,COperatorT>(mgSpaces) ),
		mgSms_(   createMGSmoothers<SmootherT>( mgOps_, Teuchos::rcpFromRef(pl->sublist("Smoother")) ) ),
		cGridSolver_(
				mgSpaces_->participating(-1)?create<CGridSolverT>( mgOps_->get(-1), Teuchos::rcpFromRef(pl->sublist("Coarse Grid Solver")) ):Teuchos::null ) {}



	/// \brief solves \f$ L y = x \f$
	/// defect correction\f$ \hat{L}u_{k+1} = f-L u_k +\hat{L}u_k \f$ and V-cylce for solving with \f$\hat{L}\f$
	/// \todo extract smooth/restrict/interpolate method???
	/// \todo template cycle method
	void apply( const DomainFieldT& x0, RangeFieldT& y ) const {

		Teuchos::RCP<MGFieldsT> x    = createMGFields<FieldT>( mgSpaces_ );
		Teuchos::RCP<MGFieldsT> temp = createMGFields<FieldT>( mgSpaces_ );
		Teuchos::RCP<MGFieldsT> b    = createMGFields<FieldT>( mgSpaces_ );

		for( int j=0; j<numCycles_; ++j ) {

			if( defectCorrection_ && !initZero_ ) {
				// defect correction rhs \hat{f}= b = x - L y
				mgOps_->get()->computeResidual( x0, y, *b->get() );

				// transfer init y and \hat{f} to coarsest coarse
				mgTrans_->getTransferOp()->apply( y, *x->get(0) );
				mgTrans_->getTransferOp()->apply( *b->get(), *b->get(0) );

				// residual temp = \hat(L) y
				mgOps_->get(0)->apply( *x->get(0), *temp->get(0) );
				// b = x - L y +\hat{L} y
				b->get(0)->add( 1., *b->get(0), 1, *temp->get(0) );

			}
			else {
				// === no defect correction
				if( initZero_ && 0==j )
					x->get(0)->init( 0. );
				else
					mgTrans_->getTransferOp()->apply( y, *x->get(0) );
				mgTrans_->getTransferOp()->apply( x0, *b->get(0) );
			}

			int i;
			for( i=0; i<mgSpaces_->getNGrids()-1; ++i ) {
				if( i>0 ) x->get(i)->init(0.); // necessary? for DivGradOp yes

				if( mgSpaces_->participating(i) ) {
					mgSms_->get(i)->apply( *b->get(i), *x->get(i) );
					mgOps_->get(i)->computeResidual( *b->get(i), *x->get(i), *temp->get(i) );
					mgTrans_->getRestrictionOp(i)->apply( *temp->get(i), *b->get(i+1) );
				}
			}

			// coarse grid solution
			i = -1;
			if( mgSpaces_->participating(i) ) {
				x->get(i)->init(0.);
				try{
					cGridSolver_->apply( *b->get(i), *x->get(i) );
				}
				catch( std::logic_error& e ) {
					std::cout << "error in MG on coarse grid:\n";
					cGridSolver_->print();
					b->get(i)->print();
					b->get(i)->write(111);
					x->get(i)->write(222);
					throw( e );
				}
				//x->get(i)->level();
			}

			for( i=-2; i>=-mgSpaces_->getNGrids(); --i ) {
				// interpolate/correct/smooth
				if( mgSpaces_->participating(i) ) {
					mgTrans_->getInterpolationOp(i)->apply( *x->get(i+1), *temp->get(i) );
					x->get(i)->add( 1., *temp->get(i), 1., *x->get(i) );
					//
					mgSms_->get( i )->apply( *b->get(i), *x->get(i) );
				}
			}

			// use temp as stopping cirterion
			mgTrans_->getTransferOp()->apply( *x->get(0), y );
		}
	}


	void assignField( const DomainFieldT& mv ) {

		Teuchos::RCP<MGFieldsT> temp = createMGFields<FieldT>( mgSpaces_ );

		mgOps_->get()->assignField( mv );
		mgTrans_->getTransferOp()->apply( mv, *temp->get(0) );

		for( int i=0; i<mgSpaces_->getNGrids()-1; ++i )  {
			if( mgSpaces_->participating(i) ) {
				mgOps_->get(i)->assignField( *temp->get(i) );
				mgTrans_->getRestrictionOp(i)->apply( *temp->get(i), *temp->get(i+1) );
			}
		}

		if( mgSpaces_->participating(-1) )
			mgOps_->get(-1)->assignField( *temp->get(-1) );
	};

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(mgSpaces_->get()); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		mgOps_->setParameter( para );
	}

	bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "MultiGrid( "+mgOps_->get()->getLabel()+" ) " ); };


	void print( std::ostream& out=std::cout ) const {
		out << "--- " << getLabel() << " ---\n";
		out << "#grids: " << mgSpaces_->getNGrids() << " numCycles: "<<numCycles_<< "init zero: "<<initZero_ << "\n";
		out << "FOperator: " << mgOps_->get()->getLabel() << " d" << FSpaceT::dimNC << "\n";
		if( mgSpaces_->participating(0) )
			out << "COperator: " << mgOps_->get(0)->getLabel()  << " d" << CSpaceT::dimNC << "\n";
		if( mgSpaces_->participating(0) )
			out << "Smoother: " << mgSms_->get(0)->getLabel() << "\n";
		if( mgSpaces_->participating(-1) )
			out << "Coarse Grid Solver: " << cGridSolver_->getLabel() << "\n";
	}


}; // end of class MultiGrid



/// \relates MultiGrid
template<
	template<class> class FieldT,
	template<class,class> class TransT,
	template<class> class RestrT,
	template<class> class InterT,
	template<class> class FOperatorT,
	template<class> class COperatorT,
	template<class> class SmootherT,
	template<class> class CGridSolverT,
	class MGSpacesT >
Teuchos::RCP< MultiGrid<MGSpacesT,FieldT,TransT,RestrT,InterT,FOperatorT,COperatorT,SmootherT,CGridSolverT> >
createMultiGrid(
		const Teuchos::RCP<const MGSpacesT>& mgSpaces,
		const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList() ) {

	return(
			Teuchos::rcp(
				new MultiGrid<MGSpacesT,FieldT,TransT,RestrT,InterT,FOperatorT,COperatorT,SmootherT,CGridSolverT>(
					mgSpaces, pl )
				)
			);
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIGRID_HPP
