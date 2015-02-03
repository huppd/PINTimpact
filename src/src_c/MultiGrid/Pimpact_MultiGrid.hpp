#pragma once
#ifndef PIMPACT_MULTIGRID_HPP
#define PIMPACT_MULTIGRID_HPP




#include "Pimpact_MGSpaces.hpp"
#include "Pimpact_MGFields.hpp"
#include "Pimpact_MGOperators.hpp"
#include "Pimpact_MGTransfers.hpp"
#include "Pimpact_MGSmoothers.hpp"


/// \defgroup MG MultiGrid
///
/// Multi Grid




namespace Pimpact {


/// \ingroup MG
/// \todo think about templeting everything MGSpaces,  MGOperators,...
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

  typedef typename MGSpacesT::FSpaceT FSpaceT;
  typedef typename MGSpacesT::CSpaceT CSpaceT;


  typedef MGTransfers<MGSpacesT,TransT,RestrT,InterT> MGTransfersT;

  typedef MGFields<MGSpacesT,FieldT> MGFieldsT;

  typedef MGOperators<MGSpacesT,FOperatorT,COperatorT>   MGOperatorsT;

  typedef MGSmoothers<MGOperatorsT, SmootherT> MGSmoothersT;

public:

  typedef FSpaceT SpaceT;

  typedef FieldT<FSpaceT>  DomainFieldT;
  typedef FieldT<FSpaceT>  RangeFieldT;

  typedef CGST< COperatorT<CSpaceT> > CGridSolverT;

protected:

  Teuchos::RCP<const MGSpacesT> mgSpaces_;

  Teuchos::RCP<const MGTransfersT> mgTrans_;

  Teuchos::RCP<const MGOperatorsT> mgOps_;

  Teuchos::RCP<const MGSmoothersT> mgSms_;

  Teuchos::RCP<MGFieldsT> x_;
  Teuchos::RCP<MGFieldsT> temp_;
  Teuchos::RCP<MGFieldsT> b_;

  Teuchos::RCP<CGridSolverT> cGridSolver_;

public:

  MultiGrid(
      const Teuchos::RCP<const MGSpacesT>& mgSpaces,
      EField type = EField::S ):
        mgSpaces_(mgSpaces),
        mgTrans_( createMGTransfers<TransT,RestrT,InterT>(mgSpaces) ),
        mgOps_( Pimpact::createMGOperators<FOperatorT,COperatorT>(mgSpaces) ),
        mgSms_( Pimpact::createMGSmoothers<SmootherT>(mgOps_) ),
        x_( createMGFields<FieldT>(mgSpaces, type ) ),
        temp_( createMGFields<FieldT>(mgSpaces, type ) ),
        b_( createMGFields<FieldT>(mgSpaces, type ) ),
        cGridSolver_( create<CGridSolverT>( mgOps_->get(-1) ) ) {}



  /// \brief solves \f$ L y = x \f$
  /// defect correction\f$ \hat{L}u_{k+1} = f-L u_k +\hat{L}u_k \f$ and V-cylce for solving with \f$\hat{L}\f$
  /// \todo extract smooth/restrict/interpolate method
  /// \todo template cycle method
  void apply( const DomainFieldT& x, RangeFieldT& y ) const {

//    int out =0;
    for( int j=0; j<1; ++j ) {
      // defect correction rhs \hat{f}= b = x - L y
      mgOps_->get()->apply( y, *b_->get() );
      b_->get()->add( 1., x, -1., *b_->get() );

      // transfer init y and \hat{f} to coarsest coarse
      mgTrans_->getTransferOp()->apply( y, *x_->get(0) );
      mgTrans_->getTransferOp()->apply( *b_->get(), *b_->get(0) );

      // residual temp = \hat(L) y
      mgOps_->get(0)->apply( *x_->get(0), *temp_->get(0) );
      // b = x - L y +\hat{L} y
      b_->get(0)->add( 1., *b_->get(0), 1, *temp_->get(0) );

			//b_->get(0)->write(78);
			//x_->get(0)->write(79);
      // smooth and restrict defect( todo extract this as method )
      int i;
      for( i=0; i<mgSpaces_->getNGrids()-1; ++i ) {
				if(i>0) x_->get(i)->init(0.); // necessary? for DivGradOp yes
				//b_->get(i)->write(98);
        mgSms_->get(i)->apply( *b_->get(i), *x_->get(i) );
				//x_->get(i)->write(98);
        mgOps_->get(i)->apply( *x_->get(i), *temp_->get(i) );
        temp_->get(i)->add( -1., *temp_->get(i), 1., *b_->get(i) );
				//temp_->get(i)->write(98);
        mgTrans_->getRestrictionOp(i)->apply( *temp_->get(i), *b_->get(i+1) );
				//b_->get(i+1)->write(98);
      }
			//b_->get(1)->write(88);
			//x_->get(1)->write(89);

      // coarse grid solution
      i = -1;
      /// \todo add level for singular stuff
      b_->get(i)->level();
      x_->get(i)->init(0.);
			//for( int j=0; j<1; ++j )
		  cGridSolver_->apply( *b_->get(i), *x_->get(i) );
			 //b_->get(i)->write(98);
			 //x_->get(i)->write(98);

      for( i=-2; i>=-mgSpaces_->getNGrids(); --i ) {
        // interpolate/correct/smooth
        mgTrans_->getInterpolationOp(i)->apply( *x_->get(i+1), *temp_->get(i) );
				//temp_->get(i)->write(98);
        x_->get(i)->add( 1., *temp_->get(i), 1., *x_->get(i) );
				//b_->get(i)->write(97);
				//x_->get(i)->write( 98 );
				mgSms_->get( i )->apply( *b_->get(i), *x_->get(i) );
				//x_->get(i)->write(99);

      }
      mgTrans_->getTransferOp()->apply( *x_->get(0), y );
      //    y.level();// only laplace

    }

  }


  void assignField( const DomainFieldT& mv ) {

    mgOps_->get()->assignField( mv );
    mgTrans_->getTransferOp()->apply( mv, *temp_->get(0) );
    mgOps_->get(0)->assignField( *temp_->get(0) );

    for( int i=0; i<mgSpaces_->getNGrids()-1; ++i )  {
      mgTrans_->getRestrictionOp(i)->apply( *temp_->get(i), *temp_->get(i+1) );
      mgOps_->get(i+1)->assignField( *temp_->get(i+1) );
    }

  };


  bool hasApplyTranspose() const { return( false ); }


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
    EField type = EField::S  ) {

  return(
      Teuchos::rcp(
          new MultiGrid<MGSpacesT,FieldT,TransT,RestrT,InterT,FOperatorT,COperatorT,SmootherT,CGridSolverT>(
              mgSpaces, type)
      )
  );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIGRID_HPP
