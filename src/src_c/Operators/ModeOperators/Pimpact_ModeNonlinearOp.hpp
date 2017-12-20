#pragma once
#ifndef PIMPACT_MODENONLINEAROP_HPP
#define PIMPACT_MODENONLINEAROP_HPP


#include "Pimpact_ModeField.hpp"




namespace Pimpact {



/// \ingroup ModeOperator
template<class OpT>
class ModeNonlinearOp {

public:

  using SpaceT = typename OpT::SpaceT;

  using InnerOpT = OpT;

  using DomainFieldT = ModeField<typename OpT::DomainFieldT >;
  using RangeFieldT = ModeField<typename OpT::RangeFieldT  >;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  Scalar mulI_;
  Scalar mulC_;
  Scalar mulL_;

  Teuchos::RCP<OpT> op_;

public:

  ModeNonlinearOp( const Teuchos::RCP<const SpaceT>& space ):
    mulI_(0.),
    mulC_(1.),
    mulL_( 1./space->getDomainSize()->getRe() ),
    op_( create<OpT>( space ) ) { };

  ModeNonlinearOp( const Teuchos::RCP<OpT>& op ):
    mulI_(0.),
    mulC_(1.),
    mulL_( 1./op->space()->getDomainSize()->getRe() ),
    op_( op ) { };


  void apply(const DomainFieldT& x, RangeFieldT& y ) const {

    // set paramters
    auto pl = Teuchos::parameterList();
    pl->set<Scalar>( "mulI", 0.    );
    pl->set<Scalar>( "mulC", mulC_ );
    pl->set<Scalar>( "mulL", mulL_ );
    op_->setParameter( pl );

    //std::cout << *pl;
    //typename OpT::RangeFieldT temp( space() );

    const B wnB = B::N;

    for( F m=F::U; m<SpaceT::sdim; ++m ) {
      x.getCField()(m).exchange();
      x.getSField()(m).exchange();

      //std::cout << "wind0: "  << op_->getConvField()->get()[static_cast<int>(m)][0].norm() << "\n";
      //std::cout << "wind1: "  << op_->getConvField()->get()[static_cast<int>(m)][1].norm() << "\n";
      //std::cout << "wind2: "  << op_->getConvField()->get()[static_cast<int>(m)][2].norm() << "\n";

      for( Ordinal k=space()->si(m,Z,wnB); k<=space()->ei(m,Z,wnB); ++k )
        for( Ordinal j=space()->si(m,Y,wnB); j<=space()->ei(m,Y,wnB); ++j )
          for( Ordinal i=space()->si(m,X,wnB); i<=space()->ei(m,X,wnB); ++i ) {
            if( 3==SpaceT::sdim ) {
              y.getCField()(m)(i,j,k) = 
                mulI_*x.getSField()(m)(i,j,k)
                +mulC_*op_->getSOp()->getConvSOp()->innerStenc3D(
                    op_->getConvField(m)[0](i,j,k),
                    op_->getConvField(m)[1](i,j,k),
                    op_->getConvField(m)[2](i,j,k),
                    x.getCField()(m), i, j, k)
                -mulL_*op_->getSOp()->getHelmOp()->innerStenc3D(
                    x.getCField()(m), m, i, j, k) ;
              y.getSField()(m)(i,j,k) = 
                -mulI_*x.getCField()(m)(i,j,k)
                +mulC_*op_->getSOp()->getConvSOp()->innerStenc3D(
                    op_->getConvField(m)[0](i,j,k),
                    op_->getConvField(m)[1](i,j,k),
                    op_->getConvField(m)[2](i,j,k),
                    x.getSField()(m), i, j, k)
                -mulL_*op_->getSOp()->getHelmOp()->innerStenc3D(
                    x.getSField()(m), m, i, j, k) ;
            }
            else {
              y.getCField()(m)(i,j,k) = 
                mulI_*x.getSField()(m)(i,j,k)
                +mulC_*op_->getSOp()->getConvSOp()->innerStenc2D(
                    op_->getConvField(m)[0](i,j,k),
                    op_->getConvField(m)[1](i,j,k),
                    x.getCField()(m), i, j, k)
                -mulL_*op_->getSOp()->getHelmOp()->innerStenc2D(
                    x.getCField()(m), m, i, j, k) ;
              y.getSField()(m)(i,j,k) = 
                -mulI_*x.getCField()(m)(i,j,k)
                +mulC_*op_->getSOp()->getConvSOp()->innerStenc2D(
                    op_->getConvField(m)[0](i,j,k),
                    op_->getConvField(m)[1](i,j,k),
                    x.getSField()(m), i, j, k)
                -mulL_*op_->getSOp()->getHelmOp()->innerStenc2D(
                    x.getSField()(m), m, i, j, k) ;
            }
          }

      // boundaries
      op_->getSOp()->getHelmOp()->applyBC( x.getCField()(m), y.getCField()(m) );
      op_->getSOp()->getHelmOp()->applyBC( x.getSField()(m), y.getSField()(m) );

      y.getCField()(m).changed();
      y.getSField()(m).changed();
    }
    //y = x;
  }

  void computeResidual( const RangeFieldT& b, const DomainFieldT& x, RangeFieldT& res ) const {
    apply( x, res );
    res.add( 1., b, -1., res );
  }

  void assignField( const DomainFieldT& mv ) {
    op_->assignField( mv.getCField() );
  };


  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return op_->space();
  };


  void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {
    if( para->name()!="Linear Solver" ) {
      mulI_ = para->get<Scalar>( "mulI" );
      mulC_ = para->get<Scalar>( "mulC" );
      mulL_ = para->get<Scalar>( "mulL" );
    }
  }


  Teuchos::RCP<OpT> getInnerOpPtr() {
    return op_;
  }


  bool hasApplyTranspose() const {
    return false;
  }


  const std::string getLabel() const {
    return "ModeNonlinearOp_"+op_->getLabel();
  };


  void print( std::ostream& out=std::cout ) const {
    out << getLabel() << ":\n";
    op_->print( out );
  }


}; // end of class ModeNonlinearOp




/// \relates ModeNonlinearOp
template< class OpT >
Teuchos::RCP< ModeNonlinearOp<OpT> >
createModeNonlinearOp( const Teuchos::RCP<OpT>& op ) {

  return Teuchos::rcp( new ModeNonlinearOp<OpT>( op ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MODENONLINEAROP_HPP
