#pragma once
#ifndef PIMPACT_ADD2OP_HPP
#define PIMPACT_ADD2OP_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Space.hpp" // just for create<>




namespace Pimpact {


/// \brief combines two operators.
///
/// Both operators are applied and the result is added.
/// the \c DomainFieldT  \c RangeFieldT has to be equal for both \c OP1 \c OP2.
/// \ingroup Operator
template< class OP1, class OP2 >
class Add2Op {

public:

  typedef typename OP1::DomainFieldT DomainFieldT;
  typedef typename OP2::RangeFieldT  RangeFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

protected:

  Teuchos::RCP<OP1> op1_;
  Teuchos::RCP<OP2> op2_;
  Teuchos::RCP<typename OP1::RangeFieldT> temp_;

public:

  Add2Op(
      const Teuchos::RCP<OP1>& op1,
      const Teuchos::RCP<OP2>& op2 ):
        op1_(op1),
        op2_(op2),
        temp_( create<typename OP1::RangeFieldT>(op1_->space()) )
        {};


  void apply(const DomainFieldT& x, RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS ) const {
    op1_->apply( x, *temp_);
    op2_->apply( x, y);
    y.add( 1., *temp_, 1., y );
  }

  void assignField( const DomainFieldT& mv ) {
    op1_->assignField( mv );
    op2_->assignField( mv );
  };

	Teuchos::RCP<const SpaceT> space() const { return(op1_->space()); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		op1_->setParameter( para );
		op2_->setParameter( para );
	}

  bool hasApplyTranspose() const { return( op1_->hasApplyTranspose() && op2_->hasApplyTranspose() ); }

}; // end of class Add2Op



/// \relates Add2Op
template<class OP1, class OP2 >
Teuchos::RCP< Add2Op<OP1, OP2> > createAdd2Op(
    const Teuchos::RCP<OP1>& op1,
    const Teuchos::RCP<OP2>& op2 ) {
  return( Teuchos::rcp( new Add2Op<OP1,OP2>(op1,op2) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_ADD2OP_HPP
