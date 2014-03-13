#pragma once
#ifndef PIMPACT_MODEOPWRAP_HPP
#define PIMPACT_MODEOPWRAP_HPP


#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_Types.hpp"

#include "Pimpact_ModeField.hpp"



namespace Pimpact {


/// \brief Operator wrapper for \c ModeField's.
/// \ingroup ModeOperator
/// wraps and \c Operator and adds the functionality of handling \c ModeField.
template<class Operator>
class ModeOpWrap  {

	Teuchos::RCP<Operator> op_;

public:

	typedef ModeField<typename Operator::DomainFieldT> DomainFieldT;
	typedef ModeField<typename Operator::RangeFieldT>  RangeFieldT;
	typedef typename Operator::OpType OpType;

	ModeOpWrap():op_( Teuchos::rcp( new Operator() ) ) {};
	ModeOpWrap( const Teuchos::RCP<Operator>& op ):op_(op) {};
	~ModeOpWrap() {op_=Teuchos::null;};


	void apply( const DomainFieldT& x, RangeFieldT& y,
			Belos::ETrans trans=Belos::NOTRANS) const {

	  op_->apply( x.getConstCField(), y.getCField() );
    op_->apply( x.getConstSField(), y.getSField() );

	}


	bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

	Teuchos::RCP<Operator> getOperatorPtr() { return( op_ ); }


}; // end of class ModeOpWrap



template<class Operator>
Teuchos::RCP< ModeOpWrap<Operator> > createModeOpWrap( const Teuchos::RCP<Operator>& op=Teuchos::null) {
  if( Teuchos::is_null(op) )
    return( Teuchos::rcp( new ModeOpWrap<Operator>( Teuchos::rcp( new Operator() ) ) ) );
  else
    return( Teuchos::rcp( new ModeOpWrap<Operator>( op ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MODEOPWRAP_HPP
