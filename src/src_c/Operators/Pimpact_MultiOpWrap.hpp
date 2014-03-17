#pragma once
#ifndef PIMPACT_MULTIOPWRAP_HPP
#define PIMPACT_MULTIOPWRAP_HPP


#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_Types.hpp"

#include "Pimpact_MultiField.hpp"



namespace Pimpact {


/// \brief Operator wrapper for \c MultiField's.
/// \ingroup MultiOperator
/// wraps and \c Operator and adds the functionality of handling \c MultiField.
template<class Operator>
class MultiOpWrap  {

	Teuchos::RCP<Operator> op_;

public:

	typedef MultiField<typename Operator::DomainFieldT> DomainFieldT;
	typedef MultiField<typename Operator::RangeFieldT> RangeFieldT;
	typedef typename Operator::OpType OpType;

	MultiOpWrap():op_( Teuchos::rcp( new Operator() ) ) {};
	MultiOpWrap( const Teuchos::RCP<Operator>& op ):op_(op) {};
	~MultiOpWrap() {op_=Teuchos::null;};


  /// \brief default apply
	void apply( const DomainFieldT& x,
			RangeFieldT& y,
			Belos::ETrans trans=Belos::NOTRANS) const {

		// Debug test: y.GetNumVecs()>=x.GetNumVecs()
		int m = x.getNumberVecs();

		for( int i=0; i<m; ++i )
			op_->apply( x.getConstField(i), y.getField(i) );
	}

  void assignField( const DomainFieldT& mv ) {
    op_->assignField( mv.getConstField(0) );
  };

	bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

	Teuchos::RCP<Operator> getOperatorPtr() { return( op_ ); }


}; // end of class MultiOpWrap



template<class Operator>
Teuchos::RCP< MultiOpWrap<Operator> > createMultiOpWrap( const Teuchos::RCP<Operator>& op=Teuchos::null) {
  if( Teuchos::is_null(op) )
    return( Teuchos::rcp( new MultiOpWrap<Operator>( Teuchos::rcp( new Operator() ) ) ) );
  else
    return( Teuchos::rcp( new MultiOpWrap<Operator>( op ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIOPWRAP_HPP
