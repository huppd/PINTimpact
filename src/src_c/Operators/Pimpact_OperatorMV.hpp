#pragma once
#ifndef PIMPACT_OPERATORMV_HPP
#define PIMPACT_OPERATORMV_HPP


#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_Types.hpp"

#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiField.hpp"



namespace Pimpact {


/// \brief Operator for \c MultiField's.
///
/// wraps and \c Operator and adds the functionality of handling \c MultiField.
/// Aditionally it can be use as Operator for Belos, wich is \deprecated,
//because it does not allow preconditionin from another operator type.
template<class Operator>
class OperatorMV  {

	typedef typename Operator::DomainFieldT DomainFieldT;
	typedef typename Operator::RangeFieldT RangeFieldT;
	typedef typename Operator::OpType opType;

	Teuchos::RCP<Operator> op_;

public:


	~OperatorMV() {op_=Teuchos::null;};


	OperatorMV( const Teuchos::RCP<Operator>& op ):op_(op) {};


  /// \brief default apply
	void apply( const MultiField<DomainFieldT>& x,
			MultiField<RangeFieldT> & y,
			Belos::ETrans trans=Belos::NOTRANS) const {

		// Debug test: y.GetNumVecs()>=x.GetNumVecs()
		int m = x.getNumberVecs();

		for( int i=0; i<m; ++i )
			op_->apply( x.getConstField(i), y.getField(i) );
	}


	/// \brief apply for \c ModeField's.
	///
	/// depending on the \c opType it will apply an \c ModeOp or an \c NonModeOP,
	//which means that the \c op is applied on each mode separetely.
	void apply(
			const MultiField<ModeField<DomainFieldT> >& x,
			MultiField<ModeField<RangeFieldT> > & y,
			Belos::ETrans trans=Belos::NOTRANS ) const {
		apply( opType(), x, y, trans );
	}


	bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

	Teuchos::RCP<Operator> getOperatorPtr() { return( op_ ); }

private:

	void apply( ModeOp,
			const MultiField<ModeField<DomainFieldT> >& x,
			MultiField<ModeField<RangeFieldT> > & y,
			Belos::ETrans trans=Belos::NOTRANS) const {

		// Debug test: y.GetNumVecs()>=x.GetNumVecs()
		int m = x.getNumberVecs();

		for( int i=0; i<m; ++i )
			op_->apply( x.getConstField(i), y.getField(i) );
	}


	void apply( NonModeOp,
			const MultiField<ModeField<DomainFieldT> >& x,
			MultiField<ModeField<RangeFieldT> > & y,
			Belos::ETrans trans=Belos::NOTRANS ) const {

		// Debug test: y.GetNumVecs()>=x.GetNumVecs()
		int m = x.getNumberVecs();

		for( int i=0; i<m; ++i ) {
			op_->apply( *x.getConstField(i).getConstFieldC(), *y.getField(i).getFieldC() );
			op_->apply( *x.getConstField(i).getConstFieldS(), *y.getField(i).getFieldS() );
		}
	}

}; // end of class OperatorMV


template<class Operator>
Teuchos::RCP< OperatorMV<Operator> > createOperatorMV() {
	return( Teuchos::rcp( new OperatorMV<Operator>( Teuchos::rcp( new Operator() ) ) ) );
}


template<class Operator>
Teuchos::RCP< OperatorMV<Operator> > createOperatorMV( const Teuchos::RCP<Operator>& op) {
	return( Teuchos::rcp( new OperatorMV<Operator>( op ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_OPERATORMV_HPP
