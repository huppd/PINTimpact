#pragma once
#ifndef PIMPACT_OPERATORMV_HPP
#define PIMPACT_OPERATORMV_HPP


/** \file Pimpact_OperatorMV.hpp
*/

//#include <Teuchos_Assert.hpp>
#include "Teuchos_RCP.hpp"
//#include "Teuchos_LocalTestingHelpers.hpp"

//#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>

#include "Pimpact_Types.hpp"
//#include "Pimpact_ScalarField.hpp"
//#include "Pimpact_VectorField.hpp"
#include "Pimpact_ModeField.hpp"
#include "Pimpact_MultiField.hpp"


namespace Pimpact {



/**
 * has multiple \c Field types, where \c Field can be a \c Pimpact:ScalarField or \c Pimpact:VectorField
 */
template<class Operator>
class OperatorMV {

	typedef typename Operator::DomainFieldT DomainFieldT;
	typedef typename Operator::RangeFieldT RangeFieldT;
	typedef typename Operator::OpType opType;

	Teuchos::RCP<Operator> op_;

public:

	OperatorMV( Teuchos::RCP<Operator> op ):op_(op) {};

	void apply( const MultiField<DomainFieldT>& x,
			MultiField<RangeFieldT> & y,
			Belos::ETrans trans=Belos::NOTRANS) const {

		// Debug test: y.GetNumVecs()>=x.GetNumVecs()
		int m = x.GetNumberVecs();

		for( int i=0; i<m; ++i )
			op_->apply( x.GetConstVec(i), y.GetVec(i) );
	}


	void apply(
			const MultiField<ModeField<DomainFieldT> >& x,
			MultiField<ModeField<RangeFieldT> > & y,
			Belos::ETrans trans=Belos::NOTRANS ) const {
		apply( opType(), x, y, trans );
	}


	bool hasApplyTranspose() const { return op_.hasApplyTranspose(); }

private:

	void apply( ModeOp,
			const MultiField<ModeField<DomainFieldT> >& x,
			MultiField<ModeField<RangeFieldT> > & y,
			Belos::ETrans trans=Belos::NOTRANS) const {

		// Debug test: y.GetNumVecs()>=x.GetNumVecs()
		int m = x.GetNumberVecs();

		for( int i=0; i<m; ++i )
			op_->apply( x.GetConstVec(i), y.GetVec(i) );
	}


	void apply( NonModeOp,
			const MultiField<ModeField<DomainFieldT> >& x,
			MultiField<ModeField<RangeFieldT> > & y,
			Belos::ETrans trans=Belos::NOTRANS ) const {

		// Debug test: y.GetNumVecs()>=x.GetNumVecs()
		int m = x.GetNumberVecs();

		for( int i=0; i<m; ++i ) {
			op_->apply( *x.GetConstVec(i).getConstFieldC(), *y.GetVec(i).getFieldC() );
			op_->apply( *x.GetConstVec(i).getConstFieldS(), *y.GetVec(i).getFieldS() );
		}
	}

}; // end of class OperatorMV


template<class Operator>
Teuchos::RCP< OperatorMV<Operator> > createOperatorMV() {
	return Teuchos::rcp( new OperatorMV<Operator>( Teuchos::rcp( new Operator() ) ) );
}

} // end of Pimpact namespace

#endif // end of file PIMPACT_OPERATORMV_HPP
