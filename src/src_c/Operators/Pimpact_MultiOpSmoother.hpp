#pragma once
#ifndef PIMPACT_MultiOpSmoother_HPP
#define PIMPACT_MultiOpSmoother_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_MultiField.hpp"




namespace Pimpact {


/// \brief Operator wrapper for \c MultiField's.
/// \ingroup MultiOperator
/// wraps and \c Operator and adds the functionality of handling \c MultiField.
template<class OT>
class MultiOpSmoother  {

public:

  typedef OT OperatorT;

  typedef typename OperatorT::SpaceT SpaceT;

  typedef MultiField<typename OperatorT::DomainFieldT> DomainFieldT;
  typedef MultiField<typename OperatorT::RangeFieldT> RangeFieldT;

protected:

  Teuchos::RCP<OperatorT> op_;

public:

//	MultiOpSmoother( const Teuchos::RCP<OperatorT>& op ):op_(op) {}

	template<class IOperatorT>
	MultiOpSmoother( const Teuchos::RCP< MultiOpSmoother<IOperatorT> >& op ):
		op_( Teuchos::rcp( new OperatorT( op->getOperatorPtr() ) ) ) {}


	void apply( const DomainFieldT& x, RangeFieldT& y,
			Belos::ETrans trans=Belos::NOTRANS) const {

    TEUCHOS_TEST_FOR_EXCEPT( x.getNumberVecs()!=y.getNumberVecs() );

    for( int i=0; i<x.getNumberVecs(); ++i )
      op_->apply( x.getConstField(i), y.getField(i) );

  }


  void assignField( const DomainFieldT& mv ) {

    op_->assignField( mv.getConstField(0) );

  };


  bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

  Teuchos::RCP<OperatorT> getOperatorPtr() { return( op_ ); }

  Teuchos::RCP<const SpaceT> space() const { return( op_->space() ); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		op_->setParameter( para );
	}

	const std::string getLabel() const { return( op_->getLabel()  ); };

  void print( std::ostream& out=std::cout ) const {
    op_->print( out );
  }

}; // end of class MultiOpSmoother



///// \relates MultiOpSmoother
///// \deprecated use create<MultiOpSmoother>( op ) instead
//template<class OperatorT>
//Teuchos::RCP< MultiOpSmoother<OperatorT> > createMultiOpSmoother( const Teuchos::RCP<OperatorT>& op ) {
//
//	 return( Teuchos::rcp( new MultiOpSmoother<OperatorT>( op ) ) );
//
//}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MultiOpSmoother_HPP
