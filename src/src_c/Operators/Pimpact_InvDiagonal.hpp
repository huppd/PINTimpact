#pragma once
#ifndef PIMPACT_INVDIAGONAL_HPP
#define PIMPACT_INVDIAGONAL_HPP


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "BelosTypes.hpp"




namespace Pimpact {


/// \brief transposes operator.
///
/// the \c DomainFieldT  \c RangeFieldT is to be equal for both \c OP.
/// \ingroup Operator
template< class OP >
class InvDiagonal {

public:

  using DomainFieldT = typename OP::DomainFieldT;
  using RangeFieldT = typename OP::RangeFieldT;

  using SpaceT = typename DomainFieldT::SpaceT;

protected:

  Teuchos::RCP<OP> op_;

public:

	InvDiagonal( const Teuchos::RCP<OP>& op ): op_(op) {
	};


	void apply(const DomainFieldT& x, RangeFieldT& y,
			Belos::ETrans trans=Belos::NOTRANS ) const {

		op_->applyInvDiag( x, y );
  }

  void assignField( const DomainFieldT& mv ) {
    op_->assignField( mv );
  };

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(op_->space()); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		op_->setParameter( para );
	}

  bool hasApplyTranspose() const { return( true ); }

	const std::string getLabel() const { return( op_->getLabel() + std::string("^T")  ); };

  void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
		op_->print( out );
  }

}; // end of class InvDiagonal


/// \relates InvDiagonal
template<class OP>
Teuchos::RCP< InvDiagonal<OP> > createInvDiagonal(
    const Teuchos::RCP<OP>& op ) {
  return( Teuchos::rcp( new InvDiagonal<OP>(op) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INVDIAGONAL_HPP
