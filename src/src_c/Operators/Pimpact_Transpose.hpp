#pragma once
#ifndef PIMPACT_TRANSPOSE_HPP
#define PIMPACT_TRANSPOSE_HPP


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "BelosTypes.hpp"




namespace Pimpact {


/// \brief transposes operator.
///
/// the \c DomainFieldT  \c RangeFieldT is to be equal for both \c OP.
/// \ingroup Operator
template< class OP >
class Transpose {

public:

  using DomainFieldT = typename OP::DomainFieldT;
  using RangeFieldT = typename OP::RangeFieldT;

  using SpaceT = typename DomainFieldT::SpaceT;

protected:

  Teuchos::RCP<OP> op_;

public:

  Transpose( const Teuchos::RCP<OP>& op ): op_(op) {
    assert( op->hasApplyTranspose() );
  };


  void apply(const DomainFieldT& x, RangeFieldT& y,
             const Belos::ETrans& trans=Belos::NOTRANS ) const {

    switch( trans ) {
    case Belos::NOTRANS : {
      op_->apply( x, y, Belos::TRANS );
      break;
    }
    case Belos::TRANS : {
      op_->apply( x, y, Belos::NOTRANS  );
      break;
    }
    case Belos::CONJTRANS : {
      op_->apply( x, y, Belos::NOTRANS  );
      break;
    }
    }
  }

  void assignField( const DomainFieldT& mv ) {
    op_->assignField( mv );
  };

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return(op_->space());
  };

  void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
    op_->setParameter( para );
  }

  bool hasApplyTranspose() const {
    return( op_->hasApplyTranspose() );
  }

  const std::string getLabel() const {
    return( op_->getLabel() + std::string("^T")  );
  };

  void print( std::ostream& out=std::cout ) const {
    out << getLabel() << ":\n";
    op_->print( out );
  }

}; // end of class Transpose


/// \relates Transpose
template<class OP>
Teuchos::RCP< Transpose<OP> > createTranspose(
  const Teuchos::RCP<OP>& op ) {
  return( Teuchos::rcp( new Transpose<OP>(op) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_TRANSPOSE_HPP
