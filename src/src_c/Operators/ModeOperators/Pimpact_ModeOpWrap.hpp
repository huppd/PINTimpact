#pragma once
#ifndef PIMPACT_MODEOPWRAP_HPP
#define PIMPACT_MODEOPWRAP_HPP


#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_ModeField.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {


/// \brief Operator wrapper for \c ModeField's.
/// \ingroup ModeOperator
/// wraps and \c Operator and adds the functionality of handling \c ModeField.
template<class Operator>
class ModeOpWrap {

public:

  using DomainFieldT = ModeField<typename Operator::DomainFieldT>;
  using RangeFieldT = ModeField<typename Operator::RangeFieldT>;

  using SpaceT = typename DomainFieldT::SpaceT;

protected:

  Teuchos::RCP<Operator> op_;

public:

  ModeOpWrap():op_( Teuchos::rcp( new Operator() ) ) {};
  ModeOpWrap( const Teuchos::RCP<Operator>& op ):op_(op) {};
  ~ModeOpWrap() {
    op_=Teuchos::null;
  };


  /// \note todo apply for helmholtx k=...
  void apply( const DomainFieldT& x, RangeFieldT& y, const Belos::ETrans
      trans=Belos::NOTRANS) const {

    op_->apply( x.getCField(), y.getCField() );
    op_->apply( x.getSField(), y.getSField() );
  }

  void assignField( const DomainFieldT& mv ) {
    op_->assignField( mv.getConstCField() );
  };


  bool hasApplyTranspose() const {
    return op_->hasApplyTranspose();
  }

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return op_->space();
  };

  void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  Teuchos::RCP<Operator> getOperatorPtr() {
    return op_;
  }

  const std::string getLabel() const {
    return "ModeOpWrap( "+op_->getLabel()+" ) ";
  };

  void print( std::ostream& out=std::cout ) const {
    out << getLabel() << ":\n";
    op_->print( out );
  }

}; // end of class ModeOpWrap




/// \relates ModeOpWrap
template<class Operator>
Teuchos::RCP< ModeOpWrap<Operator> > createModeOpWrap( const Teuchos::RCP<Operator>& op ) {
  return Teuchos::rcp( new ModeOpWrap<Operator>( op ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MODEOPWRAP_HPP
