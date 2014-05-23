#pragma once
#ifndef PIMPACT_MULTIHARMONICMLPREC_HPP
#define PIMPACT_MULTIHARMONICMLPREC_HPP


#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_Types.hpp"

#include "Pimpact_MultiField.hpp"



namespace Pimpact {


/// \brief Operator wrapper.
/// \ingroup MultiHarmonicOperator
/// wraps a \ref BaseOperator "base operator" and adds the functionality of handling \c MultiHarmonicField's.
template<class Field>
class MultiHarmonicMLPrec  {

public:

  typedef OperatorBase<MultiField<Field> > Op0;
  typedef OperatorBase<MultiField< ModeField<Field> > > Ops;

private:

  Teuchos::RCP<Op0> op0_;
  Teuchos::Array<Teuchos::RCP<Ops> > ops_;

public:

  typedef MultiHarmonicField<Field> DomainFieldT;
  typedef MultiHarmonicField<Field> RangeFieldT;

  MultiHarmonicMLPrec():
        op0_(Teuchos::null), ops_(0) {
  };

  MultiHarmonicMLPrec(
      const Teuchos::RCP<Op0>& op0,
      const Teuchos::Array<Teuchos::RCP<Ops> >& ops ):
        op0_(op0), ops_(ops) {
  };


  void apply( const DomainFieldT& x,
      RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS) const {

    op0_->apply( *createMultiField(Teuchos::rcp_const_cast<Field>(x.getConst0FieldPtr())), *createMultiField(y.get0FieldPtr()) );

    int m = x.getNumberModes();

    for( int i=0; i<m; ++i ) {
      ops_[i]->apply( *createMultiField(Teuchos::rcp_const_cast<ModeField<Field> >(x.getConstFieldPtr(i))), *createMultiField(y.getFieldPtr(i)) );
    }
  };

  void assignField( const DomainFieldT& mv ) {
//    op0_->assignField( mv.getConst0Field() );
//    for( int i=0; i<mv.getNumberVecs(); ++i ) {
//      ops_[i]->assignField( mv.getConstField() );
//    }
  };

  bool hasApplyTranspose() const { return( false ); }

//  Teuchos::RCP<Operator> getOperatorPtr() { return( op_ ); }


}; // end of class MultiHarmonicMLPrec



/// \relates MultiHarmonicMLPrec
template<class Field>
Teuchos::RCP< MultiHarmonicMLPrec<Field> > createMultiHarmonicMLPrec(
    const Teuchos::RCP<typename MultiHarmonicMLPrec<Field>::Op0>& op0=Teuchos::null,
    const Teuchos::Array<Teuchos::RCP<typename MultiHarmonicMLPrec<Field>::Ops> >& ops=Teuchos::null ) {

    return( Teuchos::rcp( new MultiHarmonicMLPrec<Field>( op0, ops ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIHARMONICMLPREC_HPP
