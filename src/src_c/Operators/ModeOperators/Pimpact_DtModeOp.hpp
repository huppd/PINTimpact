#pragma once
#ifndef PIMPACT_DTMODEOP_HPP
#define PIMPACT_DTMODEOP_HPP

#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_ModeField.hpp"



namespace Pimpact{


/// \ingroup ModeOperator
template<class ST>
class DtModeOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

protected:

	const Teuchos::RCP<const SpaceT> space_;

  Scalar alpha2_;

public:

  DtModeOp( const Teuchos::RCP<const SpaceT>& space ):
		space_(space),
		alpha2_( space->getDomain()->getDomainSize()->getAlpha2()/space->getDomain()->getDomainSize()->getRe() ) {};

  typedef ModeField<VectorField<SpaceT> >  DomainFieldT;
  typedef ModeField<VectorField<SpaceT> >  RangeFieldT;


  void apply(const DomainFieldT& x, RangeFieldT& y ) const {
    y.getCFieldPtr()->add(      0.,  x.getConstCField(), alpha2_, x.getConstSField() );
    y.getSFieldPtr()->add( -alpha2_,  x.getConstCField(),     0., x.getConstSField() );
  }

  void assignField( const DomainFieldT& mv ) {};

	Teuchos::RCP<const SpaceT> space() const { return(space_); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "DtModeOp " ); };

}; // end of class DtModeOp



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DTMODEOP_HPP
