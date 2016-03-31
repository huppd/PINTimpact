#pragma once
#ifndef PIMPACT_DTMODEOP_HPP
#define PIMPACT_DTMODEOP_HPP


#include "Pimpact_ModeField.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{


/// \ingroup ModeOperator
template<class ST>
class DtModeOp {

public:

  using SpaceT = ST;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

protected:

	const Teuchos::RCP<const SpaceT> space_;

  Scalar alpha2_;

public:

  DtModeOp( const Teuchos::RCP<const SpaceT>& space ):
		space_(space),
		alpha2_( space->getDomainSize()->getAlpha2()/space->getDomainSize()->getRe() ) {};

  using DomainFieldT = ModeField<VectorField<SpaceT> >;
  using RangeFieldT = ModeField<VectorField<SpaceT> >;


  void apply(const DomainFieldT& x, RangeFieldT& y ) const {
    y.getCFieldPtr()->add(      0.,  x.getConstCField(), alpha2_, x.getConstSField() );
    y.getSFieldPtr()->add( -alpha2_,  x.getConstCField(),     0., x.getConstSField() );
  }

  void assignField( const DomainFieldT& mv ) {};

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(space_); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "DtModeOp " ); };

  void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
		out << "alpha2: " << alpha2_ << "\n";
  }

}; // end of class DtModeOp



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DTMODEOP_HPP
