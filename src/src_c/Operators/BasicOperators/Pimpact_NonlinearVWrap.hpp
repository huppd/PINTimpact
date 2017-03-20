#pragma once
#ifndef PIMPACT_CONVECTIONVWRAP_HPP
#define PIMPACT_CONVECTIONVWRAP_HPP


#include "Pimpact_ConvectionSOp.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {



/// \brief Convection Wraper of  for Velocity fields
/// \ingroup BaseOperator
/// \relates ConvectionSOp
template<class SOpT>
class NonlinearWrap {

public:

  using SpaceT = typename SOpT::SpaceT;

  using DomainFieldT = VectorField<SpaceT>;
  using RangeFieldT = VectorField<SpaceT>;

  using FieldTensor = ScalarField<SpaceT>[3][3];

protected:

  using Scalar = typename SpaceT::Scalar;

  Teuchos::RCP<SOpT> convectionSOp_;

public:

	NonlinearWrap( const Teuchos::RCP<SOpT>& sop ): convectionSOp_( sop ) {};


  /// \note Operator's wind has to be assigned correctly
  void apply( const FieldTensor& u, const DomainFieldT& x, RangeFieldT& y, const Add& add=Add::N ) const {

    for( F i=F::U; i<SpaceT::sdim; ++i )
      convectionSOp_->apply( u[static_cast<int>(i)], x(i), y(i), add );
  }

  /// \note Operator's wind has to be assigned correctly
  void apply( const FieldTensor& u, const DomainFieldT& x, RangeFieldT& y,
			Scalar mulI, Scalar mulC, Scalar mulL, const Add& add=Add::N ) const {

    for( F i=F::U; i<SpaceT::sdim; ++i )
      convectionSOp_->apply( u[static_cast<int>(i)], x(i), y(i), mulI, mulC, mulL, add );
  }

  constexpr const Teuchos::RCP<const SpaceT>& space() const { return( convectionSOp_->space() ); }

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		convectionSOp_->setParameter( para );
	}

  constexpr const Teuchos::RCP<const SOpT> getSOp() const {
    return( convectionSOp_ );
  }

  void print( std::ostream& out=std::cout ) const {
    out << "--- NonlinearWrap(" << getLabel() << ") ---\n";
    convectionSOp_->print(out);
  }

	const std::string getLabel() const { return( convectionSOp_->getLabel() ); };

}; // end of class NonlinearWrap



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONVWRAP_HPP
