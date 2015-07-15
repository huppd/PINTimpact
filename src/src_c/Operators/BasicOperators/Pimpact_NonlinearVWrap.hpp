#pragma once
#ifndef PIMPACT_CONVECTIONVWRAP_HPP
#define PIMPACT_CONVECTIONVWRAP_HPP


#include "Pimpact_VectorField.hpp"

#include "Pimpact_ConvectionSOp.hpp"




namespace Pimpact {


/// \brief Convection Wraper of  for Velocity fields
/// \ingroup BaseOperator
/// \relates ConvectionSOp
template<class SOpT>
class NonlinearWrap {

public:

  typedef typename SOpT::SpaceT SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  static const int dimension = SpaceT::dimension;

  static const int dimNC = SpaceT::dimNC;

  typedef VectorField<SpaceT>  DomainFieldT;
  typedef VectorField<SpaceT>  RangeFieldT;

  typedef Teuchos::Tuple< Teuchos::Tuple<Teuchos::RCP<ScalarField<SpaceT> >, 3>, 3> FieldTensor;

protected:

  Teuchos::RCP<SOpT> convectionSOp_;

public:

  NonlinearWrap( const Teuchos::RCP<SOpT>& sop ):
    convectionSOp_( sop ) {};


  /// \note Operator's wind has to be assigned correctly
  void apply( const FieldTensor& u, const DomainFieldT& x, RangeFieldT& y, Scalar mul=0., Scalar mulC=1. ) const {

    for( int i=0; i<x.space()->dim(); ++i ) {
      convectionSOp_->apply( u[i], x.getConstField(i), y.getField(i), mul );
    }

  }

  /// \note Operator's wind has to be assigned correctly
  void apply( const FieldTensor& u, const DomainFieldT& x, RangeFieldT& y,
			Scalar mul, Scalar mulI, Scalar mulC, Scalar mulL ) const {

    for( int i=0; i<x.space()->dim(); ++i ) {
      convectionSOp_->apply( u[i], x.getConstField(i), y.getField(i), mul, mulI, mulC, mulL );
    }

  }

  Teuchos::RCP<const SpaceT> space() const { return( convectionSOp_->space() ); }

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		convectionSOp_->setParameter( para );
	}

  Teuchos::RCP<const SOpT> getSOp() const {
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
