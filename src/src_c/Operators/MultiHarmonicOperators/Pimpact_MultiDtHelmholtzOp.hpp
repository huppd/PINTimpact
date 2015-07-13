#pragma once
#ifndef PIMPACT_MULTIDTHELMHOLTZOP_HPP
#define PIMPACT_MULTIDTHELMHOLTZOP_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"

#include "Pimpact_DtLapOp.hpp"


namespace Pimpact{


/// \ingroup MultiHarmonicOperator
/// for generalizing this(MultiHarmonicwrapper for ModeOperator), ModeOperator needs a non Mode representation
/// \relates DtLapOp
template<class ST>
class MultiDtHelmholtz {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;

  typedef MultiHarmonicField< VectorField<SpaceT> >  DomainFieldT;
  typedef MultiHarmonicField< VectorField<SpaceT> >  RangeFieldT;

protected:

  Teuchos::RCP<DtLapOp<SpaceT> > op_;

public:

  MultiDtHelmholtz( const Teuchos::RCP<const SpaceT>& space, Scalar alpha2=1., Scalar iRe=1.):
    op_( createDtLapOp<SpaceT>( space, alpha2, iRe) ) {};

  MultiDtHelmholtz( const Teuchos::RCP<DtLapOp<SpaceT> >& op ):
    op_(op) {};


  void apply(const DomainFieldT& x, RangeFieldT& y ) const {
    op_->getInnerOpPtr()->apply( x.getConst0Field(), y.get0Field() );

    for( int i=0; i<x.getNumberModes(); ++i ) {
      op_->apply( x.getConstField(i), y.getField(i), i+1 );
      op_->apply( x.getConstField(i), y.getField(i), i+1 );
    }
  }


  void assignField( const DomainFieldT& mv ) {};


  Teuchos::RCP< HelmholtzOp<SpaceT> > getInnerOpPtr() {
    return( op_ );
  }

	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "MultiDtHelmholtz" ); };

  void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
		op_->print( out );
  }

}; // end of class MultiDtHelmholtz



/// \relates MultiDtHelmholtz
template<class SpaceT>
Teuchos::RCP< MultiDtHelmholtz<SpaceT> > createMultiDtHelmholtz(
    const Teuchos::RCP< const SpaceT>& space,
    typename SpaceT::Scalar omega=1.,
    typename SpaceT::Scalar mulL=1. ) {
  return( Teuchos::rcp( new MultiDtHelmholtz<SpaceT>( space, omega, mulL ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIDTHELMHOLTZOP_HPP
