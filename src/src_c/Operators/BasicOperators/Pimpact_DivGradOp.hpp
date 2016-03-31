#pragma once
#ifndef PIMPACT_DIVGRADOP_HPP
#define PIMPACT_DIVGRADOP_HPP


#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{



/// \brief "laplace" for pressure.
/// \ingroup BaseOperator
/// \todo not workin properly?
/// \warning does not hold test.
template<class ST>
class DivGradOp {

public:

  using SpaceT = ST;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

  Teuchos::RCP<DivOp<SpaceT> > div_;
  Teuchos::RCP<GradOp<SpaceT> > grad_;

public:

  DivGradOp( const Teuchos::RCP<const SpaceT>& space ):
    div_ ( create<DivOp>( space ) ),
    grad_( create<GradOp>( space ) ) {};

  DivGradOp(
      const Teuchos::RCP< DivOp<SpaceT> >& div,
      const Teuchos::RCP< GradOp<SpaceT> >& grad ):
    div_ (div),
    grad_(grad) {};

  void apply(const DomainFieldT& x, RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS ) const {

		Teuchos::RCP< VectorField<SpaceT> > temp = create<VectorField>( space() );

    grad_->apply( x, *temp );
    div_->apply( *temp, y );

  }

	void computeResidual( const RangeFieldT& b, const DomainFieldT& x, RangeFieldT& res ) const {
		apply( x, res );
		res.add( 1., b, -1., res );
	}

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return( div_->space() ); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

	constexpr const std::string getLabel() const { return( "DivGrad" ); };

  void print( std::ostream& out=std::cout ) const {
		out << "---" << getLabel() << "---\n";
		div_->print( out );
		grad_->print( out );

  }

}; // end of DivGradOp



/// \relates DivGradOp
template<class SpaceT>
Teuchos::RCP< DivGradOp<SpaceT> > createDivGradOp(
    const Teuchos::RCP< DivOp<SpaceT> >& div,
    const Teuchos::RCP< GradOp<SpaceT> >& grad ) {

  return(
      Teuchos::rcp( new DivGradOp<SpaceT>( div, grad ) )
  );

}


} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::DivGradOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::DivGradOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::DivGradOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::DivGradOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_DIVGRADOP_HPP
