#pragma once
#ifndef PIMPACT_CONVECTIONVSMOOTHER_HPP
#define PIMPACT_CONVECTIONVSMOOTHER_HPP


#include "Pimpact_ConvectionVWrap.hpp"
#include "Pimpact_ConvectionField.hpp"



namespace Pimpact {



/// \brief Convection Operator for Velocity fields
/// \todo make wind template parameter as well.(necessary when different winds
/// are wanted, meaning moving interpolation steps from assign to apply
/// \todo make constructor so wind can be shared by different operators.
/// \ingroup BaseOperator
/// \ingroup NonlinearOperator
template<class ConvVOpT, template<class> class ST>
class ConvectionVSmoother {

public:

  typedef typename ConvVOpT::SpaceT SpaceT;


  typedef ST<typename ConvVOpT::ConvSOpT> SSmootherT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  static const int dimension = SpaceT::dimension;

  static const int dimNC = SpaceT::dimNC;

  typedef VectorField<SpaceT>  DomainFieldT;
  typedef VectorField<SpaceT>  RangeFieldT;

protected:

  Teuchos::RCP<ConvectionVWrap<SSmootherT> > convVWrap_;

  Teuchos::RCP< ConvectionField<SpaceT> > convField_;

public:

    ConvectionVSmoother(
        const Teuchos::RCP<const ConvVOpT>& op,
        Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ):
      convVWrap_( create<ConvectionVWrap<SSmootherT> >( create<SSmootherT>( op->getSOp(), pl ) ) ),
      convField_( op->getConvField() ) {};



  /// should not be already assigned in Operators.
  void assignField( const DomainFieldT& mv ) const { convField_->assignField( mv ); };

  /// \note Operator's wind has to be assigned correctly
  /// \deprecated
  void apply( const DomainFieldT& x, RangeFieldT& y ) const {

    convVWrap_->apply( convField_->get(), x, y );

  }

  //void apply(const DomainFieldT& z, const DomainFieldT& x, RangeFieldT& y, Scalar mul=0. ) const {};

  Teuchos::RCP< ConvectionField<SpaceT> >
  getConvField() const {
    return( convField_ );
  }


  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(convVWrap_->space()); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		convVWrap_->setParameter( para );
	}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    convVWrap_->print(out);
  }

	const std::string getLabel() const { return( "ConvectionVSmoother " ); };

}; // end of class ConvectionVSmoother


} // end of namespace Pimpact



#ifdef COMPILE_ETI
#include "Pimpact_ConvectionDiffusionSORSmoother.hpp"
extern template class Pimpact::ConvectionVSmoother< Pimpact::ConvectionVOp< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,3,2> > >, Pimpact::ConvectionDiffusionSORSmoother>;
extern template class Pimpact::ConvectionVSmoother< Pimpact::ConvectionVOp< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,3,4> > >, Pimpact::ConvectionDiffusionSORSmoother>;
extern template class Pimpact::ConvectionVSmoother< Pimpact::ConvectionVOp< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,4,2> > >, Pimpact::ConvectionDiffusionSORSmoother>;
extern template class Pimpact::ConvectionVSmoother< Pimpact::ConvectionVOp< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,4,4> > >, Pimpact::ConvectionDiffusionSORSmoother>;
#endif



#endif // end of #ifndef PIMPACT_CONVECTIONVSMOOTHER_HPP
