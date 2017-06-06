#pragma once
#ifndef PIMPACT_NONLINEARSMOOTHER_HPP
#define PIMPACT_NONLINEARSMOOTHER_HPP


#include "Pimpact_ConvectionField.hpp"
#include "Pimpact_NonlinearVWrap.hpp"




namespace Pimpact {



/// \brief Convection Operator for Velocity fields
/// \note make wind template parameter as well.(necessary when different winds
/// are wanted, meaning moving interpolation steps from assign to apply
/// \note todo make constructor so wind can be shared by different operators.
/// \ingroup BaseOperator
/// \ingroup NonlinearOperator
template<class ConvVOpT, template<class> class SmootherT>
class NonlinearSmoother {

public:

  using SpaceT = typename ConvVOpT::SpaceT;

  using DomainFieldT = VectorField<SpaceT>;
  using RangeFieldT = VectorField<SpaceT>;

  using SSmootherT = SmootherT<typename ConvVOpT::ConvSOpT>;

protected:

  Teuchos::RCP< NonlinearWrap<SSmootherT> > convVWrap_;

  Teuchos::RCP< ConvectionField<SpaceT> > convField_;

public:

  NonlinearSmoother(
    const Teuchos::RCP<const ConvVOpT>& op,
    Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ):
    convVWrap_( create<NonlinearWrap<SSmootherT> >( create<SSmootherT>( op->getSOp(), pl ) ) ),
    convField_( op->getConvField() ) {};


  /// NOFX should not be already assigned in Operators.
  void assignField( const DomainFieldT& mv ) const {
  };


  /// \note Operator's wind has to be assigned correctly
  void apply( const DomainFieldT& x, RangeFieldT& y ) const {

    convVWrap_->apply( convField_->get(), x, y );
  }


  Teuchos::RCP< ConvectionField<SpaceT> > getConvField() const {
    return( convField_ );
  }


  bool hasApplyTranspose() const {
    return( false );
  }

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return(convVWrap_->space());
  };

  void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
    convVWrap_->setParameter( para );
  }

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    convVWrap_->print(out);
  }

  const std::string getLabel() const {
    return( convVWrap_->getLabel() );
  };

}; // end of class NonlinearSmoother


} // end of namespace Pimpact



#ifdef COMPILE_ETI
#include "Pimpact_ConvectionDiffusionSORSmoother.hpp"
extern template class Pimpact::NonlinearSmoother< Pimpact::NonlinearOp< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,3,2> > >, Pimpact::ConvectionDiffusionSORSmoother>;
extern template class Pimpact::NonlinearSmoother< Pimpact::NonlinearOp< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,3,4> > >, Pimpact::ConvectionDiffusionSORSmoother>;
extern template class Pimpact::NonlinearSmoother< Pimpact::NonlinearOp< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,4,2> > >, Pimpact::ConvectionDiffusionSORSmoother>;
extern template class Pimpact::NonlinearSmoother< Pimpact::NonlinearOp< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,4,4> > >, Pimpact::ConvectionDiffusionSORSmoother>;
#endif


#endif // end of #ifndef PIMPACT_NONLINEARSMOOTHER_HPP
