#pragma once
#ifndef PIMPACT_CONVECTIONVOP_HPP
#define PIMPACT_CONVECTIONVOP_HPP


#include "Pimpact_ConvectionVWrap.hpp"
#include "Pimpact_ConvectionField.hpp"



namespace Pimpact {



/// \brief Convection Operator for Velocity fields
/// \todo make wind template parameter as well.
/// \ingroup BaseOperator
/// \relates ConvectionSOp
template<class ConvVWrapT>
class ConvectionVOp {

public:

  typedef typename ConvVWrapT::SpaceT SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  static const int dimension = SpaceT::dimension;

  static const int dimNC = SpaceT::dimNC;

  typedef VectorField<SpaceT>  DomainFieldT;
  typedef VectorField<SpaceT>  RangeFieldT;

protected:

  Teuchos::RCP<const ConvVWrapT> convVWrap_;

  Teuchos::RCP< ConvectionField<SpaceT> > convField_;


public:

//  ConvectionVOp(
//      const Teuchos::RCP<const SpaceT>& space ):
//        convVWrap_( create<ConvVWrapT>( space ) ),
//        convField_( create<ConvectionField>( space ) ) {};

  ConvectionVOp(
      const Teuchos::RCP<const ConvVWrapT>& convVWrap ):
        convVWrap_( convVWrap ),
        convField_( create<ConvectionField>( convVWrap->space() ) ) {};


  void assignField( const DomainFieldT& mv ) const {

    convField_->assignField( mv );

  };

  /// \note Operator's wind has to be assigned correctly
  /// \deprecated
  void apply( const DomainFieldT& x, RangeFieldT& y, Scalar mul=0. ) const {

    if( mul<1.e-12 ) {
      y.init(0);
      mul=1.;
    }

    convVWrap_->apply( convField_->get(), x, y, mul );

  }

  void apply(const DomainFieldT& z, const DomainFieldT& x, RangeFieldT& y, Scalar mul=0. ) const {};


  bool hasApplyTranspose() const { return( false ); }


}; // end of class ConvectionVOp



/// \relates ConvectionVOp
template<class SpaceT>
Teuchos::RCP<ConvectionVOp<ConvectionVWrap<ConvectionSOp<SpaceT> > > > createConvectionVOp(
    const Teuchos::RCP<const SpaceT>& space ) {

  auto sop = Pimpact::create<Pimpact::ConvectionSOp>( space ) ;
  auto wrap = Pimpact::create<Pimpact::ConvectionVWrap>( sop );

  return( Teuchos::rcp( new ConvectionVOp<ConvectionVWrap<ConvectionSOp<SpaceT> > >( wrap ) ) );

}





} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONVOP_HPP
