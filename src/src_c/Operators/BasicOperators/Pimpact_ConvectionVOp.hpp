#pragma once
#ifndef PIMPACT_CONVECTIONVOP_HPP
#define PIMPACT_CONVECTIONVOP_HPP


#include "Pimpact_ConvectionVWrap.hpp"
#include "Pimpact_ConvectionField.hpp"



namespace Pimpact {



/// \brief Convection Operator for Velocity fields
/// \todo make wind template parameter as well.(necessary when different winds
/// are wanted, meaning moving interpolation steps from assign to apply
/// \todo make constructor so wind can be shared by different operators.
/// \todo make Smoother
/// \ingroup BaseOperator
/// \ingroup NonlinearOperator
template<class CSOPT>
class ConvectionVOp {

public:

  typedef CSOPT ConvSOpT;

  typedef typename ConvSOpT::SpaceT SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  static const int dimension = SpaceT::dimension;

  static const int dimNC = SpaceT::dimNC;

  typedef VectorField<SpaceT>  DomainFieldT;
  typedef VectorField<SpaceT>  RangeFieldT;

protected:

  Teuchos::RCP<const ConvectionVWrap<ConvSOpT> > convVWrap_;

  Teuchos::RCP< ConvectionField<SpaceT> > convField_;

public:

    ConvectionVOp( const Teuchos::RCP<const SpaceT>& space ):
      convVWrap_( createConst<ConvectionVWrap<ConvSOpT> >( createConst<ConvSOpT>(space) ) ),
      convField_( create<ConvectionField>( space ) ) {};

    template< class ConvSOpTT >
    ConvectionVOp( const Teuchos::RCP<const ConvSOpTT>& op ):
      convVWrap_( createConst<ConvectionVWrap<ConvSOpT> >( createConst<ConvSOpT>( op->getSpace() ) ) ),
      convField_( op->getConvField() ) {}


  void assignField( const DomainFieldT& mv ) const { convField_->assignField( mv ); };


  /// \note Operator's wind has to be assigned correctly
  void apply( const DomainFieldT& x, RangeFieldT& y, Scalar mul=0. ) const {

    if( mul<1.e-12 ) {
      y.init(0);
      mul=1.;
    }

    convVWrap_->apply( convField_->get(), x, y, mul );

  }


  /// \depcecated
  void apply(const DomainFieldT& z, const DomainFieldT& x, RangeFieldT& y, Scalar mul=0. ) const { std::cout << "!!!depcreated!!!\n"; };


  /// \{ \todo unifiy
  Teuchos::RCP<const SpaceT> space() const {
    return( convVWrap_->space() );
  }

  Teuchos::RCP<const SpaceT> getSpace() const {
    return( convVWrap_->space() );
  }
  /// \}

  Teuchos::RCP< ConvectionField<SpaceT> >
  getConvField() const {
    return( convField_ );
  }

  Teuchos::RCP<const ConvSOpT>
  getSOp() const {
    return( convVWrap_->getSOp() );
  }


  bool hasApplyTranspose() const { return( false ); }

  void print( std::ostream& out=std::cout ) const {
    convVWrap_->print(out);
  }


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
