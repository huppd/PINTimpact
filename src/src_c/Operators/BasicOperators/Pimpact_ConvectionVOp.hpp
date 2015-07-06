#pragma once
#ifndef PIMPACT_CONVECTIONVOP_HPP
#define PIMPACT_CONVECTIONVOP_HPP


#include "Pimpact_ConvectionVWrap.hpp"
#include "Pimpact_ConvectionField.hpp"



namespace Pimpact {



/// \brief Convection Operator for Velocity fields
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

  Teuchos::RCP<ConvectionVWrap<ConvSOpT> > convVWrap_;

  Teuchos::RCP<ConvectionField<SpaceT> > convField_;

public:

    ConvectionVOp( const Teuchos::RCP<const SpaceT>& space ):
      convVWrap_( create<ConvectionVWrap<ConvSOpT> >( create<ConvSOpT>(space) ) ),
      convField_( create<ConvectionField>( space ) ) {};

    template< class ConvSOpTT >
    ConvectionVOp( const Teuchos::RCP<ConvSOpTT>& op ):
      convVWrap_( create<ConvectionVWrap<ConvSOpT> >( create<ConvSOpT>( op->space() ) ) ),
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


  Teuchos::RCP<const SpaceT> space() const {
    return( convVWrap_->space() );
  }

  Teuchos::RCP< ConvectionField<SpaceT> >
  getConvField() const {
    return( convField_ );
  }

  Teuchos::RCP<const ConvSOpT>
  getSOp() const {
    return( convVWrap_->getSOp() );
  }


	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		convVWrap_->setParameter( para );
	}


  bool hasApplyTranspose() const { return( false ); }

  void print( std::ostream& out=std::cout ) const {
    convVWrap_->print(out);
  }

	const std::string getLabel() const { return( convVWrap_->getSOp()->getLabel() + "VOp " ); };

}; // end of class ConvectionVOp



/// \relates ConvectionVOp
template<class SpaceT>
Teuchos::RCP<ConvectionVOp<ConvectionVWrap<ConvectionSOp<SpaceT> > > > createConvectionVOp(
    const Teuchos::RCP<const SpaceT>& space ) {

  auto sop = Pimpact::create<Pimpact::ConvectionSOp>( space ) ;
//  auto wrap = Pimpact::create<Pimpact::ConvectionVWrap>( sop );

  return( Teuchos::rcp( new ConvectionVOp< ConvectionSOp<SpaceT> >( sop ) ) );


}





} // end of namespace Pimpact


#ifdef COMPILE_ETI
#include "Pimpact_ConvectionDiffusionSOp.hpp"
extern template class Pimpact::ConvectionVOp< Pimpact::ConvectionSOp< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::ConvectionVOp< Pimpact::ConvectionSOp< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::ConvectionVOp< Pimpact::ConvectionSOp< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::ConvectionVOp< Pimpact::ConvectionSOp< Pimpact::Space<double,int,4,4> > >;
extern template class Pimpact::ConvectionVOp< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::ConvectionVOp< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::ConvectionVOp< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::ConvectionVOp< Pimpact::ConvectionDiffusionSOp< Pimpact::Space<double,int,4,4> > >;
#endif


#endif // end of #ifndef PIMPACT_CONVECTIONVOP_HPP
