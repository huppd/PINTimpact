#pragma once
#ifndef PIMPACT_MODEPREC_HPP
#define PIMPACT_MODEPREC_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_ModeField.hpp"




namespace Pimpact {



/// \ingroup ModeOperator
template<class OpT>
class ModePrec {

public:

  using SpaceT = typename OpT::SpaceT;

  using DomainFieldT = ModeField<typename OpT::DomainFieldT>;
  using RangeFieldT = ModeField<typename OpT::RangeFieldT>;

protected:

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

  ST mulI_;
  ST mulC_;
  ST mulL_;

  int type_;

  Teuchos::RCP<OpT> op_;

public:

  ModePrec(
    const Teuchos::RCP<OpT>& op,
    const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList() ):

    mulI_(0.),
    mulC_(1.),
    mulL_( 1./op->space()->getDomainSize()->getRe() ),
    type_( pl->get<int>("type", -1) ),
    op_(op) {};

  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    //std:: cout << type_ << "\n";

    switch(type_) {
      case 0: {
        y = x;
        break;
      }
      case 1: {
        applyDTinv( x, y );
        break;
      }
      case 2: {
        applyCDinv( x, y );
        break;
      }
      case 3: {
        applyELinv( x, y );
        break;
      }
      case 4: {
        applyERinv( x, y );
        break;
      }
      case 5: {
        applyER2inv( x, y );
        break;
      }
      default: {
        if( mulI_>=mulC_ && mulI_>=mulL_ )
          applyERinv( x, y );
        else
          applyCDinv( x, y );
        break;
      }
    }
  }

  /// left/right same
  /// \todo fix BC
  /// \deprecated
  void applyDTinv( const DomainFieldT& x, RangeFieldT& y ) const {
    //std::cout << "applyDTinv\n";
    //y.getCField().add( 0.0,      x.getCField(), -1.0/mulI_, x.getSField(), B::N );
    //y.getSField().add( 1.0/mulI_, x.getCField(), 0.0,      x.getSField(), B::N );
  }

  void applyCDinv( const DomainFieldT& x, RangeFieldT& y ) const {
    //std::cout << "applyCDinv\n";

    //set paramters
    auto pl = Teuchos::parameterList();
    pl->set<ST>( "mulI", 0. );
    pl->set<ST>( "mulC", mulC_ );
    pl->set<ST>( "mulL", mulL_ );
    op_->setParameter( pl );

    op_->apply( x.getCField(), y.getCField() );
    op_->apply( x.getSField(), y.getSField() );
  }

  void applyELinv( const DomainFieldT& x, RangeFieldT& y ) const {

    //std::cout << "applyELinv\n";
    DomainFieldT temp( space() );

    // left
    temp = x;
    temp.getCField().add( 1.0, x.getCField(),  1.0, x.getSField(), B::N );
    temp.getSField().add( 1.0, x.getCField(), -1.0, x.getSField(), B::N );

    //// set paramters
    auto pl = Teuchos::parameterList();

    pl->set<ST>( "mulI", mulI_ );
    pl->set<ST>( "mulC", mulC_ );
    pl->set<ST>( "mulL", mulL_ );

    op_->setParameter( pl );

    op_->apply( temp.getCField(), y.getCField() );
    op_->apply( temp.getSField(), y.getSField() );

    y.scale( 0.5, B::N );
  }


  void applyERinv( const DomainFieldT& x, RangeFieldT& y ) const {

    //std::cout << "applyERinv\n";
    // right
    // set paramters
    auto pl = Teuchos::parameterList();
    pl->set<ST>( "mulI", mulI_ );
    pl->set<ST>( "mulC", mulC_ );
    pl->set<ST>( "mulL", mulL_ );

    op_->setParameter( pl );

    op_->apply( x.getCField(), y.getCField() );
    op_->apply( x.getSField(), y.getSField() );

    const B wb=B::N;

    for( F f=F::U; f<SpaceT::sdim; ++f )
      for( OT k=space()->si(f,Z,wb); k<=space()->ei(f,Z,wb); ++k )
        for( OT j=space()->si(f,Y,wb); j<=space()->ei(f,Y,wb); ++j )
          for( OT i=space()->si(f,X,wb); i<=space()->ei(f,X,wb); ++i ) {

            ST tempC =  0.5*y.getCField()(f)(i,j,k) + 0.5*y.getSField()(f)(i,j,k);
            ST tempS =  0.5*y.getCField()(f)(i,j,k) - 0.5*y.getSField()(f)(i,j,k);
            y.getCField()(f)(i,j,k) = tempC;
            y.getSField()(f)(i,j,k) = tempS;
          }
  }


  void applyER2inv( const DomainFieldT& x, RangeFieldT& y ) const {

    //std::cout << "applyERinv\n";
    // right
    // set paramters
    auto pl = Teuchos::parameterList();
    pl->set<ST>( "mulI", mulI_ );
    pl->set<ST>( "mulC", mulC_ );
    pl->set<ST>( "mulL", mulL_ );

    op_->setParameter( pl );

    op_->apply( x.getCField(), y.getCField() );

    pl->set<ST>( "mulI", -mulI_ );
    pl->set<ST>( "mulC", -mulC_ );
    pl->set<ST>( "mulL", -mulL_ );

    op_->setParameter( pl );

    op_->apply( x.getSField(), y.getSField() );
  }


  void assignField(const DomainFieldT& mv) {};

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return op_->space();
  };

  Teuchos::RCP<OpT> getOperator() const {
    return op_;
  };

  void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {

    if( para->name()!="Linear Solver" ) {
      mulI_ = para->get<ST>( "mulI" );
      mulC_ = para->get<ST>( "mulC" );
      mulL_ = para->get<ST>( "mulL" );
    }
  }

  bool hasApplyTranspose() const {
    return false;
  }

  const std::string getLabel() const {
    return "ModePrec";
  };

  void print( std::ostream& out=std::cout ) const {
    out << getLabel() << ":\n";
    op_->print( out );
  }

}; // end of class ModePrec


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MODEPREC_HPP
