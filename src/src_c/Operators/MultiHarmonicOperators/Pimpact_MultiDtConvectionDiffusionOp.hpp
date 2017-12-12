#pragma once
#ifndef PIMPACT_MULTIDTCONVECTIONDIFFUSIONOP_HPP
#define PIMPACT_MULTIDTCONVECTIONDIFFUSIONOP_HPP


#include "Pimpact_ConvectionField.hpp"
#include "Pimpact_ConvectionDiffusionSOp.hpp"
#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"
#include "Pimpact_NonlinearVWrap.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \ingroup MultiHarmonicOperator
template<class ST>
class MultiDtConvectionDiffusionOp {

public:

  using SpaceT = ST;

  using DomainFieldT = MultiHarmonicField< VectorField<SpaceT> >;
  using RangeFieldT = MultiHarmonicField< VectorField<SpaceT> >;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  Teuchos::RCP< NonlinearWrap< ConvectionDiffusionSOp<SpaceT> > > op_;

  Teuchos::RCP< ConvectionField<SpaceT> > wind0_;
  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > windc_;
  Teuchos::Array< Teuchos::RCP<ConvectionField<SpaceT> > > winds_;

  using FieldTensorT = typename ConvectionField<SpaceT>::FieldTensor;

public:

  MultiDtConvectionDiffusionOp( const Teuchos::RCP<const SpaceT>& space ):
    op_( create<NonlinearWrap>( create<ConvectionDiffusionSOp<SpaceT> >(space) ) ),
    wind0_( create<ConvectionField>(space) ),
    windc_( space->nGlo(3) ),
    winds_( space->nGlo(3) ) {

    for( Ordinal i=0; i<space->nGlo(3); ++i ) {
      windc_[i] = create<ConvectionField>( space );
      winds_[i] = create<ConvectionField>( space );
    }
  };


  void assignField( const DomainFieldT& y_ref ) {

    Teuchos::RCP<const DomainFieldT> y;

    if( y_ref.global() )
      y = Teuchos::rcpFromRef( y_ref );
    else {
      Teuchos::RCP<DomainFieldT> temp = Teuchos::rcp( new DomainFieldT( space(), true ) );
      *temp = y_ref;
      y = temp;
      //std::cout << "assign op: y->global(): " << y->global() << "\n";
    }

    y->exchange();

    wind0_->assignField( y->get0Field() );

    for( Ordinal i=1; i<=space()->nGlo(3); ++i ) {
      windc_[i-1]->assignField( y->getCField(i) );
      winds_[i-1]->assignField( y->getSField(i) );
    }
  };


  void apply( const DomainFieldT& y_ref, RangeFieldT& z, bool init_yes=true ) const {

    Teuchos::RCP<const DomainFieldT> y;
    if( y_ref.global() )
      y = Teuchos::rcpFromRef( y_ref );
    else {
      Teuchos::RCP<DomainFieldT> temp = Teuchos::rcp( new DomainFieldT( space(), true ) );
      *temp = y_ref; // needed because of const
      y = temp;
    }

    y->exchange();

    Ordinal Nf = space()->nGlo(3);
    Scalar iRe = 1./op_->space()->getDomainSize()->getRe();
    Scalar a2 = op_->space()->getDomainSize()->getAlpha2()*iRe;

    Scalar mulI;

    // computing zero mode of z
    if( 0==space()->si(F::U,3) ) {

      op_->apply( get0Wind(), y->get0Field(), z.get0Field(), 0., 1., iRe, Add::N );

      for( Ordinal i=1; i<=Nf; ++i ) {
        op_->apply( getCWind(i), y->getCField(i), z.get0Field(), 0., 0.5, 0., Add::Y );
        op_->apply( getSWind(i), y->getSField(i), z.get0Field(), 0., 0.5, 0., Add::Y );
      }
    }

    // computing cos mode of z
    for( Ordinal i=std::max(space()->si(F::U,3),1); i<=space()->ei(F::U,3); ++i ) {

      op_->apply( get0Wind( ), y->getCField(i), z.getCField(i), 0., 1., iRe, Add::N  );
      op_->apply( getCWind(i), y->get0Field( ), z.getCField(i), 0., 1., 0.,  Add::Y );

      for( Ordinal k=1; k+i<=Nf; ++k ) { // thats fine

        mulI = (k==i)?(a2*i):0;

        op_->apply( getCWind(k+i), y->getCField( k ), z.getCField(i),   0., 0.5, 0., Add::Y );
        op_->apply( getCWind( k ), y->getCField(k+i), z.getCField(i),   0., 0.5, 0., Add::Y );
        op_->apply( getSWind(k+i), y->getSField( k ), z.getCField(i), mulI, 0.5, 0., Add::Y );
        op_->apply( getSWind( k ), y->getSField(k+i), z.getCField(i),   0., 0.5, 0., Add::Y );
      }
    }

    // computing sin mode of y
    for( Ordinal i=std::max(space()->si(F::U,3),1); i<=space()->ei(F::U,3); ++i ) {

      op_->apply( get0Wind(),  y->getSField(i), z.getSField(i), 0., 1., iRe, Add::N  );
      op_->apply( getSWind(i), y->get0Field(),  z.getSField(i), 0., 1., 0. , Add::Y );

      for( Ordinal k=1; k+i<=Nf; ++k ) { // that is fine

        mulI = (k==i)?(a2*i):0;

        op_->apply( getCWind(k+i), y->getSField( k ), z.getSField(i),    0., -0.5, 0., Add::Y );
        op_->apply( getCWind( k ), y->getSField(k+i), z.getSField(i),    0.,  0.5, 0., Add::Y );
        op_->apply( getSWind(k+i), y->getCField( k ), z.getSField(i), -mulI,  0.5, 0., Add::Y );
        op_->apply( getSWind( k ), y->getCField(k+i), z.getSField(i),    0., -0.5, 0., Add::Y );
      }
    }

    // rest of time
    for( Ordinal i=std::max(space()->si(F::U,3),1); i<=space()->ei(F::U,3); ++i ) {
      if( Nf/2+1<=i && i<=Nf ) {
        mulI = a2*i;
        z.getCField(i).add( 1., z.getCField(i),  mulI, y->getSField(i), B::N );
        z.getSField(i).add( 1., z.getSField(i), -mulI, y->getCField(i), B::N );
      }
    }

    // strange terms
    Ordinal i;
    for( Ordinal k=1; k<=Nf; ++k ) {
      for( Ordinal l=1; l<=Nf; ++l ) { // that is fine
        i = k+l;
        if( i<=Nf ) { // do something here
          if( std::max(space()->si(F::U,3),1)<=i && i<=space()->ei(F::U,3) ) {
            op_->apply( getCWind(k), y->getCField(l), z.getCField(i), 0.,  0.5, 0., Add::Y );
            op_->apply( getSWind(k), y->getSField(l), z.getCField(i), 0., -0.5, 0., Add::Y );

            op_->apply( getCWind(k), y->getSField(l), z.getSField(i), 0.,  0.5, 0., Add::Y );
            op_->apply( getSWind(k), y->getCField(l), z.getSField(i), 0.,  0.5, 0., Add::Y );
          }
        }
      }
    }

    z.changed();
  }


  ST compRefRes( const DomainFieldT& y_ref, bool init_yes=true ) const {

    Teuchos::RCP<const DomainFieldT> y;
    if( y_ref.global() )
      y = Teuchos::rcpFromRef( y_ref );
    else {
      Teuchos::RCP<DomainFieldT> temp = Teuchos::rcp( new DomainFieldT( space(), true ) );
      *temp = y_ref; // needed because of const
      y = temp;
    }

    y->exchange();

    auto z = y.getField(1).clone();

    Ordinal Nf = space()->nGlo(3);
    Scalar iRe = 1./op_->space()->getDomainSize()->getRe();
    Scalar a2 = op_->space()->getDomainSize()->getAlpha2()*iRe;

    Scalar mulI;

    // computing cos mode of z
    Ordinal i = Nf+1;

    for( Ordinal k=1; k+i<=Nf; ++k ) { // thats fine

      op_->apply( getCWind(k+i), y->getCField( k ), z.getCField(), 0., 0.5, 0., Add::Y );
      op_->apply( getCWind( k ), y->getCField(k+i), z.getCField(), 0., 0.5, 0., Add::Y );
      op_->apply( getSWind(k+i), y->getSField( k ), z.getCField(), 0., 0.5, 0., Add::Y );
      op_->apply( getSWind( k ), y->getSField(k+i), z.getCField(), 0., 0.5, 0., Add::Y );
    }

    // computing sin mode of y

    op_->apply( get0Wind(),  y->getSField(i), z.getSField(i), 0., 1., iRe, Add::N  );
    op_->apply( getSWind(i), y->get0Field(),  z.getSField(i), 0., 1., 0. , Add::Y );

    for( Ordinal k=1; k+i<=Nf; ++k ) { // that is fine

      op_->apply( getCWind(k+i), y->getSField( k ), z.getSField(), 0., -0.5, 0., Add::Y );
      op_->apply( getCWind( k ), y->getSField(k+i), z.getSField(), 0.,  0.5, 0., Add::Y );
      op_->apply( getSWind(k+i), y->getCField( k ), z.getSField(), 0.,  0.5, 0., Add::Y );
      op_->apply( getSWind( k ), y->getCField(k+i), z.getSField(), 0., -0.5, 0., Add::Y );
    }

    // strange terms
    for( Ordinal k=1; k<=Nf; ++k ) {
      for( Ordinal l=1; l<=Nf; ++l ) { // that is fine
        i = k+l;
        if( i==Nf+1 ) { // do something here
          if( std::max(space()->si(F::U,3),1)<=i && i<=space()->ei(F::U,3) ) {
            op_->apply( getCWind(k), y->getCField(l), z.getCField(), 0.,  0.5, 0., Add::Y );
            op_->apply( getSWind(k), y->getSField(l), z.getCField(), 0., -0.5, 0., Add::Y );

            op_->apply( getCWind(k), y->getSField(l), z.getSField(), 0.,  0.5, 0., Add::Y );
            op_->apply( getSWind(k), y->getCField(l), z.getSField(), 0.,  0.5, 0., Add::Y );
          }
        }
      }
    }

    return z->norm( ENorm::L2 );
    //z.changed();
  }

  void applyBC( const DomainFieldT& x, RangeFieldT& y ) const {

    if( 0==space()->si(F::U,3) )
      op_->getSOp()->getHelmholtzOp()->applyBC( x.get0Field(), y.get0Field() );

    for( typename SpaceT::Ordinal i=std::max(space()->si(F::U,3),1); i<=space()->ei(F::U,3); ++i ) {
      op_->getSOp()->getHelmholtzOp()->applyBC( x.getCField(i), y.getCField(i) );
      op_->getSOp()->getHelmholtzOp()->applyBC( x.getSField(i), y.getSField(i) );
    }
  }


  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return op_->space();
  };

  void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const {
    return false;
  }

  const std::string getLabel() const {
    return "MHDtConvectionDiffusion";
  };

  void print( std::ostream& out=std::cout ) const {
    out <<  getLabel() << ":\n";
    op_->print( out );
  }

protected:

  constexpr const FieldTensorT& get0Wind() const {
    return wind0_->get();
  }

  constexpr const FieldTensorT& getCWind( const Ordinal i) const {
    return windc_[i-1]->get();
  }
  constexpr const FieldTensorT& getSWind( const Ordinal i) const {
    return winds_[i-1]->get();
  }

}; // end of class MultiDtConvectionDiffusionOp



/// \relates MultiDtConvectionDiffusionOp
template<class SpaceT>
Teuchos::RCP<MultiDtConvectionDiffusionOp<SpaceT> >
createMultiDtConvectionDiffusionOp( const Teuchos::RCP<const SpaceT>& space ) {

  return Teuchos::rcp( new MultiDtConvectionDiffusionOp<SpaceT>( space ) );
}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULTIDTCONVECTIONDIFFUSIONOP_HPP
