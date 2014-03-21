#pragma once
#ifndef PIMPACT_MULTINONLINEAROP_HPP
#define PIMPACT_MULTINONLINEAROP_HPP


#include "Pimpact_Types.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"

#include "Pimpact_NonlinearOp.hpp"



namespace Pimpact {



/// \ingroup MultiHarmonicOperator
template<class Scalar,class Ordinal>
class MultiHarmonicNonlinear {

  typedef Scalar S;
  typedef Ordinal O;

public:

  typedef MultiHarmonicField< VectorField<Scalar,Ordinal> >  DomainFieldT;
  typedef MultiHarmonicField< VectorField<Scalar,Ordinal> >  RangeFieldT;
  typedef NonModeOp OpType;

protected:

  Teuchos::RCP<VectorField<S,O> > temp_;
  Teuchos::RCP<Nonlinear<S,O> > op;

public:

  MultiHarmonicNonlinear( const Teuchos::RCP<VectorField<S,O> >& temp=Teuchos::null ):
    temp_(temp),
    op( Teuchos::rcp( new Nonlinear<S,O>() )) {};

  void assignField( const DomainFieldT& mv ) {
  };

  void apply(const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z) const {
    int Nf = x.getNumberModes();
    z.init( 0. );

    // computing zero mode of y
    op->apply( x.getConst0Field(), y.getConst0Field(), *temp_ );
    z.get0Field().add( 1., z.get0Field(), 1., *temp_ );

    for( int i=0; i<Nf; ++i ) {
      op->apply( x.getConstCField(i), y.getConstCField(i), *temp_ );
      z.get0Field().add( 1., z.get0Field(), 1., *temp_ );
    }

    for( int i=0; i<Nf; ++i ) {
      op->apply( x.getConstSField(i), y.getConstSField(i), *temp_ );
      z.get0Field().add( 1., z.get0Field(), 1., *temp_ );
    }


    // computing cos mode of y
    for( int m=0; m<Nf; ++m ) {
      op->apply( x.getConst0Field(), y.getConstCField(m), *temp_ );
      z.getCField(m).add( 1., z.getCField(m), 1., *temp_ );

      op->apply( x.getConstCField(m), y.getConst0Field(), *temp_ );
      z.getCField(m).add( 1., z.getCField(m), 1., *temp_ );

      for( int k=0; m+k<Nf; ++k ) {
        op->apply( x.getConstCField(m+k), y.getConstCField(k), *temp_ );
        z.getCField(m).add( 1., z.getCField(m), 0.5, *temp_ );

        op->apply( x.getConstCField(k), y.getConstCField(m+k), *temp_ );
        z.getCField(m).add( 1., z.getCField(m), 0.5, *temp_ );

        op->apply( x.getConstSField(m+k), y.getConstSField(k), *temp_ );
        z.getCField(m).add( 1., z.getCField(m), 0.5, *temp_ );

        op->apply( x.getConstSField(k), y.getConstSField(m+k), *temp_ );
        z.getCField(m).add( 1., z.getCField(m), 0.5, *temp_ );
      }
    }

    // computing sin mode of y
    for( int m=0; m<Nf; ++m ) {
      op->apply( x.getConst0Field(), y.getConstSField(m), *temp_ );
      z.getSField(m).add( 1., z.getSField(m), 1., *temp_ );

      op->apply( x.getConstSField(m), y.getConst0Field(), *temp_ );
      z.getSField(m).add( 1., z.getSField(m), 1., *temp_ );

      for( int k=0; k+m<Nf; ++k ) {
        op->apply( x.getConstCField(m+k), y.getConstSField(k), *temp_ );
        z.getSField(m).add( 1., z.getSField(m), -0.5, *temp_ );

        op->apply( x.getConstCField(k), y.getConstSField(m+k), *temp_ );
        z.getSField(m).add( 1., z.getSField(m), 0.5, *temp_ );

        op->apply( x.getConstSField(m+k), y.getConstCField(k), *temp_ );
        z.getSField(m).add( 1., z.getSField(m), 0.5, *temp_ );

        op->apply( x.getConstSField(k), y.getConstCField(m+k), *temp_ );
        z.getSField(m).add( 1., z.getSField(m), -0.5, *temp_ );
      }
    }

    // strange terms
    int m;
    for( int k=0; k<Nf; ++k ) {
      for( int l=0; l<Nf; ++l ) {
        m = k+l;
        if( m<Nf ) {
          op->apply( x.getConstCField(k), y.getConstCField(l), *temp_ );
          z.getCField(m).add( 1., z.getCField(m), 0.5, *temp_ );
          op->apply( x.getConstSField(k), y.getConstSField(l), *temp_ );
          z.getCField(m).add( 1., z.getCField(m), -0.5, *temp_ );

          op->apply( x.getConstCField(k), y.getConstSField(l), *temp_ );
          z.getSField(m).add( 1., z.getSField(m), 0.5, *temp_ );
          op->apply( x.getConstSField(k), y.getConstCField(l), *temp_ );
          z.getSField(m).add( 1., z.getSField(m), -0.5, *temp_ );
        }
      }
    }
  }


  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    apply(x,x,y);
  }


  bool hasApplyTranspose() const { return( false ); }

}; // end of class MultiHarmonicNonlinearOp



/// \relates MultiHarmonicNonlinearOp
template< class S, class O>
Teuchos::RCP<MultiHarmonicNonlinear<S,O> > createMultiHarmonicNonlinear(
    const Teuchos::RCP< VectorField<S,O> >& u = Teuchos::null ) {
  return( Teuchos::rcp( new MultiHarmonicNonlinear<S,O>( u ) ) );
}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULTINONLINEAROP_HPP
