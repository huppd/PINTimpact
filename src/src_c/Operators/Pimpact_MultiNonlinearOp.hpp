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
//  typedef NonModeOp OpType;

protected:

  Teuchos::RCP<VectorField<S,O> > temp_;
  Teuchos::RCP<Nonlinear<S,O> > op;

public:

  MultiHarmonicNonlinear( const Teuchos::RCP<VectorField<S,O> >& temp=Teuchos::null ):
    temp_(temp),
    op( Teuchos::rcp( new Nonlinear<S,O>() )) {};

  void assignField( const DomainFieldT& mv ) {
  };

protected:

  void apply(const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z) const {
    int Nf = x.getNumberModes();
    z.init( 0. );

    // computing zero mode of y
    op->apply( x.getConst0Field(), y.getConst0Field(), *temp_ );
    z.get0Field().add( 1., z.getConst0Field(), 1., *temp_ );

    for( int i=1; i<=Nf; ++i ) {
      op->apply( x.getConstCField(i-1), y.getConstCField(i-1), *temp_ );
      z.get0Field().add( 1., z.getConst0Field(), 0.5, *temp_ );
    }

    for( int i=1; i<=Nf; ++i ) {
      op->apply( x.getConstSField(i-1), y.getConstSField(i-1), *temp_ );
      z.get0Field().add( 1., z.getConst0Field(), 0.5, *temp_ );
    }


    // computing cos mode of y
    for( int i=1; i<=Nf; ++i ) {
      op->apply( x.getConst0Field(), y.getConstCField(i-1), *temp_ );
      z.getCField(i-1).add( 1., z.getConstCField(i-1), 1., *temp_ );

      op->apply( x.getConstCField(i-1), y.getConst0Field(), *temp_ );
      z.getCField(i-1).add( 1., z.getConstCField(i-1), 1., *temp_ );

      for( int k=1; k+i<=Nf; ++k ) {
        op->apply( x.getConstCField(k+i-1), y.getConstCField(k-1), *temp_ );
        z.getCField(i-1).add( 1., z.getConstCField(i-1), 0.5, *temp_ );

        op->apply( x.getConstCField(k-1), y.getConstCField(k+i-1), *temp_ );
        z.getCField(i-1).add( 1., z.getCField(i-1), 0.5, *temp_ );

        op->apply( x.getConstSField(k+i-1), y.getConstSField(k-1), *temp_ );
        z.getCField(i-1).add( 1., z.getConstCField(i-1), 0.5, *temp_ );

        op->apply( x.getConstSField(k-1), y.getConstSField(k+i-1), *temp_ );
        z.getCField(i-1).add( 1., z.getConstCField(i-1), 0.5, *temp_ );
      }
    }

    // computing sin mode of y
    for( int i=1; i<=Nf; ++i ) {
      op->apply( x.getConst0Field(), y.getConstSField(i-1), *temp_ );
      z.getSField(i-1).add( 1., z.getConstSField(i-1), 1., *temp_ );

      op->apply( x.getConstSField(i-1), y.getConst0Field(), *temp_ );
      z.getSField(i-1).add( 1., z.getConstSField(i-1), 1., *temp_ );

      for( int k=1; k+i<=Nf; ++k ) {
        op->apply( x.getConstCField(k+i-1), y.getConstSField(k-1), *temp_ );
        z.getSField(i-1).add( 1., z.getConstSField(i-1), -0.5, *temp_ );

        op->apply( x.getConstCField(k-1), y.getConstSField(k+i-1), *temp_ );
        z.getSField(i-1).add( 1., z.getConstSField(i-1), 0.5, *temp_ );

        op->apply( x.getConstSField(k+i-1), y.getConstCField(k-1), *temp_ );
        z.getSField(i-1).add( 1., z.getConstSField(i-1), 0.5, *temp_ );

        op->apply( x.getConstSField(k-1), y.getConstCField(k+i-1), *temp_ );
        z.getSField(i-1).add( 1., z.getConstSField(i-1), -0.5, *temp_ );
      }
    }

    // strange terms
    int i;
    for( int k=1; k<=Nf; ++k ) {
      for( int l=1; l<=Nf; ++l ) {
        i = k+l;
        if( i<=Nf ) {
          op->apply( x.getConstCField(k-1), y.getConstCField(l-1), *temp_ );
          z.getCField(i-1).add( 1., z.getConstCField(i-1), 0.5, *temp_ );
          op->apply( x.getConstSField(k-1), y.getConstSField(l-1), *temp_ );
          z.getCField(i-1).add( 1., z.getConstCField(i-1), -0.5, *temp_ );

          op->apply( x.getConstCField(k-1), y.getConstSField(l-1), *temp_ );
          z.getSField(i-1).add( 1., z.getConstSField(i-1), 0.5, *temp_ );
          op->apply( x.getConstSField(k-1), y.getConstCField(l-1), *temp_ );
          z.getSField(i-1).add( 1., z.getConstSField(i-1), 0.5, *temp_ );
        }
      }
    }
  }

public:

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
