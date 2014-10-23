#pragma once
#ifndef PIMPACT_MULTINONLINEAROP_HPP
#define PIMPACT_MULTINONLINEAROP_HPP


#include "Pimpact_Types.hpp"
#include "Pimpact_FieldFactory.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiHarmonicField.hpp"

#include "Pimpact_ConvectionVOp.hpp"



namespace Pimpact {



/// \ingroup MultiHarmonicOperator
template<class Scalar,class Ordinal>
class MultiHarmonicNonlinear {

  typedef Scalar S;
  typedef Ordinal O;

public:

  typedef MultiHarmonicField< VectorField<Scalar,Ordinal> >  DomainFieldT;
  typedef MultiHarmonicField< VectorField<Scalar,Ordinal> >  RangeFieldT;

protected:

  Teuchos::RCP<ConvectionOp<S,O> > op;

public:

  MultiHarmonicNonlinear(
      const Teuchos::RCP<const Space<S,O,3> >& space
      ):
    op( Teuchos::rcp( new ConvectionOp<S,O>( space ) ) ) {};

  void assignField( const DomainFieldT& mv ) {
  };

protected:

  void apply(const DomainFieldT& x, const DomainFieldT& y, RangeFieldT& z, bool init_yes=true ) const {

    int Nf = x.getNumberModes();
    if( init_yes )
      z.init( 0. );

    // computing zero mode of y
    op->apply( x.getConst0Field(), y.getConst0Field(), z.get0Field(), 1.);

    for( int i=1; i<=Nf; ++i ) {
      op->apply( x.getConstCField(i-1), y.getConstCField(i-1), z.get0Field(), 0.5 );
      op->apply( x.getConstSField(i-1), y.getConstSField(i-1), z.get0Field(), 0.5 );
    }


    // computing cos mode of y
    for( int i=1; i<=Nf; ++i ) {
      op->apply( x.getConst0Field(), y.getConstCField(i-1), z.getCField(i-1), 1. );

      op->apply( x.getConstCField(i-1), y.getConst0Field(), z.getCField(i-1), 1. );

      for( int k=1; k+i<=Nf; ++k ) {
        op->apply( x.getConstCField(k+i-1), y.getConstCField(k-1), z.getCField(i-1), 0.5 );

        op->apply( x.getConstCField(k-1), y.getConstCField(k+i-1), z.getCField(i-1), 0.5 );

        op->apply( x.getConstSField(k+i-1), y.getConstSField(k-1), z.getCField(i-1), 0.5 );

        op->apply( x.getConstSField(k-1), y.getConstSField(k+i-1), z.getCField(i-1), 0.5 );
      }
    }

    // computing sin mode of y
    for( int i=1; i<=Nf; ++i ) {
      op->apply( x.getConst0Field(), y.getConstSField(i-1), z.getSField(i-1), 1. );

      op->apply( x.getConstSField(i-1), y.getConst0Field(), z.getSField(i-1), 1. );

      for( int k=1; k+i<=Nf; ++k ) {
        op->apply( x.getConstCField(k+i-1), y.getConstSField(k-1), z.getSField(i-1), -0.5 );

        op->apply( x.getConstCField(k-1), y.getConstSField(k+i-1), z.getSField(i-1), 0.5 );

        op->apply( x.getConstSField(k+i-1), y.getConstCField(k-1), z.getSField(i-1), 0.5 );

        op->apply( x.getConstSField(k-1), y.getConstCField(k+i-1), z.getSField(i-1), -0.5 );
      }
    }

    // strange terms
    int i;
    for( int k=1; k<=Nf; ++k ) {
      for( int l=1; l<=Nf; ++l ) {
        i = k+l;
        if( i<=Nf ) {
          op->apply( x.getConstCField(k-1), y.getConstCField(l-1), z.getCField(i-1), 0.5 );
          op->apply( x.getConstSField(k-1), y.getConstSField(l-1), z.getCField(i-1), -0.5 );

          op->apply( x.getConstCField(k-1), y.getConstSField(l-1), z.getSField(i-1), 0.5 );
          op->apply( x.getConstSField(k-1), y.getConstCField(l-1), z.getSField(i-1), 0.5 );
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
    const Teuchos::RCP<const Space<S,O,3> >& space
//    const Teuchos::RCP< VectorField<S,O> >& u = Teuchos::null
    ) {
  return( Teuchos::rcp( new MultiHarmonicNonlinear<S,O>( space /*u*/ ) ) );
}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_MULTINONLINEAROP_HPP
