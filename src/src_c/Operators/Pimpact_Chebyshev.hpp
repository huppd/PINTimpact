#pragma once
#ifndef PIMPACT_CHEBYSHEV_HPP
#define PIMPACT_CHEBYSHEV_HPP


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Space.hpp" // just for createOstream<>




namespace Pimpact {



/// \todo make EV computation templated
/// tparam evcomputer<OperatorT>, COTOR( op, pl ), compEV(min,max)
template<class OperatorT>
class Chebyshev {

public:

  using DomainFieldT = typename OperatorT::DomainFieldT;
  using RangeFieldT = typename OperatorT::RangeFieldT;

  using SpaceT = typename DomainFieldT::SpaceT;

protected:

  using Scalar = typename SpaceT::Scalar;

  int numIters_;

  bool zeroStartingSolution_ = false;

  Scalar lamMax_;
  Scalar lamMin_;
  Scalar eigRatio_ = 1./30.;

  Teuchos::RCP<const OperatorT> op_;

  Teuchos::RCP<std::ostream> out_;

public:

  Chebyshev(
    const Teuchos::RCP<const OperatorT>& op,
    const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList() ):
    numIters_( pl->get<int>( "numIters", 12 ) ),
    lamMax_( pl->get<Scalar>( "max EV", Teuchos::ScalarTraits<Scalar>::zero() ) ),
    lamMin_( pl->get<Scalar>( "min EV", Teuchos::ScalarTraits<Scalar>::zero() ) ),
    op_(op) {

    Teuchos::RCP<DomainFieldT> x = create<DomainFieldT>( space() );
    Teuchos::RCP<RangeFieldT> r = create<RangeFieldT>( space() );

    if( pl->get<bool>( "with output", false ) )
      out_ = Pimpact::createOstream( "conv_Cheb.txt" );

    if( std::fabs( lamMax_-lamMin_ )<Teuchos::ScalarTraits<Scalar>::eps() ) {
      x->random();

      Scalar lamp = 0.;
      for( int i=0; i<500; ++i ) {
        op_->apply( *x, *r );
        Scalar lam = r->dot( *x )/x->dot( *x );
        //std::cout << "lambda: " << lam << "\t" << std::abs( lamp-lam )/std::abs(lam) << "\n";
        r->scale( 1./r->norm() );
        x.swap( r );
        if( std::fabs( lamp-lam )/std::fabs(lam) < 1.e-3 )
          break;
        else
          lamp=lam;
      }
      //x->write( 999 );

      op_->apply( *x, *r );
      lamMax_ = r->dot( *x )/x->dot( *x );


      lamMax_ *= 1.1;
      lamMin_ = lamMax_*eigRatio_;
      if( 0==space()->rankST() ) {
        std::cout << "lamMax: (" << lamMax_ << ",\t" << lamMin_ << " )\n";
      }
    }
  }


  void apply( const DomainFieldT& b, RangeFieldT& x ) const {
    //applyPIMP( b, x );
    applyBarrett93( b, x );
    //applyIFFPACK( b, x );
  }

protected:

  //void applyIFFPACK( const DomainFieldT& b, RangeFieldT& x ) const {

//#ifdef HAVE_TEUCHOS_DEBUG
  //using std::cerr;
  //using std::endl;
  //cerr << "\\|B\\|_{\\infty} = " << b.norm( Pimpact::ENorm::Inf ) << endl;
  //cerr << "\\|X\\|_{\\infty} = " << x.norm( Pimpact::ENorm::Inf ) << endl;
//#endif // HAVE_TEUCHOS_DEBUG

  //if( numIters_<=0 ) {
  //return;
  //}

  //const Scalar zero = Teuchos::as<Scalar>(0);
  //const Scalar one = Teuchos::as<Scalar>(1);
  //const Scalar two = Teuchos::as<Scalar>(2);

  //// Initialize coefficients
  //const Scalar alpha = lamMax_ / eigRatio_;
  //const Scalar beta = Teuchos::as<Scalar> (1.1) * lamMax_;
  //const Scalar delta = two / (beta - alpha);
  //const Scalar theta = (beta + alpha) / two;
  //const Scalar s1 = theta * delta;

//#ifdef HAVE_TEUCHOS_DEBUG
//#ifdef IFPACK_DETAILS_CHEBYSHEV_DEBUG
  //cerr << "alpha = " << alpha << endl
  //<< "beta = " << beta << endl
  //<< "delta = " << delta << endl
  //<< "theta = " << theta << endl
  //<< "s1 = " << s1 << endl;
//#endif // IFPACK_DETAILS_CHEBYSHEV_DEBUG
//#endif // HAVE_TEUCHOS_DEBUG

  //// Fetch cached temporary vectors.
  //Teuchos::RCP<DomainFieldT> V_ptr = b.clone( Pimpact::ECopy::Shallow );
  //Teuchos::RCP<DomainFieldT> W_ptr = b.clone( Pimpact::ECopy::Shallow );

  //// mfh 28 Jan 2013: We write V1 instead of V, so as not to confuse
  //// the multivector V with the typedef V (for Tpetra::Vector).
  ////MV V1 (B.getMap (), B.getNumVectors (), false);
  ////MV W (B.getMap (), B.getNumVectors (), false);
  //DomainFieldT& V1 = *V_ptr;
  //DomainFieldT& W = *W_ptr;

//#ifdef HAVE_TEUCHOS_DEBUG
  //cerr << "Iteration " << 1 << ":" << endl
  //<< "- \\|D\\|_{\\infty} = " << D_->norm( Pimpact::ENorm::Inf ) << endl;
//#endif // HAVE_TEUCHOS_DEBUG

  //// Special case for the first iteration.
  //if( !zeroStartingSolution_ ) {
  //op_->computeResidual( b, x, V1 ); // V1 = B - A*X

//#ifdef HAVE_TEUCHOS_DEBUG
  //cerr << "- \\|B - A*X\\|_{\\infty} = " << V1.norm( Pimpact::ENorm::Inf ) << endl;
//#endif // HAVE_TEUCHOS_DEBUG

  ////solve (W, one/theta, D_inv, V1); // W = (1/theta)*D_inv*(B-A*X)
  //W.reciprocal( *D_ );
  //W.scale( V1 );
  //W.scale( one/theta );

//#ifdef HAVE_TEUCHOS_DEBUG
  //cerr << "- \\|W\\|_{\\infty} = " << W.norm( Pimpact::ENorm::Inf ) << endl;
//#endif // HAVE_TEUCHOS_DEBUG

  //x.add( one, x, one, W ); // X = X + W
  //} else {
  ////solve (W, one/theta, D_inv, B); // W = (1/theta)*D_inv*B
  //W.reciprocal( *D_ );
  //W.scale( b );
  //W.scale( one/theta );

//#ifdef HAVE_TEUCHOS_DEBUG
  //cerr << "- \\|W\\|_{\\infty} = " << W.norm( Pimpact::ENorm::Inf ) << endl;
//#endif // HAVE_TEUCHOS_DEBUG

  ////Tpetra::deep_copy(X, W); // X = 0 + W
  //x = W;
  //}
//#ifdef HAVE_TEUCHOS_DEBUG
  //cerr << "- \\|X\\|_{\\infty} = " << x.norm( Pimpact::ENorm::Inf ) << endl;
//#endif // HAVE_TEUCHOS_DEBUG

  //// The rest of the iterations.
  //Scalar rhok = one / s1;
  //Scalar rhokp1, dtemp1, dtemp2;
  //for( int deg=1; deg<numIters_; ++deg) {

//#ifdef HAVE_TEUCHOS_DEBUG
  //cerr << "Iteration " << deg+1 << ":" << endl;
  //cerr << "- \\|D\\|_{\\infty} = " << D_->norm( Pimpact::ENorm::Inf ) << endl;
  //cerr << "- \\|B\\|_{\\infty} = " << b.norm( Pimpact::ENorm::Inf ) << endl;
  ////cerr << "- \\|A\\|_{\\text{frob}} = " << A_->getFrobeniusNorm () << endl;
  //cerr << "- rhok = " << rhok << endl;
  //V1.init( zero );
//#endif // HAVE_TEUCHOS_DEBUG

  ////computeResidual (V1, B, A, X); // V1 = B - A*X
  //op_->computeResidual( b, x, V1 ); // V1 = B - A*X

//#ifdef HAVE_TEUCHOS_DEBUG
  //cerr << "- \\|B - A*X\\|_{\\infty} = " << V1.norm( Pimpact::ENorm::Inf ) << endl;
//#endif // HAVE_TEUCHOS_DEBUG

  //rhokp1 = one / (two * s1 - rhok);
  //dtemp1 = rhokp1 * rhok;
  //dtemp2 = two * rhokp1 * delta;
  //rhok = rhokp1;

//#ifdef HAVE_TEUCHOS_DEBUG
  //cerr << "- dtemp1 = " << dtemp1 << endl
  //<< "- dtemp2 = " << dtemp2 << endl;
//#endif // HAVE_TEUCHOS_DEBUG

  ////W.scale( dtemp1 );
  ////W.elementWiseMultiply( dtemp2, D_inv, V1, one ); // W = one*W + dtemp2*D_inv.*V1
  //{
  //auto temp = W.clone( Pimpact::ECopy::Shallow );
  //temp->reciprocal( *D_ );
  //temp->scale( V1 );
  //W.add( dtemp1, W, dtemp2, *temp );
  //}

  ////X.update (one, W, one);
  //x.add( one, x, one, W ); // X = X + W

//#ifdef HAVE_TEUCHOS_DEBUG
  //cerr << "- \\|W\\|_{\\infty} = " << W.norm( Pimpact::ENorm::Inf ) << endl;
  //cerr << "- \\|X\\|_{\\infty} = " << x.norm( Pimpact::ENorm::Inf ) << endl;
//#endif // HAVE_TEUCHOS_DEBUG

  //}
  //}

  void applyBarrett93( const DomainFieldT& b, RangeFieldT& x0 ) const {

    Scalar d = ( lamMax_ + lamMin_ )/2.;
    //Scalar c = std::abs( lamMax_ - lamMin_ )/2;
    Scalar c = std::abs( lamMax_ - lamMin_ )/2.;

    Scalar alpha = 1./d;
    Scalar beta = 0.;

    Teuchos::RCP<DomainFieldT> x = x0.clone( ECopy::Deep );

    // r_0 = B-Ax_0
    Teuchos::RCP<RangeFieldT> r = x0.clone( ECopy::Shallow );
    Teuchos::RCP<RangeFieldT> z = x0.clone( ECopy::Shallow );
    Teuchos::RCP<RangeFieldT> p = x0.clone( ECopy::Shallow );

    op_->computeResidual( b, *x, *r );

    if( !out_.is_null() )
      *out_ << r->norm() << "\n";

    for( int n=0; n<numIters_; ++n ) {

      // precondition r
      *z = *r;

      if( n==1 ) {
        *p = *z;
      } else {
        beta = 1./( d - beta/alpha );
        p->add( 1., *z, beta, *p );
      }
      x->add( 1., *x, alpha, *p );
      op_->computeResidual( b, *x, *r );
      if( !out_.is_null() )
        *out_ << r->norm() << "\n";
    }

    x0 = *x;
  }



  void applyPIMP( const DomainFieldT& b, RangeFieldT& x0 ) const {
    // after Gutknecht: three-term Chebyshev iteration
    //
    Teuchos::RCP<RangeFieldT> r = b.clone( ECopy::Shallow );

    Teuchos::RCP<DomainFieldT> p = x0.clone( ECopy::Shallow ); // ~ init(0)
    Teuchos::RCP<DomainFieldT> x = x0.clone( ECopy::Deep ); // x->assign( x0 )


    Scalar alpha = ( lamMax_ + lamMin_ )/2;
    Scalar c = std::abs( lamMax_ - lamMin_ )/2;
    Scalar eta = -alpha/c;

    // r-1 = o
    //rp->init( 0. );
    //p->init( 0. );

    // r0 = B-Ax_0
    //x->assign( x );
    op_->computeResidual( b, *x, *r );

    if( !out_.is_null() )
      *out_ << r->norm() << "\n";

    Scalar beta = 0.;
    Scalar gamma = -alpha;

    for( int n=0; n<numIters_; ++n ) {

      if( n==1 )
        beta = -c*c/alpha/2.;

      if( n>= 2 )
        beta = c*c/4./gamma;

      if( n>=1 )
        gamma = -( alpha + beta );

      //x[n+1] = -( r + alpha*x[n] + beta*x[n-1] )/gamma;
      p->add( -beta/gamma, *p, -alpha/gamma, *x );
      p->add( 1., *p, -1./gamma, *r );

      //x->level();
      //r[n+1] = ( A*r[n] - alpha*r[n] - beta*r[n-1] )/gamma;
      // three-ter recursion, explicitly computed residuals
      op_->computeResidual( b, *p, *r );

      p.swap( x );

      if( !out_.is_null() )
        *out_ << r->norm() << "\n";
    }
    x0 = *x;
  }

public:

  void assignField( const DomainFieldT& mv ) {
    //op_->assignField( mv );
  };

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return op_->space();
  };

  void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
    op_->setParameter( para );
  }

  bool hasApplyTranspose() const {
    return op_->hasApplyTranspose();
  }

  const std::string getLabel() const {
    return "Chebyshev";
  };

  void print( std::ostream& out=std::cout ) const {
    out << getLabel() << ":\n";
    //op_->print( out );
  }


}; // end of class Chebyshev


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CHEBYSHEV_HPP
