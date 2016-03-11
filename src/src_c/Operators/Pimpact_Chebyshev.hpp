#pragma once
#ifndef PIMPACT_CHEBYSHEV_HPP
#define PIMPACT_CHEBYSHEV_HPP


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Space.hpp" // just for create<>




namespace Pimpact {



template<class OperatorT>
class Chebyshev {

public:

	using DomainFieldT = typename OperatorT::DomainFieldT;
	using RangeFieldT = typename OperatorT::RangeFieldT;

	using SpaceT = typename DomainFieldT::SpaceT;

	using Scalar = typename SpaceT::Scalar;

protected:

	int numIters_;

	bool zeroStartingSolution_ = false;

	Scalar lamMax_;
	Scalar lamMin_;
	Scalar eigRatio_ = 1./30.;

	Teuchos::RCP<RangeFieldT> r_;

	mutable Teuchos::RCP<DomainFieldT> xp_;
	mutable Teuchos::RCP<DomainFieldT> x_;

	Teuchos::RCP<const OperatorT> op_;

	mutable Teuchos::RCP<DomainFieldT> D_;
	mutable Teuchos::RCP<DomainFieldT> D_inv_;

	Teuchos::RCP<std::ostream> out_;

public:

	Chebyshev(
			const Teuchos::RCP<const OperatorT>& op,
			const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList() ):
    numIters_( pl->get<int>( "numIters", 8 ) ),
		lamMax_( pl->get<Scalar>( "max EV", 0. ) ),
		lamMin_( pl->get<Scalar>( "min EV", 0. ) ),
    r_(  Teuchos::rcp( new RangeFieldT( op->space()  ) ) ),
    xp_( Teuchos::rcp( new DomainFieldT( op->space() ) ) ),
    x_(  Teuchos::rcp( new DomainFieldT( op->space() ) ) ),
    op_(op) {

			D_ = r_->clone();
			D_inv_ = r_->clone();
			D_->init( 1. );
			D_inv_->init( 1. );

			if( pl->get<bool>( "with output", false ) )
				out_ = Pimpact::createOstream( "conv_Cheb.txt" );

			if( std::abs( lamMax_-lamMin_ )<1.e-32 ) {
				x_->random();

				Scalar lamp = 0.;
				for( int i=0; i<500; ++i ) {
					op_->apply( *x_, *xp_ );
					Scalar lam = xp_->dot( *x_ )/x_->dot( *x_ );
					//std::cout << "lambda: " << lam << "\t" << std::abs( lamp-lam )/std::abs(lam) << "\n";
					xp_->scale( 1./xp_->norm() );
					x_.swap( xp_ );
					if( std::abs( lamp-lam )/std::abs(lam) < 1.e-3 )
						break;
					else
						lamp=lam;
				}
				//x_->write( 999 );

				op_->apply( *x_, *xp_ );
				lamMax_ = xp_->dot( *x_ )/x_->dot( *x_ );
				std::cout << "lamMax: " << lamMax_ << "\n";

				lamMax_ *= 1.1;
				lamMin_ = lamMax_*eigRatio_;

			}

		}

	void apply( const DomainFieldT& b, RangeFieldT& x ) const { 
		//applyPIMP( b, x );
		applyIFFPACK( b, x );
	}

protected:

	void applyIFFPACK( const DomainFieldT& b, RangeFieldT& x ) const { 

#ifdef HAVE_TEUCHOS_DEBUG
		using std::cerr;
		using std::endl;
		cerr << "\\|B\\|_{\\infty} = " << b.norm( Belos::InfNorm ) << endl;
		cerr << "\\|X\\|_{\\infty} = " << x.norm( Belos::InfNorm ) << endl;
#endif // HAVE_TEUCHOS_DEBUG

		if( numIters_<=0 ) {
			return;
		}

		const Scalar zero = Teuchos::as<Scalar>(0);
		const Scalar one = Teuchos::as<Scalar>(1);
		const Scalar two = Teuchos::as<Scalar>(2);

		// Initialize coefficients
		const Scalar alpha = lamMax_ / eigRatio_;
		const Scalar beta = Teuchos::as<Scalar> (1.1) * lamMax_;
		const Scalar delta = two / (beta - alpha);
		const Scalar theta = (beta + alpha) / two;
		const Scalar s1 = theta * delta;

#ifdef HAVE_TEUCHOS_DEBUG
#ifdef IFPACK_DETAILS_CHEBYSHEV_DEBUG
		cerr << "alpha = " << alpha << endl
			<< "beta = " << beta << endl
			<< "delta = " << delta << endl
			<< "theta = " << theta << endl
			<< "s1 = " << s1 << endl;
#endif // IFPACK_DETAILS_CHEBYSHEV_DEBUG
#endif // HAVE_TEUCHOS_DEBUG

		// Fetch cached temporary vectors.
		Teuchos::RCP<DomainFieldT> V_ptr = b.clone( Pimpact::ShallowCopy );
		Teuchos::RCP<DomainFieldT> W_ptr = b.clone( Pimpact::ShallowCopy );

		// mfh 28 Jan 2013: We write V1 instead of V, so as not to confuse
		// the multivector V with the typedef V (for Tpetra::Vector).
		//MV V1 (B.getMap (), B.getNumVectors (), false);
		//MV W (B.getMap (), B.getNumVectors (), false);
		DomainFieldT& V1 = *V_ptr;
		DomainFieldT& W = *W_ptr;

#ifdef HAVE_TEUCHOS_DEBUG
		cerr << "Iteration " << 1 << ":" << endl
			<< "- \\|D\\|_{\\infty} = " << D_->norm( Belos::InfNorm ) << endl;
#endif // HAVE_TEUCHOS_DEBUG

		// Special case for the first iteration.
		if( !zeroStartingSolution_ ) {
			op_->computeResidual( b, x, V1 ); // V1 = B - A*X

#ifdef HAVE_TEUCHOS_DEBUG
			cerr << "- \\|B - A*X\\|_{\\infty} = " << V1.norm( Belos::InfNorm ) << endl;
#endif // HAVE_TEUCHOS_DEBUG

			//solve (W, one/theta, D_inv, V1); // W = (1/theta)*D_inv*(B-A*X)
			W.reciprocal( *D_ );
			W.scale( V1 );
			W.scale( one/theta );

#ifdef HAVE_TEUCHOS_DEBUG
			cerr << "- \\|W\\|_{\\infty} = " << W.norm( Belos::InfNorm ) << endl;
#endif // HAVE_TEUCHOS_DEBUG

			x.add( one, x, one, W ); // X = X + W
		} else {
			//solve (W, one/theta, D_inv, B); // W = (1/theta)*D_inv*B
			W.reciprocal( *D_ );
			W.scale( b );
			W.scale( one/theta );

#ifdef HAVE_TEUCHOS_DEBUG
			cerr << "- \\|W\\|_{\\infty} = " << W.norm( Belos::InfNorm ) << endl;
#endif // HAVE_TEUCHOS_DEBUG

			//Tpetra::deep_copy(X, W); // X = 0 + W
			x.assign( W );
		}
#ifdef HAVE_TEUCHOS_DEBUG
		cerr << "- \\|X\\|_{\\infty} = " << x.norm( Belos::InfNorm ) << endl;
#endif // HAVE_TEUCHOS_DEBUG

		// The rest of the iterations.
		Scalar rhok = one / s1;
		Scalar rhokp1, dtemp1, dtemp2;
		for( int deg=1; deg<numIters_; ++deg) {

#ifdef HAVE_TEUCHOS_DEBUG
			cerr << "Iteration " << deg+1 << ":" << endl;
			cerr << "- \\|D\\|_{\\infty} = " << D_->norm( Belos::InfNorm ) << endl;
			cerr << "- \\|B\\|_{\\infty} = " << b.norm( Belos::InfNorm ) << endl;
			//cerr << "- \\|A\\|_{\\text{frob}} = " << A_->getFrobeniusNorm () << endl;
			cerr << "- rhok = " << rhok << endl;
			V1.init( zero );
#endif // HAVE_TEUCHOS_DEBUG

			//computeResidual (V1, B, A, X); // V1 = B - A*X
			op_->computeResidual( b, x, V1 ); // V1 = B - A*X

#ifdef HAVE_TEUCHOS_DEBUG
			cerr << "- \\|B - A*X\\|_{\\infty} = " << V1.norm( Belos::InfNorm ) << endl;
#endif // HAVE_TEUCHOS_DEBUG

			rhokp1 = one / (two * s1 - rhok);
			dtemp1 = rhokp1 * rhok;
			dtemp2 = two * rhokp1 * delta;
			rhok = rhokp1;

#ifdef HAVE_TEUCHOS_DEBUG
			cerr << "- dtemp1 = " << dtemp1 << endl
				<< "- dtemp2 = " << dtemp2 << endl;
#endif // HAVE_TEUCHOS_DEBUG

			//W.scale( dtemp1 );
			//W.elementWiseMultiply( dtemp2, D_inv, V1, one ); // W = one*W + dtemp2*D_inv.*V1
			{
				auto temp = W.clone( Pimpact::ShallowCopy );
				temp->reciprocal( *D_ );
				temp->scale( V1 );
				W.add( dtemp1, W, dtemp2, *temp );
			}

			//X.update (one, W, one);
			x.add( one, x, one, W ); // X = X + W

#ifdef HAVE_TEUCHOS_DEBUG
			cerr << "- \\|W\\|_{\\infty} = " << W.norm( Belos::InfNorm ) << endl;
			cerr << "- \\|X\\|_{\\infty} = " << x.norm( Belos::InfNorm ) << endl;
#endif // HAVE_TEUCHOS_DEBUG

		}
	}

	void applyPIMP( const DomainFieldT& b, RangeFieldT& x ) const {
		// after Gutknecht: three-term Chebyshev iteration

		Scalar alpha = ( lamMax_ + lamMin_ )/2;
		Scalar c = std::abs( lamMax_ - lamMin_ )/2;
		Scalar eta = -alpha/c;

		// r_-1 = o
		//rp->init( 0. );
		xp_->init( 0. );

		// r_0 = B-Ax_0
		x_->assign( x );
		op_->computeResidual( b, *x_, *r_ );


		if( !out_.is_null() )
			*out_ << r_->norm() << "\n";
		
		Scalar beta = 0.;
		Scalar gamma = -alpha;

		for( int n=0; n<numIters_; ++n ) {

			if( n==1 )
				beta = -c*c/alpha/2.;

			if( n>= 2 )
				beta = c*c/4./gamma;

			if( n>=1 )
				gamma = -( alpha + beta );

			//x[n+1] = -( r_ + alpha*x[n] + beta*x[n-1] )/gamma;
			xp_->add( -beta/gamma, *xp_, -alpha/gamma, *x_ );
			xp_->add( 1., *xp_, -1./gamma, *r_ );

			//x_->level();
			//r[n+1] = ( A*r[n] - alpha*r[n] - beta*r[n-1] )/gamma;
			// three-ter recursion, explicitly computed residuals
			op_->computeResidual( b, *xp_, *r_ );

			xp_.swap( x_ );

			if( !out_.is_null() )
				*out_ << r_->norm() << "\n";
		}
		x.assign( *x_ );
	}

public:

	void assignField( const DomainFieldT& mv ) {
		op_->assignField( mv );
	};

	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		op_->setParameter( para );
	}

	bool hasApplyTranspose() const { return( op_->hasApplyTranspose() ); }

	const std::string getLabel() const { return( "Chebyshev( " + op_->getLabel() + " )" ); };

	void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
		//op_->print( out );
	}


}; // end of class Chebyshev


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CHEBYSHEV_HPP
