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

	Scalar lamMax_;
	Scalar lamMin_;

	Teuchos::RCP<RangeFieldT> r_;

	mutable Teuchos::RCP<DomainFieldT> xp_;
	mutable Teuchos::RCP<DomainFieldT> x_;

	Teuchos::RCP<const OperatorT> op_;

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
				lamMin_ = 0.3*lamMax_;

			}

		}


	void apply( const DomainFieldT& b, RangeFieldT& x ) const {
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
