#pragma once
#ifndef PIMPACT_STENCIL_HPP
#define PIMPACT_STENCIL_HPP


#include <iostream>

#include "Teuchos_ScalarTraits.hpp"


namespace Pimpact {


/// \brief 
///
/// \tparam Scalar
/// \tparam Ordinal
/// \tparam ss start index
/// \tparam lb lower stencil bound
/// \tparam ub_ upper stencil bound
template<class Scalar, class Ordinal, int ss, int lb, int ub>
class Stencil {

protected:

	using ScalarArray = Scalar*;

	ScalarArray c_;

	Ordinal nn_;
	static const int w_ = ub-lb+1;

public:

	Stencil() : c_(nullptr),nn_(ss) {}

	Stencil( const Ordinal& nn ) : nn_(nn) {

		static_assert( lb<=ub, "Stencil width cannot be negative" );
		assert( ss<=nn_ );

		Ordinal nTemp = ( nn_ - ss + 1 )*w_;

		c_ = new Scalar[ nTemp ];

		std::fill_n( c_, nTemp, Teuchos::ScalarTraits<Scalar>::zero() );
	}

	Stencil( const Stencil& that ) : Stencil(that.nn_) { 

		Ordinal nTemp = ( nn_ - ss + 1 )*w_;

		std::copy_n( that.c_, nTemp, c_ );
	}

	friend void swap( Stencil& one, Stencil& two ) {
		using std::swap;
		swap( one.nn_, two.nn_ );
		swap( one.c_, two.c_ );
	}

	Stencil( Stencil&& that ) : Stencil() { 
		swap( *this, that );
	}

	Stencil& operator=( Stencil that ) {
		swap( *this, that );
		return *this;
	}

	~Stencil() { delete[] c_; }

	inline constexpr ScalarArray get() const { return( c_ ); }

	inline constexpr Scalar& operator()( const Ordinal& index, const int& offset ) {
		assert( offset>=lb );
		assert( offset<=ub );
		assert( index>=ss );
		assert( index<=nn_ );
		return( c_[ offset-lb + (index-ss)*w_ ] );
	};


	static inline constexpr int bl() {
		return( lb );
	}

	static inline constexpr int bu() {
		return( ub );
	}

	void print( std::ostream& out=std::cout ) const {
		out << std::setw(8) << "bl: ";
		out << std::scientific;
		out << std::setprecision( 3 );
		for( int k=lb; k<=ub; ++k ) 
			out << std::setw(12) << k ;
		out << " :bu\n";

		for( int bla=0; bla<12*w_+5; bla++ )
			out << "-";
		out << "\n";

		for( int i=ss; i<=nn_; ++i ) {
			out << "i: " << std::setw(3) << i << " (";
			for( int ii=lb; ii<=ub; ++ii ) 
				out << std::setw(12) << (*this)(i,ii);
			out << ")\n";
		}
	}

}; // end of class Stencil 


} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_STENCIL_HPP
