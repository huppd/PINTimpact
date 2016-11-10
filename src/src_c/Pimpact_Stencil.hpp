#pragma once
#ifndef PIMPACT_STENCIL_HPP
#define PIMPACT_STENCIL_HPP




namespace Pimpact {


template<class Scalar, class Ordinal>
class Stencil {

protected:

	using ScalarArray = Scalar*;

	ScalarArray c_;

	const Ordinal ss_; /// < \todo make template parameter
	const Ordinal nn_;
	const int bl_;	/// < \todo make template parameter
	const int bu_;	/// < \todo make template parameter
	const int w_;

public:

	Stencil( const Ordinal& ss, const Ordinal& nn, const int& bl, const int& bu ) : 
		ss_(ss), nn_(nn), bl_(bl), bu_(bu), w_(bu_-bl_+1) {

			Ordinal nTemp = ( nn_ - ss_ + 1 )*w_;

			c_ = new Scalar[ nTemp ];

			for( Ordinal i=0; i<nTemp; ++i )
				c_[i] = 0.;
		}

	~Stencil() {
		delete[] c_;
	}

	ScalarArray get() {
		return( c_ );
	}

	constexpr Scalar& operator()( const Ordinal& index, const int& offset ) {
		return( c_[ offset-bl_ + (index-ss_)*w_ ] );
	};

	constexpr const Scalar& at( const Ordinal& index, const int& offset ) const {
		return( c_[ offset-bl_ + (index-ss_)*w_ ] );
	}

	constexpr const int& bl() const {
		return( bl_ );
	}

	constexpr const int& bu() const {
		return( bu_ );
	}

	void print( std::ostream& out=std::cout ) const {
		out << std::setw(8) << "bl: ";
		for( int k=bl_; k<=bu_; ++k ) 
			out << std::setw(9) << k ;
		out << " :bu\n";

		for( int bla=0; bla<9*w_+5; bla++ )
			out << "-";
		out << "\n";

		for( int i=ss_; i<=nn_; ++i ) {
			out << "i: " << std::setw(3) << i << " (";
			for( int k=bl_; k<=bu_; ++k ) 
				out << std::setw(9) << at(i,k);
			out << ")\n";
		}
	}

}; // end of class Stencil 


} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_STENCIL_HPP
