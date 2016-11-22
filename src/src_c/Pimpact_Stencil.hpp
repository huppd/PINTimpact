#pragma once
#ifndef PIMPACT_STENCIL_HPP
#define PIMPACT_STENCIL_HPP




namespace Pimpact {


/// \brief 
///
/// \tparam Scalar
/// \tparam Ordinal
/// \tparam ss
/// \tparam bl_
/// \tparam bu_
template<class Scalar, class Ordinal, int ss, int bl_, int bu_>
class Stencil {

protected:

	using ScalarArray = Scalar*;

	ScalarArray c_;

	const Ordinal nn_;
	//const int bl_;	/// < \todo make template parameter
	//const int bu_;	/// < \todo make template parameter
	static const int w_ = bu_-bl_+1;

public:

	Stencil( const Ordinal& nn ) : 
		nn_(nn) {

			assert( bl_<=bu_ );
			assert( ss<=nn_ );
			
			Ordinal nTemp = ( nn_ - ss + 1 )*w_;

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
		assert( offset>=bl_ );
		assert( offset<=bu_ );
		assert( index>=ss );
		assert( index<=nn_ );
		return( c_[ offset-bl_ + (index-ss)*w_ ] );
	};


	constexpr const Scalar& at( const Ordinal& index, const int& offset ) const {
		assert( offset>=bl_ );
		assert( offset<=bu_ );
		assert( index>=ss );
		assert( index<=nn_ );
		return(
				c_[ offset-bl_ + (index-ss)*w_ ] );
	}

	static inline constexpr int bl() {
		return( bl_ );
	}

	static inline constexpr int bu() {
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

		for( int i=ss; i<=nn_; ++i ) {
			out << "i: " << std::setw(3) << i << " (";
			for( int ii=bl_; ii<=bu_; ++ii ) 
				out << std::setw(9) << at(i,ii);
			out << ")\n";
		}
	}

}; // end of class Stencil 


} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_STENCIL_HPP
