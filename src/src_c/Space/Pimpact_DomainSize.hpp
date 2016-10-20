#pragma once
#ifndef PIMPACT_DOMAINSIZE_HPP
#define PIMPACT_DOMAINSIZE_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"




namespace Pimpact{



/// \brief Domain or physical set up would be better names
/// \ingroup SpaceObject
template<class ScalarT>
class DomainSize {

	template<class ST>
	friend Teuchos::RCP<const DomainSize<ST> > createDomainSize(
			int dim, ST re, ST alpha2,
			ST L1, ST L2, ST L3,
			ST x1, ST x2, ST x3 );

public:

	using TS3 = const Teuchos::Tuple<ScalarT,3>;

protected:

	const int dim_;

	const ScalarT re_;

	const ScalarT alpha2_;

	TS3 domainSize_;

	TS3 origin_;

	DomainSize(
			int dim, ScalarT re, ScalarT alpha2,
			ScalarT L1, ScalarT L2, ScalarT L3,
			ScalarT x1, ScalarT x2, ScalarT x3 ):
		dim_(dim),re_(re),alpha2_(alpha2),
		domainSize_( Teuchos::tuple(L1, L2, L3) ),
		origin_( Teuchos::tuple( x1, x2, x3 ) ) {};


public:

	/// \name getter
	/// @{ 

	/// \todo make dim_ template parameter(e.g. spaceDim) here or somewhere else?
	constexpr const int& getDim() const { return( dim_ ); }

	constexpr const ScalarT* getSize() const { return( domainSize_.getRawPtr() ); }

	constexpr const ScalarT& getSize( const int& i) const { return( domainSize_[i] ); }

	constexpr const ScalarT* getOrigin() const { return( origin_.getRawPtr() ); }

	constexpr const ScalarT& getOrigin( const int& i) const { return( origin_[i] ); }

	constexpr const ScalarT& getRe() const { return( re_ ); }

	constexpr const ScalarT& getAlpha2() const { return( alpha2_ ); }

	///  @} 

	void print( std::ostream& out=std::cout ) const {
		out << "\tspatial dim: " << dim_ << "\n"
			<< "\tRe= "      << re_ << "\n"
			<< "\talpha^2= " << alpha2_ << "\n"
			<< "\tlx= "      << domainSize_[0]
			<< "\tly= "      << domainSize_[1]
			<< "\tlz= "      << domainSize_[2] << "\n"
			<< "\tox= "      << origin_[0]
			<< "\toy= "      << origin_[1]
			<< "\toz= "      << origin_[2] << "\n";
	};

}; // end of DomainSize




/// \relates DomainSize
template<class ScalarT>
Teuchos::RCP<const DomainSize<ScalarT> >
createDomainSize(
		int dim, ScalarT re, ScalarT alpha2,
		ScalarT L1, ScalarT L2, ScalarT L3,
		ScalarT x1, ScalarT x2, ScalarT x3 ) {

	return( Teuchos::rcp(
          new DomainSize<ScalarT>( dim, re, alpha2, L1, L2, L3, x1, x2, x3 ) ) );
}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_DOMAINSIZE_HPP
