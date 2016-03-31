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
	friend Teuchos::RCP<const DomainSize<ST> > createDomainSize();

	template<class ST>
	friend Teuchos::RCP<const DomainSize<ST> > createDomainSize( int dim, ST L1, ST L2, ST L3 );

	template<class ST>
	friend Teuchos::RCP<const DomainSize<ST> > createDomainSize( int dim, ST re, ST alpha2, ST L1, ST L2, ST L3 );

public:

	using TS3 = const Teuchos::Tuple<ScalarT,3>;

protected:

	const int dim_;

	const ScalarT re_;

	const ScalarT alpha2_;

	TS3 domainSize_;

	TS3 origin_;

	DomainSize( int dim, ScalarT L1, ScalarT L2, ScalarT L3 ):
		dim_(dim),re_(1.),alpha2_(1.),
		domainSize_( Teuchos::tuple(L1, L2, L3) ),
		origin_( Teuchos::tuple( 0.,0.,0.) ) {};

	DomainSize( int dim, ScalarT re, ScalarT alpha2, ScalarT L1, ScalarT L2, ScalarT L3 ):
		dim_(dim),re_(re),alpha2_(alpha2),
		domainSize_( Teuchos::tuple(L1, L2, L3) ),
		origin_( Teuchos::tuple( 0.,0.,0.) ) {};

	DomainSize( int dim, TS3 domainSize ):
		dim_(dim),re_(1.),alpha2_(1.),
		domainSize_( domainSize ),
		origin_( Teuchos::tuple( 0.,0.,0.) ) {};

public:

	/// \name getter
	/// @{ 

	/// \todo make dim_ template parameter(e.g. spaceDim) here or somewhere else?
	inline constexpr const int& getDim() const { return( dim_ ); }

	inline constexpr const ScalarT* getSize() const { return( domainSize_.getRawPtr() ); }

	inline constexpr const ScalarT& getSize( int i) const { return( domainSize_[i] ); }
	inline constexpr const ScalarT& getSize( ECoord i) const { return( domainSize_[ (int)i ] ); }

	inline constexpr const ScalarT* getOrigin() const { return( origin_.getRawPtr() ); }

	inline constexpr const ScalarT& getOrigin( int i) const { return( origin_[i] ); }

	inline constexpr const ScalarT& getRe() const { return( re_ ); }

	inline constexpr const ScalarT& getAlpha2() const { return( alpha2_ ); }

	///  @} 

	void print( std::ostream& out=std::cout ) const {
		out << "\tspatial dim: " << dim_ << "\n"
			<< "\tRe= "      << re_ << "\n"
			<< "\talpha^2= " << alpha2_ << "\n"
			<< "\tlx= "      << domainSize_[0]
			<< "\tly= "      << domainSize_[1]
			<< "\tlz= "      << domainSize_[2] << "\n";
	};

}; // end of DomainSize



/// \relates DomainSize
/// \deprecated
template<class ScalarT>
Teuchos::RCP<const DomainSize<ScalarT> > createDomainSize( int dim, ScalarT L1, ScalarT L2, ScalarT L3 ) {
  return(
      Teuchos::rcp(
          new DomainSize<ScalarT>( dim, 1., 1., L1, L2, L3 ) ) );
}



/// \relates DomainSize
template<class ScalarT>
Teuchos::RCP<const DomainSize<ScalarT> > createDomainSize( int dim, ScalarT re, ScalarT alpha2, ScalarT L1, ScalarT L2, ScalarT L3 ) {
  return(
      Teuchos::rcp(
          new DomainSize<ScalarT>( dim, re, alpha2, L1, L2, L3 ) ) );
}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_DOMAINSIZE_HPP
