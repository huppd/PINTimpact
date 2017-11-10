#pragma once
#ifndef PIMPACT_DOMAINSIZE_HPP
#define PIMPACT_DOMAINSIZE_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \brief Domain or physical set up would be better names
/// \ingroup SpaceObject
template<class ScalarT, int sd>
class DomainSize {

  template<class ST, int sd_>
  friend Teuchos::RCP<const DomainSize<ST,sd_> > createDomainSize(
    ST re, ST alpha2,
    ST L1, ST L2, ST L3,
    ST x1, ST x2, ST x3 );

public:

  using TS3 = const Teuchos::Tuple<ScalarT,3>;

protected:

  const ScalarT re_;

  const ScalarT alpha2_;

  TS3 domainSize_;

  TS3 origin_;

  DomainSize(
    ScalarT re, ScalarT alpha2,
    ScalarT L1, ScalarT L2, ScalarT L3,
    ScalarT x1, ScalarT x2, ScalarT x3 ):
    re_(re),alpha2_(alpha2),
    domainSize_( Teuchos::tuple(L1, L2, L3) ),
    origin_( Teuchos::tuple( x1, x2, x3 ) ) {

    static_assert( sd!=2 || sd!=3, "spatial dimension not valid" );
  };


public:

  /// \name getter
  /// @{

  constexpr const ScalarT* getSize() const {
    return domainSize_.getRawPtr();
  }

  constexpr const ScalarT& getSize( const int i) const {
    return domainSize_[i];
  }

  constexpr const ScalarT* getOrigin() const {
    return origin_.getRawPtr();
  }

  constexpr const ScalarT& getOrigin( const int i) const {
    return origin_[i];
  }

  constexpr const ScalarT& getRe() const {
    return re_;
  }

  constexpr const ScalarT& getAlpha2() const {
    return alpha2_;
  }

  ///  @}

  void print( std::ostream& out=std::cout ) const {
    out << "\tspatial dim: " << sd << "\n"
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
template<class ScalarT, int sd>
Teuchos::RCP<const DomainSize<ScalarT,sd> >
createDomainSize(
  ScalarT re, ScalarT alpha2,
  ScalarT L1, ScalarT L2, ScalarT L3,
  ScalarT x1, ScalarT x2, ScalarT x3 ) {

  return Teuchos::rcp( new DomainSize<ScalarT,sd>( re, alpha2, L1, L2, L3, x1, x2, x3 ) );
}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_DOMAINSIZE_HPP
