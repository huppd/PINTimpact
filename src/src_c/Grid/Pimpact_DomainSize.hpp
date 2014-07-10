#pragma once
#ifndef PIMPACT_DOMAINSIZE_HPP
#define PIMPACT_DOMAINSIZE_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"



extern "C" {
void fsetDS(const double& L1, const double& L2, const double& L3 );
}



namespace Pimpact{



/// \brief Domain or physical set up would be better names
template< class Scalar=double >
class DomainSize {

public:

  typedef const Teuchos::Tuple<Scalar,3> TS3;

protected:

  Scalar re_;

  Scalar alpha2_;

  TS3 domainSize_;

public:

  DomainSize( Scalar L1, Scalar L2, Scalar L3 ):
    re_(1.),alpha2_(1.),
    domainSize_( Teuchos::tuple(L1, L2, L3) ) {
    set_Impact();
  };

  DomainSize( Scalar re, Scalar alpha2, Scalar L1, Scalar L2, Scalar L3 ):
    re_(re),alpha2_(alpha2),
    domainSize_( Teuchos::tuple(L1, L2, L3) ) {
    set_Impact();
  };

  DomainSize( TS3 domainSize ):
    re_(1.),alpha2_(1.),
    domainSize_( domainSize ) {};

  void set_Impact(){
    fsetDS( domainSize_[0], domainSize_[1], domainSize_[2] );
  };

  void print( std::ostream& out=std::cout ) {
    out << "\tRe= "      << re_ << "\n"
        << "\talpha^2= " << alpha2_ << "\n"
        << "\tlx= "      << domainSize_[0]
        << "\tly= "      << domainSize_[1]
        << "\tlz= "      << domainSize_[2] << "\n";
  };

}; // end of DomainSize



///// \relates DomainSize
//template<class Scalar=double>
//Teuchos::RCP<DomainSize<Scalar> > createDomainSize( Scalar L1=1, Scalar L2=1, Scalar L3=1 ) {
//  return(
//      Teuchos::rcp(
//          new DomainSize<Scalar>( L1, L2, L3 ) ) );
//}

/// \relates DomainSize
template<class Scalar=double>
Teuchos::RCP<DomainSize<Scalar> > createDomainSize( Scalar re=1., Scalar alpha2=1., Scalar L1=1, Scalar L2=1, Scalar L3=1 ) {
  return(
      Teuchos::rcp(
          new DomainSize<Scalar>( re, alpha2, L1, L2, L3 ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_DOMAINSIZE_HPP
