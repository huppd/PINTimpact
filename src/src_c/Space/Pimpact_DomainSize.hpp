#pragma once
#ifndef PIMPACT_DOMAINSIZE_HPP
#define PIMPACT_DOMAINSIZE_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"



//extern "C" {
//void fsetDS( const double& L1, const double& L2, const double& L3 );
//}



namespace Pimpact{



/// \brief Domain or physical set up would be better names
/// \ingroup domain
template<class Scalar>
class DomainSize {

  template<class ST>
  friend Teuchos::RCP<const DomainSize<ST> > createDomainSize();

  template<class ST>
  friend Teuchos::RCP<const DomainSize<ST> > createDomainSize( int dim, ST L1, ST L2, ST L3 );

  template<class ST>
  friend Teuchos::RCP<const DomainSize<ST> > createDomainSize( int dim, ST re, ST alpha2, ST L1, ST L2, ST L3 );

public:

  typedef const Teuchos::Tuple<Scalar,3> TS3;

protected:

  const int dim_;

  const Scalar re_;

  const Scalar alpha2_;

  TS3 domainSize_;

  TS3 origin_;


  DomainSize( int dim, Scalar L1, Scalar L2, Scalar L3 ):
    dim_(dim),re_(1.),alpha2_(1.),
    domainSize_( Teuchos::tuple(L1, L2, L3) ),
    origin_( Teuchos::tuple( 0.,0.,0.) ) {};

  DomainSize( int dim, Scalar re, Scalar alpha2, Scalar L1, Scalar L2, Scalar L3 ):
    dim_(dim),re_(re),alpha2_(alpha2),
    domainSize_( Teuchos::tuple(L1, L2, L3) ),
    origin_( Teuchos::tuple( 0.,0.,0.) ) {};

  DomainSize( int dim, TS3 domainSize ):
    dim_(dim),re_(1.),alpha2_(1.),
    domainSize_( domainSize ),
    origin_( Teuchos::tuple( 0.,0.,0.) ) {};

public:

  const int& getDim() const { return( dim_ ); }

  const Scalar* getSize() const { return( domainSize_.getRawPtr() ); }

  const Scalar& getSize( int i) const { return( domainSize_[i] ); }
  const Scalar& getSize( ECoord i) const { return( domainSize_[ (int)i ] ); }

  const Scalar* getOrigin() const { return( origin_.getRawPtr() ); }

  const Scalar& getOrigin( int i) const { return( origin_[i] ); }

  const Scalar& getRe() const { return( re_ ); }

  const Scalar& getAlpha2() const { return( alpha2_ ); }

//  void set_Impact() const {
//    fsetDS( domainSize_[0], domainSize_[1], domainSize_[2] );
//  };

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
template<class Scalar=double>
Teuchos::RCP<const DomainSize<Scalar> > createDomainSize( int dim, Scalar L1, Scalar L2, Scalar L3 ) {
  return(
      Teuchos::rcp(
          new DomainSize<Scalar>( dim, 1., 1., L1, L2, L3 ) ) );
}



/// \relates DomainSize
template<class Scalar>
Teuchos::RCP<const DomainSize<Scalar> > createDomainSize( int dim, Scalar re, Scalar alpha2, Scalar L1, Scalar L2, Scalar L3 ) {
  return(
      Teuchos::rcp(
          new DomainSize<Scalar>( dim, re, alpha2, L1, L2, L3 ) ) );
}



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_DOMAINSIZE_HPP