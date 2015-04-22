#pragma once
#ifndef PIMPACT_MGTRANSFERS_HPP
#define PIMPACT_MGTRANSFERS_HPP


#include "Pimpact_MGSpaces.hpp"
#include "Pimpact_TransferOp.hpp"
#include "Pimpact_RestrictionOp.hpp"
#include "Pimpact_InterpolationOp.hpp"




namespace Pimpact {



/// \ingroup MG
/// \todo make Interpolation/Restriction templated
template<class MGSpacesT,
  template<class,class> class TransT,
  template<class> class RestrT,
  template<class> class InterT >
class MGTransfers {

public:

  typedef typename MGSpacesT::FSpaceT FSpaceT;
  typedef typename MGSpacesT::CSpaceT CSpaceT;

  typedef typename FSpaceT::Scalar  Scalar;
  typedef typename FSpaceT::Ordinal Ordinal;

  static const int dimension = FSpaceT::dimension;

  static const int dimNCF = FSpaceT::dimNC;
  static const int dimNCC = CSpaceT::dimNC;


  typedef TransT<FSpaceT,CSpaceT> TransferOpT;
  typedef RestrT<CSpaceT> RestrictionOpT;
  typedef InterT<CSpaceT> InterpolationOpT;

//  template<
//    class MGSpacesTT,
//    template<class,class> class TransTT,
//    template<class> class RestrTT,
//    template<class> class InterTT >
//  friend
//  Teuchos::RCP<const MGTransfers<MGSpacesTT,TransTT,RestrTT,InterTT> >
//  createMGTransfers(
//      const Teuchos::RCP<const MGSpacesTT>& space );

protected:

  Teuchos::RCP<const MGSpacesT> mgSpaces_;

  Teuchos::RCP<const TransferOpT> transferOp_;

  std::vector< Teuchos::RCP<const RestrictionOpT> >   restrictionOps_;
  std::vector< Teuchos::RCP<const InterpolationOpT> > interpolationOps_;

public:

  MGTransfers(
      const Teuchos::RCP<const MGSpacesT>& mgSpaces ):
        mgSpaces_(mgSpaces),
        transferOp_(),
        restrictionOps_(),
        interpolationOps_() {

    transferOp_ = create<TransferOpT>( mgSpaces_->get(), mgSpaces_->get(0) );

		for( unsigned i=0; i < mgSpaces_->getNGrids()-1; ++i ) {
			if( mgSpaces_->participating(i) ) {
				restrictionOps_.push_back(
						Teuchos::rcp(
							new RestrT<CSpaceT>(
								mgSpaces_->get(i),
								mgSpaces_->get(i+1),
								mgSpaces_->get()->getProcGridSize()->getTuple()
								)
							)
						);
				interpolationOps_.push_back( create<InterT>( mgSpaces_->get(i+1), mgSpaces_->get(i) ) );
			}
		}

	// not working on brutus
    //interpolationOps_.shrink_to_fit();

  }

public:

  Teuchos::RCP<const TransferOpT>       getTransferOp     (       ) const { return( transferOp_ ); }

  /// \brief gets ith RestrictionOp, similar to python i=-1 is gets you the coarses space
  Teuchos::RCP<const RestrictionOpT>    getRestrictionOp  ( int i ) const {
    if( i<0 )
      return( restrictionOps_[mgSpaces_->getNGrids()+i] );
    else
      return( restrictionOps_[i] );
  }

  /// \brief gets ith InterpolationOp, similar to python i=-1 is gets you the coarses space
  Teuchos::RCP<const InterpolationOpT>  getInterpolationOp( int i ) const {
    if( i<0 )
      return( interpolationOps_[mgSpaces_->getNGrids()+i] );
    else
      return( interpolationOps_[i] );
  }

  void print(  std::ostream& out=std::cout ) const {

    transferOp_->print(out);

    for( int i = 0; i<restrictionOps_.size(); ++i ) {
			if( mgSpaces_->participating(i) ) {
					out << "-------- restrictor: "<< i << "--------\n";
					restrictionOps_[i]->print(out);
			}
    }
    for( int i = 0; i<interpolationOps_.size(); ++i ) {
			if( mgSpaces_->participating(i) ) {
				out << "-------- interpolator: "<< i << "--------\n";
				interpolationOps_[i]->print(out);
			}
    }
  }

}; // end of class MGTransfers



/// \relates MGTransfers
template<
  template<class,class> class TransT,
  template<class> class RestrT,
  template<class> class InterT,
  class MGSpacesT >
Teuchos::RCP<const MGTransfers<MGSpacesT,TransT,RestrT,InterT> >
createMGTransfers(
    const Teuchos::RCP<const MGSpacesT>& mgSpaces ) {

  return(
      Teuchos::rcp(
          new MGTransfers<MGSpacesT,TransT,RestrT,InterT>( mgSpaces )
      )
  );

}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MGTRANSFERS_HPP
