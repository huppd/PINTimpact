#pragma once
#ifndef PIMPACT_RESTRICTIONTIMEOP_HPP
#define PIMPACT_RESTRICTIONTIMEOP_HPP

#include "Teuchos_RCP.hpp"
#include <BelosTypes.hpp>
#include "Pimpact_Types.hpp"
#include "Pimpact_TimeField.hpp"



namespace Pimpact {
template<class OperatorT>
class RestrictionTimeOp  {

public:

  typedef TimeField<typename OperatorT::DomainFieldT> DomainFieldT;
  typedef TimeField<typename OperatorT::RangeFieldT> RangeFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

  typedef typename SpaceT::Ordinal Ordinal;

protected:

  Teuchos::RCP<OperatorT> op_; // restriction in space
  Teuchos::RCP<typename OperatorT::RangeFieldT> temp_; 

public:
    
    RestrictionTimeOp( const Teuchos::RCP<const SpaceT>& spaceF,
                       const Teuchos::RCP<const SpaceT>& spaceC):
    		        op_(Teuchos::rcp(new OperatorT(spaceF,spaceC))),
			temp_( create<typename OperatorT::RangeFieldT>(spaceC) ){};


// x is fn in this case
void apply( const DomainFieldT& x, RangeFieldT& y) const {
      
      Ordinal d = spaceF()->nLoc(3)/spaceC()->nLoc(3);

      x.exchange();

      temp_->init(0.);

	std::cout << "CR space" << spaceC()->nLoc(3) << std::endl;
	
	for( int i=0; i < spaceF()->nLoc(3); ++i ) {
	
		if ( i%d>0 ){
		//	op_->apply( x.getConstField(i), *temp_ );
                  //   	y.getFieldPtr((i-1)/d)->add(1.,y.getField((i-1)/d),0.25,*temp_);
			
		//	if (i == spaceF()->nLoc(3)-1)
		//		 y.getFieldPtr(0)->add(0.25,*temp_,1.,y.getField(0));
		}
		else {

			std::cout << "i/d" << i/d << std::endl;
			op_->apply( x.getConstField(i), y.getField(i/d) );	
		//	y.getFieldPtr(i/d)->add(0.25,*temp_,0.5,y.getField(i/d));
		}		
      }

      //if (d > 1) {
          
          //}
      //}
   y.changed();
 }

    Teuchos::RCP<const SpaceT> spaceC() const { return(op_->get_spaceC_()); };
    Teuchos::RCP<const SpaceT> spaceF() const { return(op_->get_spaceF_()); };

 Teuchos::RCP<OperatorT> getOperatorPtr() { return( op_ ); }
	
};

template< class OpT >
Teuchos::RCP< RestrictionTimeOp<OpT> > createRestrictionTimeOp(
    const Teuchos::RCP<const typename OpT::FSpaceT>& spaceF, const Teuchos::RCP<const typename OpT::CSpaceT>& spaceC ) {
	return( Teuchos::rcp( new RestrictionTimeOp<OpT>( spaceF,spaceC ) ) );
	}

} // end of namespace Pimpact


#endif


