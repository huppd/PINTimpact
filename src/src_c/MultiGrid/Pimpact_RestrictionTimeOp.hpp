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
      for( int i=0; i <= spaceF()->nLoc(3); ++i ) {
	
	if (i==0 || i==spaceF()->nLoc(3))
        	
		op_->apply( x.getConstField(i), y.getField(i/d) );
		
	else {
		if ( i%d>0 )
			op_->apply( x.getConstField(i), temp_ );
			
			if (i>1)
                        	y.getFieldPtr(i/d)->add(1.0,y.getField((i-1)/d),0.25,temp_);      
		else	
			op_->apply( x.getConstField(i), y.getField(i/d) );
			y.getFieldPtr(i/d)->add(0.25,temp_,0.5,y.getField(i/d));
					
	}
      }

      //if (d > 1) {
          
        //  for( int i=1; i < spaceF()->nLoc(3); ++i ) {
              
          //    y.getFieldPtr(i)->add(0.5,y.getField(i-1),0.5,y.getField(i+1));
          //}
      //}
   y.changed();
 }


 Teuchos::RCP<const SpaceT> spaceC() const { return(op_->spaceC()); };
 Teuchos::RCP<const SpaceT> spaceF() const { return(op_->spaceF()); };


 Teuchos::RCP<OperatorT> getOperatorPtr() { return( op_ ); }
	
};

template< class OpT >
Teuchos::RCP< RestrictionTimeOp<OpT> > createRestrictionTimeOp(
    const Teuchos::RCP<const typename OpT::FSpaceT>& spaceF, const Teuchos::RCP<const typename OpT::CSpaceT>& spaceC ) {
	return( Teuchos::rcp( new RestrictionTimeOp<OpT>( spaceF,spaceC ) ) );
	}

} // end of namespace Pimpact


#endif


