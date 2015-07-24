#pragma once
#ifndef PIMPACT_INTERPOLATIONTIMEOP_HPP
#define PIMPACT_INTERPOLATIONTIMEOP_HPP

#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_Types.hpp"

#include "Pimpact_TimeField.hpp"



namespace Pimpact {


/// \brief Operator wrapper for \c TimeField's.
/// \ingroup TimeOperator
/// wraps and \c OperatorT and adds the functionality of handling \c TimeField.
template<class OperatorT>
class InterpolationTimeOp  {

public:

  typedef TimeField<typename OperatorT::DomainFieldT> DomainFieldT;
  typedef TimeField<typename OperatorT::RangeFieldT> RangeFieldT;

  typedef typename DomainFieldT::SpaceT SpaceT;

  typedef SpaceT FSpaceT;
  typedef SpaceT CSpaceT;

  typedef typename SpaceT::Ordinal Ordinal;

protected:

  Teuchos::RCP<OperatorT> op_; //  interpolation operator in space
    
public:
    
    InterpolationTimeOp(
                        const Teuchos::RCP<const SpaceT>& spaceC,
                        const Teuchos::RCP<const SpaceT>& spaceF):
    op_(Teuchos::rcp(new OperatorT(spaceC,spaceF)))
    {};

  /// \brief default apply
  void apply( const DomainFieldT& x,
      RangeFieldT& y) const {
      
      Ordinal d = (spaceF()->nLoc(3)) / (spaceC()->nLoc(3));
	
	//std::cout << "cr = " << spaceC()->nLoc(3) << std::endl;
	//std::cout << "fn = " << spaceF()->nLoc(3) << std::endl;		
	//std::cout << "d = " << d << std::endl;
      
     x.exchange();
      
     for( int i=0; i <= spaceC()->nLoc(3); ++i ) { 

        op_->apply( x.getConstField(i), y.getField(d*i) );
      }
     
      if (d > 1) { 
	
	for (int j = d; j > 1; j = j/2) {          
		
          for( int i=j/2; i <= spaceF()->nLoc(3) - j/2; i=i+j ) { 
              
              y.getFieldPtr(i)->add(0.5,y.getField(i-1),0.5,y.getField(i+1));
          }
      }
    } 
    y.changed(); 
  }


    Teuchos::RCP<const SpaceT> spaceC() const { return(op_->get_spaceC_()); };
    Teuchos::RCP<const SpaceT> spaceF() const { return(op_->get_spaceF_()); };


  Teuchos::RCP<OperatorT> getOperatorPtr() { return( op_ ); }
	


}; // end of class InterpolationTimeOp

    
/// \relates
template< class OpT >
Teuchos::RCP< InterpolationTimeOp<OpT> > createInterpolationTimeOp(
    const Teuchos::RCP<const typename OpT::CSpaceT>& spaceC, const Teuchos::RCP<const typename OpT::FSpaceT>& spaceF ) {

	return( Teuchos::rcp( new InterpolationTimeOp<OpT>( spaceC,spaceF ) ) );

}

} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INTERPOLATIONTIMEOP_HPP
