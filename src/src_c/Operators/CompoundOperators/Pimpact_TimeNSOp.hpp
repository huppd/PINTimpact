#pragma once
#ifndef PIMPACT_TIMENSOp_HPP
#define PIMPACT_TIMENSOp_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_CompoundField.hpp"
#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_TimeDTConvectionDiffusionOp.hpp"




namespace Pimpact {


/// \ingroup CompoundOperator
///
// \f[ \begin{bmatrix} opV2V & opS2V \\ opV2S & 0 \end{bmatrix} \mathbf{x} = \mathbf{y} \f]
template<class ST>
class TimeNSOp {

public:

  typedef typename OpV2V::DomainFieldT  VF;
  typedef typename OpS2V::DomainFieldT  SF;

  typedef typename ST SpaceT;

  typedef CompoundField< TimeField<VectorField<ST> >, TimeField<ScalarField<ST> > >  DomainFieldT;
  typedef CompoundField< TimeField<VectorField<ST> >, TimeField<ScalarField<ST> > >  RangeFieldT;


protected:

  Teuchos::RCP<const TimeDtConvectionDiffusionOp<ST> > opV2V_;
  Teuchos::RCP<const GradOp<ST> > opS2V_;
  Teuchos::RCP<const DivOp<ST> > opV2S_;

public:

  /// \todo constructor from space
  TimeNSOp(
      const Teuchos::RCP<const SpaceT>& space ):
        opV2V_( create< TimeDtConvectionDiffusionOp<ST> >( space ) ),
        opS2V_( create< DivOp <ST>(space),
        opV2S_( create< GradOp<ST>(space) ) {};

  void apply(const DomainFieldT& x, RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS  ) const {
    // H-blockz
    opV2V_->apply( x.getConstVField(), y.getVField() );
    // ~grad
    opS2V_->apply( x.getConstSField(), *temp_ );
    y.getVField().add( 1., y.getConstVField(), 1., *temp_ );
    // ~div
    opV2S_->apply( x.getConstVField(), y.getSField() );
//		y.getSField().level();
  }

  void assignField( const DomainFieldT& mv ) {
    opV2V_->assignField( mv.getConstVField() );
  };

	Teuchos::RCP<const SpaceT> space() const { return(opV2V_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "TimeNSOp " ); };

}; // end of class TimeNSOp



} // end of namespace Pimpact



#endif // end of #ifndef PIMPACT_TIMENSOp_HPP
