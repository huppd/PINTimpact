#pragma once
#ifndef PIMPACT_SCALINGOP_HPP
#define PIMPACT_SCALINGOP_HPP


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "BelosTypes.hpp"




namespace Pimpact {


/// \brief transposes operator.
///
/// the \c DomainFieldT  \c RangeFieldT is to be equal for both \c OP.
/// \ingroup Operator
template<class FT>
class ScalingOp {

public:

  using DomainFieldT = FT;
  using RangeFieldT = FT;

  using SpaceT = typename DomainFieldT::SpaceT;

protected:

  Teuchos::RCP<FT> field_;

public:

  ScalingOp( const Teuchos::RCP<FT>& field ): field_(field) {};


  void apply(const DomainFieldT& x, RangeFieldT& y, const Belos::ETrans
      trans=Belos::NOTRANS ) const {

    y = x;
    y.scale(*field_);
  }

  void assignField( const DomainFieldT& mv ) {};

  constexpr const Teuchos::RCP<const SpaceT>& space() const {
    return field_->space();
  };

  void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {}

  bool hasApplyTranspose() const {
    return true;
  }

  const std::string getLabel() const {
    return "ScalingOp";
  };

  void print( std::ostream& out=std::cout ) const {
    out << getLabel() << ":\n";
    //op_->print( out );
  }

}; // end of class ScalingOp



/// \relates ScalingOp
template<class FT>
Teuchos::RCP< ScalingOp<FT> > createScalingOp(
    const Teuchos::RCP<FT>& field ) {

  return Teuchos::rcp( new ScalingOp<FT>(field) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SCALINGOP_HPP