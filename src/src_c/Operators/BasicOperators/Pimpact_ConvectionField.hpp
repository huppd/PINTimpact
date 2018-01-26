#pragma once
#ifndef PIMPACT_CONVECTIONFIELD_HPP
#define PIMPACT_CONVECTIONFIELD_HPP


#include "Pimpact_InterpolateS2VOp.hpp"
#include "Pimpact_InterpolateV2SOp.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {



/// \brief Stores the wind on differnt grid types.
/// should become template parameter for others such that interplating can be moved from assign to get(different storage needed)
/// \relates NonlinearOp
/// \todo remove RCP
template<class ST>
class ConvectionField {

public:

  using SpaceT = ST;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  static const int sdim = SpaceT::sdim;
  static const int dimension = SpaceT::dimension;

  static const int dimNC = SpaceT::dimNC;

  using DomainFieldT = VectorField<SpaceT>;

  using FieldTensor = ScalarField<SpaceT>[3][3];

protected:

  const Teuchos::RCP<const InterpolateS2V<SpaceT> > interpolateS2V_;
  const Teuchos::RCP<const InterpolateV2S<Scalar,Ordinal,sdim,dimension,dimNC> > interpolateV2S_;

  FieldTensor u_;

public:


  ConvectionField( const Teuchos::RCP<const SpaceT>& space  ):
    interpolateS2V_( create<InterpolateS2V>(space) ),
    interpolateV2S_( createInterpolateV2S( space ) ),
    u_{
    { { space, Owning::Y, F::U } ,
      { space, Owning::Y, F::U } ,
      { space, Owning::Y, F::U } },{
      { space, Owning::Y, F::V } ,
      { space, Owning::Y, F::V } ,
      { space, Owning::Y, F::V } },{
      { space, Owning::Y, F::W } ,
      { space, Owning::Y, F::W } ,
      { space, Owning::Y, F::W } }
  } {};

  ConvectionField(
    const Teuchos::RCP<const SpaceT>& space,
    const Teuchos::RCP< InterpolateS2V<SpaceT> >& interpolateS2V,
    const Teuchos::RCP< InterpolateV2S<Scalar,Ordinal,sdim,dimension,dimNC> >& interpolateV2S ):
    interpolateS2V_(interpolateS2V),
    interpolateV2S_(interpolateV2S),
    u_{
    { space,true,F::U } ,
    { space,true,F::U } ,
    { space,true,F::U } ,
    { space,true,F::V } ,
    { space,true,F::V } ,
    { space,true,F::V } ,
    { space,true,F::W } ,
    { space,true,F::W } ,
    { space,true,F::W }
  } {};


  void assignField( const VectorField<SpaceT>& mv ) {

    ScalarField<SpaceT> temp( mv.space() );

    for( int i=0; i<SpaceT::sdim; ++i ) {
      interpolateV2S_->apply( mv(static_cast<F>(i)), temp );
      for( int j=0; j<SpaceT::sdim; ++j ) {
        interpolateS2V_->apply( temp, u_[j][i] );
      }
    }
  };


  constexpr const FieldTensor& get() const {
    return u_;
  }

  //FieldTensor& get() { return u_; }


}; // end of class ConvectionField


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONFIELD_HPP
