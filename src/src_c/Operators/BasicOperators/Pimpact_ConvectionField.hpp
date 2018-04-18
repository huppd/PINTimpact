/// Pimpact 
/// \author huppd
/// \date 2018


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

  using GridT = ST;

  using Scalar = typename GridT::Scalar;
  using Ordinal = typename GridT::Ordinal;

  static const int sdim = GridT::sdim;
  static const int dimension = GridT::dimension;

  static const int dimNC = GridT::dimNC;

  using DomainFieldT = VectorField<GridT>;

  using FieldTensor = ScalarField<GridT>[3][3];

protected:

  const Teuchos::RCP<const InterpolateS2V<GridT> > interpolateS2V_;
  const Teuchos::RCP<const InterpolateV2S<Scalar, Ordinal, sdim, dimension, dimNC> > interpolateV2S_;

  FieldTensor u_;

public:


  ConvectionField(const Teuchos::RCP<const GridT>& grid):
    interpolateS2V_(create<InterpolateS2V>(grid)),
    interpolateV2S_(createInterpolateV2S(grid)),
    u_{
    { { grid, Owning::Y, F::U } ,
      { grid, Owning::Y, F::U } ,
      { grid, Owning::Y, F::U } }, {
      { grid, Owning::Y, F::V } ,
      { grid, Owning::Y, F::V } ,
      { grid, Owning::Y, F::V } }, {
      { grid, Owning::Y, F::W } ,
      { grid, Owning::Y, F::W } ,
      { grid, Owning::Y, F::W } }
  } {};

  ConvectionField(
    const Teuchos::RCP<const GridT>& grid,
    const Teuchos::RCP<InterpolateS2V<GridT> >& interpolateS2V,
    const Teuchos::RCP<InterpolateV2S<Scalar, Ordinal, sdim, dimension, dimNC> >& interpolateV2S):
    interpolateS2V_(interpolateS2V),
    interpolateV2S_(interpolateV2S),
    u_{
    { grid, true, F::U } ,
    { grid, true, F::U } ,
    { grid, true, F::U } ,
    { grid, true, F::V } ,
    { grid, true, F::V } ,
    { grid, true, F::V } ,
    { grid, true, F::W } ,
    { grid, true, F::W } ,
    { grid, true, F::W }
  } {};


  void assignField(const VectorField<GridT>& mv) {

    ScalarField<GridT> temp(mv.grid());

    for(int i=0; i<GridT::sdim; ++i) {
      interpolateV2S_->apply(mv(static_cast<F>(i)), temp);
      for(int j=0; j<GridT::sdim; ++j) {
        interpolateS2V_->apply(temp, u_[j][i]);
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
