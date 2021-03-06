/// Pimpact 
/// \author huppd
/// \date 2018


#pragma once
#ifndef PIMPACT_VECTORFIELD_HPP
#define PIMPACT_VECTORFIELD_HPP


#include <iostream>
#include <vector>

#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_AbstractField.hpp"
#include "Pimpact_extern_ScalarField.hpp"
#include "Pimpact_extern_VectorField.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \brief important basic Vector class  it wraps three ScalarFields.
/// \ingroup Field
/// \relates ScalarField
template<class GridType>
class VectorField : protected AbstractField<GridType> {

public:

  using GridT = GridType;

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  using ScalarArray = ST*;

  using SF = ScalarField<GridT>;

  ScalarArray s_;

  const Owning owning_; /// <not template parameter

  SF sFields_[3];

  void allocate() {
    OT n = getStorageSize();
    setStoragePtr(new ST[n]);
    std::uninitialized_fill_n(s_, n , 0.);
  }

public:

  VectorField(const Teuchos::RCP<const GridT>& grid, const Owning owning=Owning::Y):
    AbstractField<GridT>(grid),
    owning_(owning),
    sFields_{
      {grid, Owning::N, F::U},
      {grid, Owning::N, F::V},
      {grid, Owning::N, F::W} } {

    if(owning_==Owning::Y) allocate();
  };


  /// \brief copy constructor.
  ///
  /// shallow copy, because of efficiency and conistency with \c Pimpact::MultiField
  /// \param vF
  /// \param copyType by default a ECopy::Shallow is done but allows also to deepcopy the field
  VectorField(const VectorField& vF, const ECopy copyType=ECopy::Deep):
    AbstractField<GridT>(vF.grid()),
    owning_(vF.owning_),
    sFields_{ {vF(F::U), copyType}, {vF(F::V), copyType}, {vF(F::W), copyType} } {

    if(owning_==Owning::Y) {

      allocate();

      switch(copyType) {
      case ECopy::Shallow:
        break;
      case ECopy::Deep:
        *this = vF;
        break;
      }
    }
  };


  ~VectorField() {
    if(owning_==Owning::Y) delete[] s_;
  }

  Teuchos::RCP<VectorField> clone(const ECopy copyType=ECopy::Deep) const {

    Teuchos::RCP<VectorField> vf = Teuchos::rcp(new VectorField(grid()));

    switch(copyType) {
      case ECopy::Shallow:
        break;
      case ECopy::Deep:
        *vf = *this;
        break;
    }

    return vf;
  }

  /// \name Attribute methods
  /// @{


  /// \brief returns the length of Field.
  ///
  /// the vector length is with regard to the inner points such that
  /// \f[ N_u = (N_x-1)(N_y-2)(N_z-2) \f]
  /// \f[ N_v = (N_x-2)(N_y-1)(N_z-2) \f]
  /// \f[ N_w = (N_x-2)(N_y-2)(N_z-1) \f]
  /// \return vect length \f[= N_u+N_v+N_w\f]
  constexpr OT getLength() {
    OT n = 0;
    for(F i=F::U; i<GridT::sdim; ++i)
      n += at(i).getLength();

    return n;
  }


  /// @}
  /// \name Update methods
  /// @{

  /// \brief Replace \c this with \f$\alpha a + \beta b\f$.
  ///
  /// only inner points
  void add(const ST alpha, const VectorField& a, const ST beta, const
            VectorField& b, const B wB=B::Y) {

    // add test for consistent VectorGrids in debug mode
    for(F i=F::U; i<GridT::sdim; ++i)
      at(i).add(alpha, a(i), beta, b(i), wB);

    changed();
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1, \dots, n \f]
  /// \return Reference to this object
  void abs(const VectorField& y, const B bcYes=B::Y) {
    for(F i=F::U; i<GridT::sdim; ++i)
      at(i).abs(y(i), bcYes);
    changed();
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1, \dots, n  \f]
  /// \return Reference to this object
  void reciprocal(const VectorField& y, const B bcYes=B::Y) {
    // add test for consistent VectorGrids in debug mode
    for(F i=F::U; i<GridT::sdim; ++i)
      at(i).reciprocal(y(i), bcYes);
    changed();
  }


  /// \brief Scale each element of the vectors in \c this with \c alpha.
  void scale(const ST alpha, const B bcYes=B::Y) {
    for(F i=F::U; i<GridT::sdim; ++i)
      at(i).scale(alpha, bcYes);
    changed();
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot a_i \quad \mbox{for } i=1, \dots, n \f]
  /// \return Reference to this object
  void scale(const VectorField& a, const B bcYes=B::Y) {
    // add test for consistent VectorGrids in debug mode
    for(F i=F::U; i<GridT::sdim; ++i)
      at(i).scale(a(i), bcYes);
    changed();
  }


  /// \brief Compute a scalar \c b, which is the dot-product of \c a and \c this, i.e.\f$b = a^H this\f$.
  constexpr ST dotLoc (const VectorField& a, const B bcYes=B::Y) const {
    ST b = 0.;

    for(F i=F::U; i<GridT::sdim; ++i)
      b += at(i).dotLoc(a(i), bcYes);

    return b;
  }


  /// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
  constexpr ST dot(const VectorField& y, const B bcYes=B::Y) const {

    return this->reduce(comm(), dotLoc(y, bcYes));
  }


  /// @}
  /// \name Norm method
  /// @{

  constexpr ST normLoc(ENorm type = ENorm::Two, const B bcYes=B::Y) const {

    ST normvec = 0.;

    for(F i=F::U; i<GridT::sdim; ++i)
      normvec =
        (type==ENorm::Inf)?
        std::max(at(i).normLoc(type, bcYes), normvec):
        (normvec+at(i).normLoc(type, bcYes));

    return normvec;
  }


/// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
  constexpr ST norm(const ENorm type = ENorm::Two, const B bcYes=B::Y) const {

    ST normvec = this->reduce(comm(), normLoc(type, bcYes),
                   (ENorm::Inf==type)?MPI_MAX:MPI_SUM);

    normvec = (ENorm::Two==type||ENorm::L2==type) ?
      std::sqrt(normvec) :
      normvec;

    return normvec;
  }


  /// \brief Weighted 2-Norm.
  ///
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr ST normLoc(const VectorField& weights, const B bcYes=B::Y) const {
    ST normvec = 0.;

    for(F i=F::U; i<GridT::sdim; ++i)
      normvec += at(i).normLoc(weights(i), bcYes);

    return normvec;
  }


  /// \brief Weighted 2-Norm.
  ///
  /// \warning untested
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr ST norm(const VectorField& weights, const B bcYes=B::Y) const {
    return std::sqrt(this->reduce(comm(), normLoc(weights, bcYes)));
  }


  /// @}
  /// \name Initialization methods
  /// @{


  /// \brief *this := a
  ///
  /// Assign (deep copy) a into mv.
  VectorField& operator=(const VectorField& a) {

    for(F i=F::U; i<GridT::sdim; ++i)
      at(i) = a(i);

    return *this;
  }


  /// \brief Replace the vectors with a random vectors.
  ///
  /// depending on Fortrans \c Random_number implementation, with always same
  /// seed => not save, if good randomness is required
  void random(bool useSeed=false, const B bcYes=B::Y, int seed=1) {

    for(F i=F::U; i<GridT::sdim; ++i)
      at(i).random(useSeed, bcYes, seed);

    changed();
  }

  /// \brief Replace each element of the vector  with \c alpha.
  void init(const ST alpha = Teuchos::ScalarTraits<ST>::zero(), const B bcYes=B::Y) {
    for(F i=F::U; i<GridT::sdim; ++i)
      at(i).init(alpha, bcYes);
    changed();
  }

private:

  /// \brief kind of VectorField profile
  /// \relates VectorField::initField
  enum EVectorField {
    ZeroFlow,
    PoiseuilleFlow2D_inX,
    PoiseuilleFlow2D_inY,
    PoiseuilleFlow2D_inZ,
    Pulsatile2D_inXC,
    Pulsatile2D_inXS,
    Pulsatile2D_inYC,
    Pulsatile2D_inYS,
    Circle2D,
    Circle2D_inXZ,
    RankineVortex2D,
    GaussianForcing1D,
    BoundaryFilter1D,
    GaussianForcing2D,
    BoundaryFilter2D,
    Streaming2DC,
    Streaming2DS,
    VPoint2D,
    Disc2D,
    RotationDisc2D,
    ConstFlow,
    SweptHiemenzFlow,
    Disturbance,
    ScalarFields,
    Couette,
    Cavity
  };

  /// \brief helper function getting \c EVectorField for switch statement from name
  ///
  /// \param name input name
  /// \return according int number
  EVectorField string2enum(const std::string& name) {

    std::string lcName = name;
    std::transform(lcName.begin(), lcName.end(), lcName.begin(), ::tolower);

    if("zero" == lcName) return ZeroFlow;
    else if("constant" == lcName) return ConstFlow;
    else if("poiseuille" == lcName) return PoiseuilleFlow2D_inX;
    else if("poiseuille in x" == lcName) return PoiseuilleFlow2D_inX;
    else if("poiseuille in y" == lcName) return PoiseuilleFlow2D_inY;
    else if("poiseuille in z" == lcName) return PoiseuilleFlow2D_inZ;
    else if("pulsatile in x cos" == lcName) return Pulsatile2D_inXC;
    else if("pulsatile in x sin" == lcName) return Pulsatile2D_inXS;
    else if("pulsatile in y cos" == lcName) return Pulsatile2D_inYC;
    else if("pulsatile in y sin" == lcName) return Pulsatile2D_inYS;
    else if("streaming" == lcName) return Streaming2DS;
    else if("circle" == lcName) return Circle2D;
    else if("circle xy" == lcName) return Circle2D;
    else if("circle xz" == lcName) return Circle2D_inXZ;
    else if("rankine vortex" == lcName) return RankineVortex2D;
    else if("gaussian forcing 1d" == lcName) return GaussianForcing1D;
    else if("boundary filter 1d" == lcName) return BoundaryFilter1D;
    else if("gaussian forcing 2d" == lcName) return GaussianForcing2D;
    else if("boundary filter 2d" == lcName) return BoundaryFilter2D;
    else if("streaming 2d cos" == lcName) return Streaming2DC;
    else if("streaming 2d sin" == lcName) return Streaming2DS;
    else if("v point 2d" == lcName) return VPoint2D;
    else if("disc 2d" == lcName) return Disc2D;
    else if("rotation disc 2d" == lcName) return RotationDisc2D;
    else if("shbl" == lcName) return SweptHiemenzFlow;
    else if("swept hiemenz flow" == lcName) return SweptHiemenzFlow;
    else if("disturbance" == lcName) return Disturbance;
    else if("scalar" == lcName) return ScalarFields;
    else if("couette" == lcName) return Couette;
    else if("cavity" == lcName) return Cavity;
    else {
      const bool& Flow_Type_not_known = true;
      TEUCHOS_TEST_FOR_EXCEPT(Flow_Type_not_known);
    }
    return ZeroFlow; // just to please the compiler
  }

public:

  /// \brief initializes including boundaries to zero
  /// \todo rm FORTRAN
  /// \todo make it init or addable
  void initField(Teuchos::ParameterList& para, const Add add=Add::N) {

    EVectorField type =
      string2enum(para.get<std::string>("Type", "zero"));

    switch(type) {
      case ZeroFlow : {
        for(F i=F::U; i<GridT::sdim; ++i)
          if(Add::N==add) at(i).init();
        break;
      }
      case ConstFlow : {
        ST u = para.get<ST>("U", 1.);
        ST v = para.get<ST>("V", 1.);
        ST w = para.get<ST>("W", 1.);

        at(F::U).initFromFunction([&u](ST x, ST y, ST z)->ST{
            return u; },
            add);
        at(F::V).initFromFunction([&v](ST x, ST y, ST z)->ST{
            return v; },
            add);
        at(F::W).initFromFunction([&w](ST x, ST y, ST z)->ST{
            return w; },
            add);
        break;
      }
      case PoiseuilleFlow2D_inX : {
        for(F i=F::U; i<GridT::sdim; ++i)
          if(F::U==i)
            at(i).initFromFunction(
                [] (ST x, ST y, ST z)->ST { return 4.*y*(1.-y); },
                add);
          else if(Add::N==add) at(i).init();
        break;
      }
      case PoiseuilleFlow2D_inY : {
        for(F i=F::U; i<GridT::sdim; ++i)
          if(F::V==i)
            at(i).initFromFunction(
                [] (ST x, ST y, ST z)->ST { return 4.*x*(1.-x); },
                add);
          else if(Add::N==add) at(i).init();
        break;
      }
      case PoiseuilleFlow2D_inZ : {
        for(F i=F::U; i<GridT::sdim; ++i)
          if(F::W==i)
            at(i).initFromFunction(
                [] (ST x, ST y, ST z)->ST { return 4.*x*(1.-x); },
                add);
          else if(Add::N==add) at(i).init();
        break;
      }
      case Pulsatile2D_inXC : {
        VF_init_2DPulsatileXC(
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->sIndB(F::U),
            grid()->eIndB(F::U),
            grid()->sIndB(F::V),
            grid()->eIndB(F::V),
            grid()->sIndB(F::W),
            grid()->eIndB(F::W),
            grid()->getDomainSize()->getSize(1),
            grid()->getCoordinatesLocal()->getX(F::S, Y),
            grid()->getDomainSize()->getRe(),     // TODO: verify
            grid()->getDomainSize()->getAlpha2(), // TODO: verify
            para.get<ST>("px", 1.),          // TODO: verify
            at(F::U).getRawPtr(),
            at(F::V).getRawPtr(),
            at(F::W).getRawPtr());
        break;
      }
      case Pulsatile2D_inYC : {
        VF_init_2DPulsatileYC(
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->sIndB(F::U),
            grid()->eIndB(F::U),
            grid()->sIndB(F::V),
            grid()->eIndB(F::V),
            grid()->sIndB(F::W),
            grid()->eIndB(F::W),
            grid()->getDomainSize()->getSize(0),
            grid()->getCoordinatesLocal()->getX(F::S, X),
            grid()->getDomainSize()->getRe(),     // TODO: verify
            grid()->getDomainSize()->getAlpha2(), // TODO: verify
            para.get<ST>("px", 1.),          // TODO: verify
            at(F::U).getRawPtr(),
            at(F::V).getRawPtr(),
            at(F::W).getRawPtr());
        break;
      }
      case Pulsatile2D_inXS : {
        VF_init_2DPulsatileXS(
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->sIndB(F::U),
            grid()->eIndB(F::U),
            grid()->sIndB(F::V),
            grid()->eIndB(F::V),
            grid()->sIndB(F::W),
            grid()->eIndB(F::W),
            grid()->getDomainSize()->getSize(1),
            grid()->getCoordinatesLocal()->getX(F::S, Y),
            grid()->getDomainSize()->getRe(),     // TODO: verify
            grid()->getDomainSize()->getAlpha2(), // TODO: verify
            para.get<ST>("px", 1.),          // TODO: verify
            at(F::U).getRawPtr(),
            at(F::V).getRawPtr(),
            at(F::W).getRawPtr());
        break;
      }
      case Pulsatile2D_inYS : {
        VF_init_2DPulsatileYS(
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->sIndB(F::U),
            grid()->eIndB(F::U),
            grid()->sIndB(F::V),
            grid()->eIndB(F::V),
            grid()->sIndB(F::W),
            grid()->eIndB(F::W),
            grid()->getDomainSize()->getSize(0),
            grid()->getCoordinatesLocal()->getX(F::S, X),
            grid()->getDomainSize()->getRe(),     // TODO: verify
            grid()->getDomainSize()->getAlpha2(), // TODO: verify
            para.get<ST>("px", 1.),          // TODO: verify
            at(F::U).getRawPtr(),
            at(F::V).getRawPtr(),
            at(F::W).getRawPtr());
        break;
      }
      case Streaming2DC : {
        ST amp = grid()->getDomainSize()->getRe();
        ST pi = 4.*std::atan(1.);
        //ST L1 = grid()->getDomainSize()->getSize(X);
        ST om = 2.*pi;

        if(Add::N==add) at(F::U).init();
        at(F::V).initFromFunction(
            [&amp, &om](ST x, ST y, ST z) -> ST {
            return amp*std::cos(om*x); },
            add);
        if(Add::N==add) at(F::W).init();
        break;
      }
      case Streaming2DS: {
        //ST amp = grid()->getDomainSize()->getRe();
        ST amp = 1.;
        ST pi = 4.*std::atan(1.);
        //ST L1 = grid()->getDomainSize()->getSize(X);
        ST om = 2.*pi;

        if(Add::N==add) at(F::U).init();
        at(F::V).initFromFunction(
            [&amp, &om](ST x, ST y, ST z) -> ST {
            return amp*std::sin(om*x); },
            add);
        if(Add::N==add) at(F::W).init();
        break;
      }
      case Circle2D : {
        at(F::U).initField(Grad2D_inY, -1.);
        at(F::V).initField(Grad2D_inX,  1.);
        if(Add::N==add) at(F::W).init();
        break;
      }
      case Circle2D_inXZ : {
        at(F::U).initField(Grad2D_inY, -1.);
        if(Add::N==add) at(F::V).init();
        at(F::W).initField(Grad2D_inX,  1.);
        break;
      }
      case RankineVortex2D : {
        VF_init_RankineVortex(
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->sIndB(F::U),
            grid()->eIndB(F::U),
            grid()->sIndB(F::V),
            grid()->eIndB(F::V),
            grid()->sIndB(F::W),
            grid()->eIndB(F::W),
            grid()->getDomainSize()->getSize(),
            grid()->getCoordinatesLocal()->getX(F::S, X),
            grid()->getCoordinatesLocal()->getX(F::S, Y),
            grid()->getCoordinatesLocal()->getX(F::U, X),
            grid()->getCoordinatesLocal()->getX(F::V, Y),
            at(F::U).getRawPtr(),
            at(F::V).getRawPtr(),
            at(F::W).getRawPtr());
        break;
      }
      case GaussianForcing1D : {
        VF_init_GaussianForcing1D(
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->sIndB(F::U),
            grid()->eIndB(F::U),
            grid()->sIndB(F::V),
            grid()->eIndB(F::V),
            grid()->sIndB(F::W),
            grid()->eIndB(F::W),
            grid()->getDomainSize()->getSize(0),
            grid()->getCoordinatesLocal()->getX(F::U, X),
            at(F::U).getRawPtr(),
            at(F::V).getRawPtr(),
            at(F::W).getRawPtr());
        break;
      }
      case BoundaryFilter1D : {
        VF_init_BoundaryFilter1D(
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->sIndB(F::U),
            grid()->eIndB(F::U),
            grid()->sIndB(F::V),
            grid()->eIndB(F::V),
            grid()->sIndB(F::W),
            grid()->eIndB(F::W),
            grid()->getDomainSize()->getSize(0),
            grid()->getCoordinatesLocal()->getX(F::U, X),
            at(F::U).getRawPtr(),
            at(F::V).getRawPtr(),
            at(F::W).getRawPtr());
        break;
      }
      case GaussianForcing2D : {
        VF_init_GaussianForcing2D(
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->sIndB(F::U),
            grid()->eIndB(F::U),
            grid()->sIndB(F::V),
            grid()->eIndB(F::V),
            grid()->sIndB(F::W),
            grid()->eIndB(F::W),
            grid()->getDomainSize()->getSize(),
            grid()->getCoordinatesLocal()->getX(F::S, X),
            grid()->getCoordinatesLocal()->getX(F::S, Y),
            grid()->getCoordinatesLocal()->getX(F::U, X),
            grid()->getCoordinatesLocal()->getX(F::V, Y),
            at(F::U).getRawPtr(),
            at(F::V).getRawPtr(),
            at(F::W).getRawPtr());
        break;
      }
      case BoundaryFilter2D : {
        VF_init_BoundaryFilter2D(
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->sIndB(F::U),
            grid()->eIndB(F::U),
            grid()->sIndB(F::V),
            grid()->eIndB(F::V),
            grid()->sIndB(F::W),
            grid()->eIndB(F::W),
            grid()->getDomainSize()->getSize(),
            grid()->getCoordinatesLocal()->getX(F::S, X),
            grid()->getCoordinatesLocal()->getX(F::S, Y),
            grid()->getCoordinatesLocal()->getX(F::U, X),
            grid()->getCoordinatesLocal()->getX(F::V, Y),
            at(F::U).getRawPtr(),
            at(F::V).getRawPtr(),
            at(F::W).getRawPtr());
        break;
      }
      case VPoint2D : {
        VF_init_Vpoint(
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->sIndB(F::U),
            grid()->eIndB(F::U),
            grid()->sIndB(F::V),
            grid()->eIndB(F::V),
            grid()->sIndB(F::W),
            grid()->eIndB(F::W),
            grid()->getDomainSize()->getSize(),
            grid()->getCoordinatesLocal()->getX(F::U, X),
            grid()->getCoordinatesLocal()->getX(F::S, Y),
            grid()->getDomainSize()->getRe(),     // TODO: verify
            at(F::U).getRawPtr(),
            at(F::V).getRawPtr(),
            at(F::W).getRawPtr());
        break;
      }
      case Disc2D : {
        VF_init_Disc(
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->sIndB(F::U),
            grid()->eIndB(F::U),
            grid()->sIndB(F::V),
            grid()->eIndB(F::V),
            grid()->sIndB(F::W),
            grid()->eIndB(F::W),
            grid()->getCoordinatesLocal()->getX(F::S, X),
            grid()->getCoordinatesLocal()->getX(F::S, Y),
            grid()->getCoordinatesLocal()->getX(F::S, Z),
            grid()->getCoordinatesLocal()->getX(F::U, X),
            grid()->getCoordinatesLocal()->getX(F::V, Y),
            //re, om, px, sca,
            para.get<ST>("center x", 1.),
            para.get<ST>("center y", 1.),
            para.get<ST>("radius", 1.),
            para.get<ST>("sca", 0.1),
            at(F::U).getRawPtr(),
            at(F::V).getRawPtr(),
            at(F::W).getRawPtr());
        break;
      }
      case RotationDisc2D : {
        VF_init_RotatingDisc(
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->sIndB(F::U),
            grid()->eIndB(F::U),
            grid()->sIndB(F::V),
            grid()->eIndB(F::V),
            grid()->sIndB(F::W),
            grid()->eIndB(F::W),
            grid()->getCoordinatesLocal()->getX(F::S, X),
            grid()->getCoordinatesLocal()->getX(F::S, Y),
            para.get<ST>("center x", 1.),
            para.get<ST>("center y", 1.),
            para.get<ST>("omega", 1.),
            at(F::U).getRawPtr(),
            at(F::V).getRawPtr(),
            at(F::W).getRawPtr());
        break;
      }
      case SweptHiemenzFlow : {
        ST pi = 4.*std::atan(1.);

        OT nTemp = grid()->gu(X) - grid()->gl(X) + 1;
        //for(OT i=0; i<=grid()->nLoc(X); ++i) {
        //std::cout << "i: " << i<< " (\t";
        //for(OT ii=0; ii<nTemp; ++ii)
        //std::cout << grid()->getInterpolateV2S()->getC(X)[ i*nTemp + ii ] << ", \t";
        //std::cout << ")\n";
        //}
        //ST c[6] = {
        //grid()->getInterpolateV2S()->getC(X)[ nTemp + 0 ],
        //grid()->getInterpolateV2S()->getC(X)[ nTemp + 1 ],
        //grid()->getInterpolateV2S()->getC(X)[ nTemp + 2 ],
        //grid()->getInterpolateV2S()->getC(X)[ nTemp + 3 ],
        //grid()->getInterpolateV2S()->getC(X)[ nTemp + 4 ],
        //grid()->getInterpolateV2S()->getC(X)[ nTemp + 5 ] };
        //std::cout << "c\n";
        //for(OT ii=0; ii<nTemp; ++ii)
        //std::cout << c[  ii ] << ", \t";
        //std::cout << "c\n";

        VF_init_SHBF(
            //1,
            grid()->rankST(),
            grid()->getShift(0),
            grid()->getProcGrid()->getIB(0),
            grid()->nGlo(),
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->dl(),
            grid()->du(),
            grid()->sIndB(F::U),
            grid()->eIndB(F::U),
            grid()->sIndB(F::V),
            grid()->eIndB(F::V),
            grid()->sIndB(F::W),
            grid()->eIndB(F::W),
            grid()->getCoordinatesGlobal()->getX(F::S, X),
            grid()->getCoordinatesGlobal()->getX(F::U, X),
            grid()->getCoordinatesLocal()->getX(F::W, Z),
            grid()->getInterpolateV2S()->getC(X)+2*nTemp-1, /// \todo rm dirty hack
            grid()->getDomainSize()->getRe(),
            para.get<int>("nonDim", 0),
            para.get<ST>("kappa", 0.),
            para.get<ST>("sweep angle", 0.),
            at(F::U).getRawPtr(),
            at(F::V).getRawPtr(),
            at(F::W).getRawPtr());

        //std::cout << "hello\n" << grid()->getInterpolateV2S()->getC(X)[9] <<
        //"\n";

        break;
      }
      case Disturbance : {
        ST pi = 4.*std::atan(1.);
        ST xc = para.get<ST>("xc", 3.0);
        ST zc = para.get<ST>("zc", 1.5);
        ST b  = para.get<ST>("b",  3.);
        ST A  = para.get<ST>("A",  0.1);

        Teuchos::RCP<const DomainSize<ST, GridT::sdim> > domain = grid()->getDomainSize();
        Teuchos::RCP<const CoordinatesLocal<ST, OT, GridT::dimension, GridT::dimNC> > coord =
          grid()->getCoordinatesLocal();

        const B bY = B::Y;
        if(0<grid()->bcl(Y)) {
          OT j=grid()->si(F::U, Y, bY);
          for(OT k=grid()->si(F::U, Z, bY); k<=grid()->ei(F::U, Z, bY); ++k)
            for(OT i=grid()->si(F::U, X, bY); i<=grid()->ei(F::U, X, bY); ++i) {
              ST x = coord->getX(F::U, X, i);
              ST z = coord->getX(F::U, Z, k);

              ST u = 0.;
              if(std::fabs(x-xc)<b && std::fabs(z-zc)<b)
                u += 0.5*A*std::sin(pi*(z-zc)/b)*(1. + std::cos(pi*(x-xc)/b));
              if(std::fabs(x-xc)<b && std::fabs(z+zc)<b)
                u -= 0.5*A*std::sin(pi*(z+zc)/b)*(1. + std::cos(pi*(x-xc)/b));
              if(Add::Y==add)
                at(F::U)(i, j, k) += u;
              else
                at(F::U)(i, j, k) = u;
            }
        }
        at(F::U).changed();

        if(0<grid()->bcl(Y)) {
          OT j=grid()->si(F::W, Y, bY);
          for(OT k=grid()->si(F::W, Z, bY); k<=grid()->ei(F::W, Z, bY); ++k)
            for(OT i=grid()->si(F::W, X, bY); i<=grid()->ei(F::W, X, bY); ++i) {
              ST x = coord->getX(F::W, X, i);
              ST z = coord->getX(F::W, Z, k);

              ST w = 0.;
              if(std::fabs(x-xc)<b && std::fabs(z-zc)<b)
                w -= 0.5*A*std::sin(pi*(x-xc)/b)*(1. + std::cos(pi*(z-zc)/b));
              if(std::fabs(x-xc)<b && std::fabs(z+zc)<b)
                w += 0.5*A*std::sin(pi*(x-xc)/b)*(1. + std::cos(pi*(z+zc)/b));
              if(Add::Y==add)
                at(F::W)(i, j, k) += w;
              else
                at(F::W)(i, j, k) = w;
            }
        }
        at(F::W).changed();

        //at(F::W).initFromFunction(
            //[=](ST x_, ST y_, ST z_) -> ST {

              //ST x = x_*domain->getSize(X) + domain->getOrigin(X);
              //ST y = y_*domain->getSize(Y) + domain->getOrigin(Y);
              //ST z = z_*domain->getSize(Z) + domain->getOrigin(Z);

              //ST w = 0.;
              //if(y<=Teuchos::ScalarTraits<ST>::eps() && std::fabs(x-xc)<b && std::fabs(z-zc)<b)
                //w -= 0.5*A*std::sin(pi*(x-xc)/b)*(1. + std::cos(pi*(z-zc)/b));
              //if(y<=Teuchos::ScalarTraits<ST>::eps() && std::fabs(x-xc)<b && std::fabs(z+zc)<b)
                //w += 0.5*A*std::sin(pi*(x-xc)/b)*(1. + std::cos(pi*(z+zc)/b));
              //return w; },
              //add);
        //at(F::U).initFromFunction(
            //[=](ST x_, ST y_, ST z_) -> ST {

              //ST x = x_*domain->getSize(X) + domain->getOrigin(X);
              //ST y = y_*domain->getSize(Y) + domain->getOrigin(Y);
              //ST z = z_*domain->getSize(Z) + domain->getOrigin(Z);

              //ST u = 0.;
              //if(y<=Teuchos::ScalarTraits<ST>::eps() && std::fabs(x-xc)<b && std::fabs(z-zc)<b)
                //u += 0.5*A*std::sin(pi*(z-zc)/b)*(1. + std::cos(pi*(x-xc)/b));
              //if(y<=Teuchos::ScalarTraits<ST>::eps() && std::fabs(x-xc)<b && std::fabs(z+zc)<b)
                //u -= 0.5*A*std::sin(pi*(z+zc)/b)*(1. + std::cos(pi*(x-xc)/b));
              //return u; },
            //add);
        if(Add::N==add) at(F::V).init();
        break;
      }
      case ScalarFields : {
        at(F::U).initField(para.sublist("U"), add);
        at(F::V).initField(para.sublist("V"), add);
        at(F::W).initField(para.sublist("W"), add);
        break;
      }
      case Couette : {
        at(F::U).initFromFunction(
            [](ST x, ST y, ST z) -> ST {
            return y; },
            add);
        if(Add::N==add) at(F::V).init();
        if(Add::N==add) at(F::W).init();
        break;
      }
      case Cavity : {
        at(F::U).initFromFunction(
            [](ST x, ST y, ST z) -> ST {
            if(std::fabs(x-1.)<0.5)
            return 1.;
            else
            return 0.; },
            add);
        if(Add::N==add) at(F::V).init();
        if(Add::N==add) at(F::W).init();
        break;
      }
    }
    changed();
  }


  /// \brief extrapolates on the boundaries such that it is zero
  /// \note dirty hack(necessary for TripleCompostion)
  void extrapolateBC(const Belos::ETrans trans=Belos::NOTRANS) {
    for(F i=F::U; i<GridT::sdim; ++i)
      at(i).extrapolateBC(trans);
  }


  /// dirty hack(necessary for MG)
  void level() const {}

  /// @}


  /// \brief Print the vector.  To be used for debugging only.
  void print(std::ostream& out=std::cout) const {
    for(F i=F::U; i<GridT::sdim; ++i)
      at(i).print(out);
  }


  /// \brief writes vector field to hdf5 file
  ///
  /// \param count
  /// \param restart decides if velocity is interpolated to pressure points
  /// \todo implement/test restart and read
  void write(const int count=0, const bool restart=false) const {

    if(0==grid()->rankS()) {
      Teuchos::Tuple<OT, 3> N;
      for(int i=0; i<3; ++i) {
        N[i] = grid()->nGlo(i);
        if(grid()->getBCGlobal()->getBCL(i)==Pimpact::BC::Periodic)
          N[i] = N[i]-1;
      }
      if(!restart) {
        std::ofstream xfile;
        std::ostringstream ss;
        ss << std::setw(5) << std::setfill('0') << count;
        //        std::string fname = "v_"+ss.str();
        xfile.open("vel_"+ ss.str() +".xmf", std::ofstream::out);
        xfile<< "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.1\">\n";
        xfile << "\t<Domain>\n";
        xfile << "\t\t<Grid Name=\"3DRectMesh\" GridType=\"Uniform\">\n";
        xfile << "\t\t\t<Topology TopologyType=\"3DRectMesh\" Dimensions=\""<< N[2] << " " << N[1] << " " << N[0] << "\"/>\n";
        xfile << "\t\t\t<Geometry GeometryType=\"VXVYVZ\">\n";
        xfile << "\t\t\t\t<DataItem ItemType=\"Uniform\"\n";
        xfile << "\t\t\t\t\tDimensions=\""<< N[0] << "\"\n";
        xfile << "\t\t\t\t\tNumberType=\"Float\"\n";
        xfile << "\t\t\t\t\tPrecision=\"8\"\n";
        xfile << "\t\t\t\t\tFormat=\"HDF\">\n";
        xfile << "\t\t\t\t\tvelX_"<< ss.str() << ".h5:/VectorX\n";
        xfile << "\t\t\t\t</DataItem>\n";
        xfile << "\t\t\t\t<DataItem ItemType=\"Uniform\"\n";
        xfile << "\t\t\t\t\tDimensions=\""<< N[1] << "\"\n";
        xfile << "\t\t\t\t\tNumberType=\"Float\"\n";
        xfile << "\t\t\t\t\tPrecision=\"8\"\n";
        xfile << "\t\t\t\t\tFormat=\"HDF\">\n";
        xfile << "\t\t\t\t\tvelX_"<< ss.str() << ".h5:/VectorY\n";
        xfile << "\t\t\t\t</DataItem>\n";
        xfile << "\t\t\t\t<DataItem ItemType=\"Uniform\"\n";
        xfile << "\t\t\t\t\tDimensions=\""<< N[2] << "\"\n";
        xfile << "\t\t\t\t\tNumberType=\"Float\"\n";
        xfile << "\t\t\t\t\tPrecision=\"8\"\n";
        xfile << "\t\t\t\t\tFormat=\"HDF\">\n";
        xfile << "\t\t\t\t\tvelX_"<< ss.str() << ".h5:/VectorZ\n";
        xfile << "\t\t\t\t</DataItem>\n";
        xfile << "\t\t\t</Geometry>\n";
        xfile << "\t\t\t<Attribute Name=\"VelX\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        xfile << "\t\t\t\t<DataItem Dimensions=\""<< N[2] << " " << N[1] << " " << N[0] << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
        xfile << "\t\t\t\t\tvelX_"<< ss.str() << ".h5:/velX\n";
        xfile << "\t\t\t\t</DataItem>\n";
        xfile << "\t\t\t</Attribute>\n";
        xfile << "\t\t\t<Attribute Name=\"VelY\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        xfile << "\t\t\t\t<DataItem Dimensions=\""<< N[2] << " " << N[1] << " " << N[0] << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
        xfile << "\t\t\t\t\tvelY_"<< ss.str() << ".h5:/velY\n";
        xfile << "\t\t\t\t</DataItem>\n";
        xfile << "\t\t\t</Attribute>\n";
        xfile << "\t\t\t<Attribute Name=\"VelZ\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        xfile << "\t\t\t\t<DataItem Dimensions=\""<< N[2] << " " << N[1] << " " << N[0] << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
        xfile << "\t\t\t\t\tvelZ_"<< ss.str() << ".h5:/velZ\n";
        xfile << "\t\t\t\t</DataItem>\n";
        xfile << "\t\t\t</Attribute>\n";
        xfile << "\t\t</Grid>\n";
        xfile << "\t</Domain>\n";
        xfile << "</Xdmf>\n";
        xfile.close();
      }
    }
    for(F i=F::U; i<GridT::sdim; ++i)
      at(i).write(count, restart);
  }

  void read(const int count=0) {
    for(F i=F::U; i<GridT::sdim; ++i)
      at(i).read(count);
  }


public:

  /// \name storage methods.
  /// \brief highly dependent on underlying storage should only be used by
  /// Operator or on top field implementer.
  /// @{

  constexpr OT getStorageSize() const {
    return sFields_[0].getStorageSize()*3;
  }

  void setStoragePtr(ST*  array) {
    s_ = array;
    OT n = sFields_[0].getStorageSize();
    for(int i=0; i<3; ++i)
      sFields_[i].setStoragePtr(s_+i*n);
  }

  constexpr ST* getRawPtr() {
    return s_;
  }

  constexpr const ST* getConstRawPtr() const {
    return s_;
  }

  /// \deprecated
  ST* getRawPtr (int i)       {
    return sFields_[i].getRawPtr();
  }

  /// \deprecated
  constexpr const ST* getConstRawPtr (int i)  const  {
    return sFields_[i].getConstRawPtr();
  }

  ///  @}


  SF& operator()(F i) {
    return at(i);
  }

  constexpr const SF& operator()(F i) const {
    return at(i);
  }

protected:

  SF& at(F i) {
    assert(!(F::S==i));
    return sFields_[ static_cast<int>(i) ];
  }

  constexpr const SF& at(F i) const {
    assert(!(F::S==i));
    return sFields_[ static_cast<int>(i) ];
  }

public:

  constexpr const Teuchos::RCP<const GridT>& grid() {
    return AbstractField<GridT>::grid_;
  }

  constexpr const MPI_Comm& comm() {
    return grid()->comm();
  }


  constexpr ST allReduce(ST local, const MPI_Op& op=MPI_SUM) {
    return this->reduce(comm(), local);
  }

  /// \name comunication methods.
  /// \brief highly dependent on underlying storage should only be used by Operator or on top field implementer.
  ///
  /// @{

  void changed(const int vel_dir, const int dir) const {
    sFields_[vel_dir].changed(dir);
  }

  void changed() const {
    for(int vel_dir=0; vel_dir<GridT::sdim; ++vel_dir)
      for(int dir=0; dir<GridT::sdim; ++dir) {
        changed(vel_dir, dir);
      }
  }

  void exchange(const int vel_dir, const int dir) const {
    sFields_[vel_dir].exchange(dir);
  }

  void exchange() const {
    for(int vel_dir=0; vel_dir<GridT::sdim; ++vel_dir)
      for(int dir=0; dir<GridT::sdim; ++dir)
        exchange(vel_dir, dir);
  }

  void setExchanged(const int vel_dir, const int dir) const {
    sFields_[vel_dir].setExchanged(dir);
  }

  void setExchanged() const {
    for(int vel_dir=0; vel_dir<GridT::sdim; ++vel_dir)
      for(int dir=0; dir<GridT::sdim; ++dir) {
        setExchanged(vel_dir, dir);
      }
  }

  /// @}


}; // end of class VectorField


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_VECTORFIELD_HPP
