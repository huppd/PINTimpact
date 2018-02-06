#pragma once
#ifndef PIMPACT_SCALARFIELD_HPP
#define PIMPACT_SCALARFIELD_HPP


#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

#include "mpi.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Tuple.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_AbstractField.hpp"
#include "Pimpact_extern_ScalarField.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {


/// \brief important basic Vector class
/// vector for a scalar field, e.g.: pressure,
/// \ingroup Field
template<class GridType>
class ScalarField : private AbstractField<GridType > {

public:

  using GridT = GridType;

protected:

  using ST = typename GridT::Scalar;
  using OT = typename GridT::Ordinal;

  using ScalarArray = ST*;
  using State = Teuchos::Tuple<bool, 3>;

  using SW = typename GridT::SW;

  ScalarArray s_;

  const Owning owning_;

  State exchangedState_;

  const F fType_; /// <make template parameter (default:=S)

  const OT stride1_;
  const OT stride2_;

  void allocate() {
    OT n = getStorageSize();
    setStoragePtr(new ST[n]);
    std::uninitialized_fill_n(s_, n , 0.);
  }

public:

  ScalarField(const Teuchos::RCP<const GridT>& grid, const Owning owning=Owning::Y,
      const F fType=F::S):
    AbstractField<GridT>(grid),
    owning_(owning),
    exchangedState_(Teuchos::tuple(true, true, true)),
    fType_(fType),
    stride1_(grid->nLoc(0)+SW::BU(0)-SW::BL(0)+1),
    stride2_((grid()->nLoc(0)+SW::BU(0)-SW::BL(0)+1)*(
          grid->nLoc(1)+SW::BU(1)-SW::BL(1)+1)) {

      if(owning_ == Owning::Y) allocate();
    };


  /// \brief copy constructor.
  ///
  /// \note copyType is
  /// \param sF ScalarField which is copied
  /// \param copyType by default a ECopy::Deep is done but also allows to ECopy::Shallow
  ScalarField(const ScalarField& sF, const ECopy copyType=ECopy::Deep):
    AbstractField<GridT>(sF.grid()),
    owning_(sF.owning_),
    exchangedState_(sF.exchangedState_),
    fType_(sF.fType_),
    stride1_(sF.stride1_),
    stride2_(sF.stride2_) {

      if(owning_ == Owning::Y) {

        allocate();

        switch(copyType) {
          case ECopy::Shallow:
            break;
          case ECopy::Deep:
            *this = sF;
            break;
        }
      }
    };


  ~ScalarField() {
    if(owning_ == Owning::Y) delete[] s_;
  }


  Teuchos::RCP<ScalarField> clone(const ECopy copyType=ECopy::Deep) const {

    Teuchos::RCP<ScalarField> mv = Teuchos::rcp(new ScalarField(grid(), Owning::Y, this->fType_));

    switch(copyType) {
    case ECopy::Shallow:
      break;
    case ECopy::Deep:
      *mv = *this;
      break;
    }
    return mv;
  }

  /// \name Attribute methods
  /// \{

  /// \brief returns the length of Field.
  constexpr OT getLength() {

    Teuchos::RCP<const BoundaryConditionsGlobal<GridT::dimension> > bc =
      grid()->getBCGlobal();

    OT vl = 1;

    for(int dir = 0; dir<GridT::sdim; ++dir) {
      vl *= grid()->nGlo(dir) +
            ((BC::Periodic == bc->getBCL(dir))?
              -1:
              (fType_ == dir)?1:0);
    }

    return vl;
  }



  /// @}
  /// \name Update methods
  /// @{

  /// \brief Replace \c this with \f$\alpha a + \beta B\f$.
  /// \todo make checks for grids and k
  void add(const ST alpha, const ScalarField& a, const ST beta, const
            ScalarField& b, const B wb=B::Y) {

    assert(a.getType() == b.getType());
    assert(getType() == b.getType());
#ifndef NDEBUG
    for(int dir=0; dir<3; ++dir) {
      bool same_space = grid()->nLoc(dir)>a.grid()->nLoc(dir) ||
                        grid()->nLoc(dir)>b.grid()->nLoc(dir);
      assert(!same_space);
      bool consistent_space = (
                                (a.grid()->nLoc(dir)-1)%(grid()->nLoc(dir)-1))!=0 || (
                                (b.grid()->nLoc(dir)-1)%(grid()->nLoc(dir)-1))!=0 ;
      assert(!consistent_space);
    }
#endif
    Teuchos::Tuple<OT, 3> da;
    Teuchos::Tuple<OT, 3> db;

    bool with_d_yes = false;
    for(int dir=0; dir<3; ++dir) {
      da[dir] = (a.grid()->nLoc(dir)-1)/(grid()->nLoc(dir)-1);
      db[dir] = (b.grid()->nLoc(dir)-1)/(grid()->nLoc(dir)-1);
      if(1!=da[dir]) with_d_yes=true;
      if(1!=db[dir]) with_d_yes=true;
    }

    if(with_d_yes) {
      for(OT k=grid()->si(fType_, Z, wb); k<=grid()->ei(fType_, Z, wb); ++k)
        for(OT j=grid()->si(fType_, Y, wb); j<=grid()->ei(fType_, Y, wb); ++j)
          for(OT i=grid()->si(fType_, X, wb); i<=grid()->ei(fType_, X, wb); ++i)
            at(i, j, k) = alpha*a.at((i-1)*da[0]+1, (j-1)*da[1]+1, (k-1)*da[2]+1)
              + beta*b.at((i-1)*db[0]+1, (j-1)*db[1]+1, (k-1)*db[2]+1);
    }
    else {
      for(OT k=grid()->si(fType_, Z, wb); k<=grid()->ei(fType_, Z, wb); ++k)
        for(OT j=grid()->si(fType_, Y, wb); j<=grid()->ei(fType_, Y, wb); ++j)
          for(OT i=grid()->si(fType_, X, wb); i<=grid()->ei(fType_, X, wb); ++i)
            at(i, j, k) = alpha*a.at(i, j, k) + beta*b.at(i, j, k);
    }

    changed();
  }


  /// \brief Put element-wise absolute values of source vector \c y into this
  /// vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = | y_i | \quad \mbox{for } i=1, \dots, n \f]
  /// \return Reference to this object
  void abs(const ScalarField& y, const B bcYes=B::Y) {

    for(int dir=0; dir<3; ++dir)
      assert(grid()->nLoc(dir) == y.grid()->nLoc(dir));

    for(OT k=grid()->si(fType_, Z, bcYes); k<=grid()->ei(fType_, Z, bcYes); ++k)
      for(OT j=grid()->si(fType_, Y, bcYes); j<=grid()->ei(fType_, Y, bcYes); ++j)
        for(OT i=grid()->si(fType_, X, bcYes); i<=grid()->ei(fType_, X, bcYes); ++i)
          at(i, j, k) = std::fabs(y.at(i, j, k));

    changed();
  }


  /// \brief Put element-wise reciprocal of source vector \c y into this vector.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i =  \frac{1}{y_i} \quad \mbox{for } i=1, \dots, n  \f]
  /// \return Reference to this object
  void reciprocal(const ScalarField& y, const B bcYes=B::Y) {

#ifndef NDEBUG
    for(int dir=0; dir<3; ++dir) {
      bool same_space = grid()->nLoc(dir)!=y.grid()->nLoc(dir);
      assert(!same_space);
    }
#endif

    for(OT k=grid()->si(fType_, Z, bcYes); k<=grid()->ei(fType_, Z, bcYes); ++k)
      for(OT j=grid()->si(fType_, Y, bcYes); j<=grid()->ei(fType_, Y, bcYes); ++j)
        for(OT i=grid()->si(fType_, X, bcYes); i<=grid()->ei(fType_, X, bcYes); ++i)
          at(i, j, k) = Teuchos::ScalarTraits<ST>::one()/ y.at(i, j, k);

    changed();
  }


  /// \brief Scale each element of the vector with \c alpha.
  void scale(const ST alpha, const B wB=B::Y) {

    for(OT k=grid()->si(fType_, Z, wB); k<=grid()->ei(fType_, Z, wB); ++k)
      for(OT j=grid()->si(fType_, Y, wB); j<=grid()->ei(fType_, Y, wB); ++j)
        for(OT i=grid()->si(fType_, X, wB); i<=grid()->ei(fType_, X, wB); ++i)
          at(i, j, k) *= alpha;

    changed();
  }


  /// \brief Scale this vector <em>element-by-element</em> by the vector a.
  ///
  /// Here x represents this vector, and we update it as
  /// \f[ x_i = x_i \cdot y_i \quad \mbox{for } i=1, \dots, n \f]
  void scale(const ScalarField& y, const B bcYes=B::Y) {

#ifndef NDEBUG
    for(int dir=0; dir<3; ++dir) {
      bool same_space = grid()->nLoc(dir)!=y.grid()->nLoc(dir);
      assert(!same_space);
    }
#endif

    for(OT k=grid()->si(fType_, Z, bcYes); k<=grid()->ei(fType_, Z, bcYes); ++k)
      for(OT j=grid()->si(fType_, Y, bcYes); j<=grid()->ei(fType_, Y, bcYes); ++j)
        for(OT i=grid()->si(fType_, X, bcYes); i<=grid()->ei(fType_, X, bcYes); ++i)
          at(i, j, k) *= y.at(i, j, k);
    changed();
  }

  /// @}
  /// \name Norm method(reductions)
  /// @{

  /// \brief Compute a local scalar \c b, which is the dot-product of \c y and \c this, i.e.\f$b = y^H this\f$.
  constexpr ST dotLoc(const ScalarField& y, const B bcYes=B::Y) {

#ifndef NDEBUG
    for(int dir=0; dir<3; ++dir) {
      bool same_space = grid()->nLoc(dir)!=y.grid()->nLoc(dir);
      assert(!same_space);
    }
#endif

    ST b = Teuchos::ScalarTraits<ST>::zero();
    //auto coord = grid()->getCoordinatesLocal();

    for(OT k=grid()->si(fType_, Z, bcYes); k<=grid()->ei(fType_, Z, bcYes); ++k)
      for(OT j=grid()->si(fType_, Y, bcYes); j<=grid()->ei(fType_, Y, bcYes); ++j)
        for(OT i=grid()->si(fType_, X, bcYes); i<=grid()->ei(fType_, X, bcYes); ++i) {
          //ST volume = coord->dx(fType_, X, i) * coord->dx(fType_, Y, j) * coord->dx(fType_, Z, k);
          b += /*volume**/at(i, j, k)*y.at(i, j, k);
        }

    return b;
  }

  /// \brief Compute/reduces a scalar \c b, which is the dot-product of \c y
  /// and \c this, i.e.\f$b = y^H this\f$.
  constexpr ST dot(const ScalarField& y, const B bcYes=B::Y) {
    return this->reduce(comm(), dotLoc(y, bcYes));
  }

  constexpr ST normLoc1(const B bcYes=B::Y) {

    ST normvec = Teuchos::ScalarTraits<ST>::zero();

    for(OT k=grid()->si(fType_, Z, bcYes); k<=grid()->ei(fType_, Z, bcYes); ++k)
      for(OT j=grid()->si(fType_, Y, bcYes); j<=grid()->ei(fType_, Y, bcYes); ++j)
        for(OT i=grid()->si(fType_, X, bcYes); i<=grid()->ei(fType_, X, bcYes); ++i)
          normvec += std::fabs(at(i, j, k));

    return normvec;
  }


  constexpr ST normLoc2(const B bcYes=B::Y) {

    //return normLocL2(bcYes);
    ST normvec = Teuchos::ScalarTraits<ST>::zero();

    for(OT k=grid()->si(fType_, Z, bcYes); k<=grid()->ei(fType_, Z, bcYes); ++k)
      for(OT j=grid()->si(fType_, Y, bcYes); j<=grid()->ei(fType_, Y, bcYes); ++j)
        for(OT i=grid()->si(fType_, X, bcYes); i<=grid()->ei(fType_, X, bcYes); ++i)
          normvec += std::pow(at(i, j, k), 2);

    return normvec;
  }


  constexpr ST normLocInf(const B bcYes=B::Y) {

    ST normvec = Teuchos::ScalarTraits<ST>::zero();

    for(OT k=grid()->si(fType_, Z, bcYes); k<=grid()->ei(fType_, Z, bcYes); ++k)
      for(OT j=grid()->si(fType_, Y, bcYes); j<=grid()->ei(fType_, Y, bcYes); ++j)
        for(OT i=grid()->si(fType_, X, bcYes); i<=grid()->ei(fType_, X, bcYes); ++i)
          normvec = std::fmax(std::fabs(at(i, j, k)), normvec);

    return normvec;
  }


  constexpr ST normLocL2(const B bcYes=B::Y) {

    ST normvec = Teuchos::ScalarTraits<ST>::zero();

    auto coord = grid()->getCoordinatesLocal();

    for(OT k=grid()->si(fType_, Z, bcYes); k<=grid()->ei(fType_, Z, bcYes); ++k)
      for(OT j=grid()->si(fType_, Y, bcYes); j<=grid()->ei(fType_, Y, bcYes); ++j)
        for(OT i=grid()->si(fType_, X, bcYes); i<=grid()->ei(fType_, X, bcYes); ++i) {
          ST volume = coord->dx(fType_, X, i) * coord->dx(fType_, Y, j) * coord->dx(fType_, Z, k);
          normvec += volume*std::pow(at(i, j, k), 2);
        }

    return normvec;
  }


  constexpr ST normLoc(const ENorm type=ENorm::Two, const B bcYes=B::Y) {

    return (ENorm::One == type)?
      normLoc1(bcYes):
      (ENorm::Two == type)?
        normLoc2(bcYes): (ENorm::Inf == type)?
          normLocInf(bcYes): normLocL2(bcYes);
  }


  /// \brief compute the norm
  /// \return by default holds the value of \f$||this||_2\f$, or in the specified norm.
  constexpr ST norm(const ENorm type=ENorm::Two, const B bcYes=B::Y) {

    ST normvec = this->reduce(comm(), normLoc(type, bcYes),
        (ENorm::Inf == type)?MPI_MAX:MPI_SUM);

    normvec = (ENorm::Two == type||ENorm::L2 == type) ?
      std::sqrt(normvec) :
      normvec;

    return normvec;
  }


  /// \brief Weighted 2-Norm.
  ///
  /// \warning untested
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr ST normLoc(const ScalarField& weights, const B bcYes=B::Y) {

    for(int dir=0; dir<3; ++dir)
      assert(grid()->nLoc(dir) == weights.grid()->nLoc(dir));

    ST normvec = Teuchos::ScalarTraits<ST>::zero();

    for(OT k=grid()->si(fType_, Z, bcYes); k<=grid()->ei(fType_, Z, bcYes); ++k)
      for(OT j=grid()->si(fType_, Y, bcYes); j<=grid()->ei(fType_, Y, bcYes); ++j)
        for(OT i=grid()->si(fType_, X, bcYes); i<=grid()->ei(fType_, X, bcYes); ++i)
          normvec += at(i, j, k)*at(i, j, k)*weights.at(i, j, k)*weights.at(i, j, k);

    return normvec;
  }


  /// \brief Weighted 2-Norm.
  ///
  /// \warning untested
  /// Here x represents this vector, and we compute its weighted norm as follows:
  /// \f[ \|x\|_w = \sqrt{\sum_{i=1}^{n} w_i \; x_i^2} \f]
  /// \return \f$ \|x\|_w \f$
  constexpr ST norm(const ScalarField& weights, const B bcYes=B::Y) {
    return std::sqrt(this->reduce(comm(), normLoc(weights, bcYes)));
  }


  //\}
  /// \name Initialization methods
  //\{

  /// \brief *this := a
  ///
  /// Assign (deep copy) \c a into \c this.
  /// total deep, boundaries and everything.
  /// \note the \c StencilWidths is not take care of assuming every field is generated with one
  ScalarField& operator=(const ScalarField& a) {

    assert(getType() == a.getType());
    assert(getStorageSize() == a.getStorageSize());

    std::copy_n(a.s_, getStorageSize(), s_);

    for(int dir=0; dir<GridT::sdim; ++dir)
      exchangedState_[dir] = a.exchangedState_[dir];
    return *this;
  }

  /// \brief Replace the vectors with a random vectors.
  /// Depending on Fortrans \c Random_number implementation, with always same seed => not save, if good randomness is required
  void random(bool useSeed = false, const B bcYes=B::Y , int seed = 1) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.5, 0.5);

    for(OT k=grid()->si(fType_, Z, bcYes); k<=grid()->ei(fType_, Z, bcYes); ++k)
      for(OT j=grid()->si(fType_, Y, bcYes); j<=grid()->ei(fType_, Y, bcYes); ++j)
        for(OT i=grid()->si(fType_, X, bcYes); i<=grid()->ei(fType_, X, bcYes); ++i)
          at(i, j, k) = dis(gen);

    if(!grid()->getProcGrid()->participating())
      for(OT k=grid()->si(fType_, Z, B::Y); k<=grid()->ei(fType_, Z, B::Y); ++k)
        for(OT j=grid()->si(fType_, Y, B::Y); j<=grid()->ei(fType_, Y, B::Y); ++j)
          for(OT i=grid()->si(fType_, X, B::Y); i<=grid()->ei(fType_, X, B::Y); ++i)
            at(i, j, k) = Teuchos::ScalarTraits<ST>::zero();
    changed();
  }


  /// \brief Replace each element of the vector  with \c alpha.
  /// \param alpha init value
  /// \param bcYes also initializing the boundary values
  void init(const ST alpha = Teuchos::ScalarTraits<ST>::zero(), const B bcYes=B::Y) {

    if(B::Y == bcYes) {
      std::fill_n(s_, getStorageSize(), alpha);
      exchangedState_[X] = true;
      exchangedState_[Y] = true;
      exchangedState_[Z] = true;
    } else {
      for(OT k=grid()->si(fType_, Z, bcYes); k<=grid()->ei(fType_, Z, bcYes); ++k)
        for(OT j=grid()->si(fType_, Y, bcYes); j<=grid()->ei(fType_, Y, bcYes); ++j)
          for(OT i=grid()->si(fType_, X, bcYes); i<=grid()->ei(fType_, X, bcYes); ++i)
            at(i, j, k) = alpha;

      if(!grid()->getProcGrid()->participating())
        for(OT k=grid()->si(fType_, Z, B::Y); k<=grid()->ei(fType_, Z, B::Y); ++k)
          for(OT j=grid()->si(fType_, Y, B::Y); j<=grid()->ei(fType_, Y, B::Y); ++j)
            for(OT i=grid()->si(fType_, X, B::Y); i<=grid()->ei(fType_, X, B::Y); ++i)
              at(i, j, k) = Teuchos::ScalarTraits<ST>::zero();
      changed();
    }
  }

protected:

  /// \brief helper function getting \c EScalarField for switch statement from name
  ///
  /// \param name input name
  /// \return according int number
  EScalarField string2enum(const std::string& name) {

    std::string lcName = name;
    std::transform(lcName.begin(), lcName.end(), lcName.begin(), ::tolower);

    if("constant" == lcName) return ConstField;
    else if("grad in x" == lcName) return Grad2D_inX;
    else if("grad in y" == lcName) return Grad2D_inY;
    else if("grad in z" == lcName) return Grad2D_inZ;
    else if("poiseuille" == lcName) return Poiseuille2D_inX;
    else if("poiseuille in x" == lcName) return Poiseuille2D_inX;
    else if("poiseuille in y" == lcName) return Poiseuille2D_inY;
    else if("poiseuille in z" == lcName) return Poiseuille2D_inZ;
    else if("point" == lcName) return FPoint;
    else {
#ifndef NDEBUG
      const bool& Flow_Type_not_known = false;
      assert(Flow_Type_not_known);
#endif
    }
    return ConstField; // just to please the compiler
  }

public:

  /// \brief initializes field from a lambda function
  ///
  /// \tparam Functor need a to have an Operator (x, y, z)->u
  /// \param func  note that x, y, z, should be defined between 0 and one the scaling happens here
  template<typename Functor>
  void initFromFunction(Functor&& func, const Add add=Add::N) {

    Teuchos::RCP<const CoordinatesLocal<ST, OT, GridT::dimension, GridT::dimNC> > coord =
      grid()->getCoordinatesLocal();
    Teuchos::RCP<const DomainSize<ST, GridT::sdim> > domain = grid()->getDomainSize();

    const B bY = B::Y;

    for(OT k=grid()->si(fType_, Z, bY); k<=grid()->ei(fType_, Z, bY); ++k)
      for(OT j=grid()->si(fType_, Y, bY); j<=grid()->ei(fType_, Y, bY); ++j)
        for(OT i=grid()->si(fType_, X, bY); i<=grid()->ei(fType_, X, bY); ++i) {
          ST val = func(
                (coord->getX(fType_, X, i)-domain->getOrigin(X))/domain->getSize(X),
                (coord->getX(fType_, Y, j)-domain->getOrigin(Y))/domain->getSize(Y),
                (coord->getX(fType_, Z, k)-domain->getOrigin(Z))/domain->getSize(Z));
          if(Add::Y == add)
            at(i, j, k) += val;
          else
            at(i, j, k) = val;
        }
    changed();
  }


  ///  \brief initializes including boundaries to zero
  void initField(Teuchos::ParameterList& para, const Add add=Add::N) {

    EScalarField type =
      string2enum(para.get<std::string>("Type", "constant"));

    switch(type) {
      case ConstField : {
        if(Add::N == add) init();
        break;
      }
      case Grad2D_inX : {
        ST a = para.get<ST>("dx", Teuchos::ScalarTraits<ST>::one());
        initFromFunction(
            [&a] (ST x, ST y, ST z)->ST { return a*(x-0.5); },
            add);
        break;
      }
      case Grad2D_inY : {
        ST a = para.get<ST>("dy", Teuchos::ScalarTraits<ST>::one());
        initFromFunction(
            [&a] (ST x, ST y, ST z)->ST { return a*(y-0.5); },
            add);
        break;
      }
      case Grad2D_inZ : {
        ST a = para.get<ST>("dz", Teuchos::ScalarTraits<ST>::one());
        initFromFunction(
            [&a] (ST x, ST y, ST z)->ST { return a*(z-0.5); },
            add);
        break;
      }
      case Poiseuille2D_inX : {
        initFromFunction(
            [] (ST x, ST y, ST z)->ST { return 4.*x*(1.-x); },
            add);
        break;
      }
      case Poiseuille2D_inY : {
        initFromFunction(
            [] (ST x, ST y, ST z)->ST { return 4.*y*(1.-y); },
            add);
        break;
      }
      case Poiseuille2D_inZ : {
        initFromFunction(
            [] (ST x, ST y, ST z)->ST { return 4.*z*(1.-z); },
            add);
        break;
      }
      case FPoint : {
        ST xc[3] = {
          para.get<ST>("c_x", Teuchos::ScalarTraits<ST>::one()),
          para.get<ST>("c_y", grid()->getDomainSize()->getSize(Y)/2.),
          para.get<ST>("c_z", grid()->getDomainSize()->getSize(Z)/2.)
        };
        ST amp = para.get<ST>("amp", Teuchos::ScalarTraits<ST>::one());
        ST sig[3] = {
          para.get<ST>("sig_x", 0.2),
          para.get<ST>("sig_y", 0.2),
          para.get<ST>("sig_z", 0.2)
        };

        Teuchos::RCP<const DomainSize<ST, GridT::sdim> > domain = grid()->getDomainSize();

        initFromFunction(
            [&xc, &amp, &sig, &domain] (ST x_, ST y_, ST z_)->ST {

            ST x = x_*domain->getSize(X) + domain->getOrigin(X);
            ST y = y_*domain->getSize(Y) + domain->getOrigin(Y);
            ST z = z_*domain->getSize(Z) + domain->getOrigin(Z);
            return amp*std::exp(
              -std::pow((x-xc[0])/sig[0], 2)
              -std::pow((y-xc[1])/sig[1], 2)
              -std::pow((z-xc[2])/sig[2], 2)); },
            add);
        break;
      }
    }

    if(!grid()->getProcGrid()->participating()) // not sure why?
      init(Teuchos::ScalarTraits<ST>::zero());

    changed();
  }


  /// \brief initializes VectorField with the initial field defined in Fortran
  /// \deprecated
  void initField(const EScalarField fieldType, const ST alpha=Teuchos::ScalarTraits<ST>::zero()) {

    switch(fieldType) {
      case ConstField : {
        init(alpha);
        break;
      }
      case Grad2D_inX : {
        ST a = (std::fabs(alpha)<Teuchos::ScalarTraits<ST>::eps())?Teuchos::ScalarTraits<ST>::one():alpha;
        initFromFunction([&a] (ST x, ST y, ST z)->ST { return a*(x-0.5); });
        break;
      }
      case Grad2D_inY : {
        ST a = (std::fabs(alpha)<Teuchos::ScalarTraits<ST>::eps())?Teuchos::ScalarTraits<ST>::one():alpha;
        initFromFunction([&a] (ST x, ST y, ST z)->ST { return a*(y-0.5); });
        break;
      }
      case Grad2D_inZ : {
        ST a = (std::fabs(alpha)<Teuchos::ScalarTraits<ST>::eps())?Teuchos::ScalarTraits<ST>::one():alpha;
        initFromFunction([&a] (ST x, ST y, ST z)->ST { return a*(z-0.5); });
        break;
      }
      case Poiseuille2D_inX : {
        initFromFunction([] (ST x, ST y, ST z)->ST { return 4.*x*(1.-x); });
        break;
      }
      case Poiseuille2D_inY : {
        initFromFunction([] (ST x, ST y, ST z)->ST { return 4.*y*(1.-y); });
        break;
      }
      case Poiseuille2D_inZ : {
        initFromFunction([] (ST x, ST y, ST z)->ST { return 4.*z*(1.-z); });
        break;
      }
      case FPoint : {
        ST xc[3] = { 0.5, 0.5, 0.5 };
        ST amp = alpha;
        ST sig[3] = { 0.2, 0.2, 0.2 };
        initFromFunction([&xc, &amp, &sig] (ST x, ST y, ST z)->ST {
            return amp*std::exp(
              -std::pow((x-xc[0])/sig[0], 2)
              -std::pow((x-xc[1])/sig[1], 2)
              -std::pow((x-xc[2])/sig[2], 2)); }
            );
        break;
      }
    }

    if(!grid()->getProcGrid()->participating())
      init(Teuchos::ScalarTraits<ST>::zero());
    changed();
  }


  /// \brief for Dirichlet BC extrapolate the velocity points outside the domain such that the
  ///  interpolated value on the boundary is zero
  ///
  /// \test Neumann BC
  /// \param trans transposed
  void extrapolateBC(const Belos::ETrans trans=Belos::NOTRANS) {

    switch(trans) {
    case(Belos::NOTRANS): {
      switch(fType_) {
      case(F::U):  {
        using StencD = Stencil<ST, OT, 0, SW::DL(0), SW::DU(0) >;

        StencD c_(grid()->nLoc(X));

        if(grid()->bcl(X) == BC::Neumann || grid()->bcu(X) == BC::Neumann)
          FD_getDiffCoeff(
            1,
            grid()->nLoc(X),
            grid()->bl(X),
            grid()->bu(X),
            grid()->dl(X),
            grid()->du(X),
            grid()->getBCLocal()->getBCL(X),
            grid()->getBCLocal()->getBCU(X),
            grid()->getShift(X),
            3,
            1,
            1,
            0,
            false, // mapping
            grid()->getStencilWidths()->getDimNcbD(X),
            grid()->getStencilWidths()->getNcbD(X),
            grid()->getCoordinatesLocal()->getX(F::U, X),
            grid()->getCoordinatesLocal()->getX(F::S, X),
            c_.get());

        if(0 <grid()->bcl(X)) {
          OT i = grid()->si(fType_, X, B::Y);
          for(OT k=grid()->si(fType_, Z, B::Y); k<=grid()->ei(fType_, Z, B::Y); ++k)
            for(OT j=grid()->si(fType_, Y, B::Y); j<=grid()->ei(fType_, Y, B::Y); ++j) {
              at(i, j, k) = 0.;
              if(BC::Dirichlet == grid()->bcl(X)) {
                for(OT ii=0; ii<=SW::DU(X); ++ii)
                  at(i, j, k) -= at(1+ii, j, k)*grid()->getInterpolateV2S()->getC(X, 1, ii)/grid()->getInterpolateV2S()->getC(X, 1, -1);
              } else if(BC::Neumann == grid()->bcl(X)) {
                for(OT ii=0; ii<=SW::DU(X); ++ii)
                  at(i, j, k) -= at(1+ii, j, k)*c_(1, ii)/c_(1, -1);
              }
            }
        }
        if(0 <grid()->bcu(X)) {

          OT i = grid()->ei(fType_, X, B::Y);
          for(OT k=grid()->si(fType_, Z, B::Y); k<=grid()->ei(fType_, Z, B::Y); ++k)
            for(OT j=grid()->si(fType_, Y, B::Y); j<=grid()->ei(fType_, Y, B::Y); ++j) {
              at(i, j, k) = 0.;
              if(BC::Dirichlet == grid()->bcu(X)) {
                for(OT ii=SW::DL(X); ii<=-1; ++ii)
                  at(i, j, k) -= grid()->getInterpolateV2S()->getC(X, i, ii)*at(i+ii, j, k)/grid()->getInterpolateV2S()->getC(X, i, 0);
              } else if(BC::Neumann == grid()->bcu(X)) {
                for(OT ii=SW::DL(X); ii<=-1; ++ii)
                  at(i, j, k) -= c_(i, ii)*at(i+ii, j, k)/c_(i, 0);
              }
            }
        }
        break;
      }
      case(F::V) : {
        using StencD = Stencil<ST, OT, 0, SW::DL(0), SW::DU(0) >;

        StencD c_(grid()->nLoc(Y));

        if(grid()->bcl(Y) == BC::Neumann || grid()->bcu(Y) == BC::Neumann)
          FD_getDiffCoeff(
            1,
            grid()->nLoc(Y),
            grid()->bl(Y),
            grid()->bu(Y),
            grid()->dl(Y),
            grid()->du(Y),
            grid()->getBCLocal()->getBCL(Y),
            grid()->getBCLocal()->getBCU(Y),
            grid()->getShift(Y),
            3,
            1,
            1,
            0,
            false, // mapping
            grid()->getStencilWidths()->getDimNcbD(Y),
            grid()->getStencilWidths()->getNcbD(Y),
            grid()->getCoordinatesLocal()->getX(F::V, Y),
            grid()->getCoordinatesLocal()->getX(F::S, Y),
            c_.get());

        if(0 <grid()->bcl(Y)) {
          OT j = grid()->si(fType_, Y, B::Y);
          for(OT k=grid()->si(fType_, Z, B::Y); k<=grid()->ei(fType_, Z, B::Y); ++k)
            for(OT i=grid()->si(fType_, X, B::Y); i<=grid()->ei(fType_, X, B::Y); ++i) {
              at(i, j, k) = 0.;
              if(BC::Dirichlet == grid()->bcl(Y)) {
                for(OT jj=0; jj<=SW::DU(Y); ++jj)
                  at(i, j, k) -= at(i, 1+jj, k)*grid()->getInterpolateV2S()->getC(Y, 1, jj)/grid()->getInterpolateV2S()->getC(Y, 1, -1);
              } else if(BC::Neumann == grid()->bcl(Y)) {
                for(OT jj=0; jj<=SW::DU(Y); ++jj)
                  at(i, j, k) -= at(i, 1+jj, k)*c_(1, jj)/c_(1, -1);
              }
            }
        }
        if(0 <grid()->bcu(Y)) {
          OT j = grid()->ei(fType_, Y, B::Y);
          for(OT k=grid()->si(fType_, Z, B::Y); k<=grid()->ei(fType_, Z, B::Y); ++k)
            for(OT i=grid()->si(fType_, X, B::Y); i<=grid()->ei(fType_, X, B::Y); ++i) {
              at(i, j, k) = 0.;
              if(BC::Dirichlet == grid()->bcu(Y)) {
                for(OT jj=SW::DL(Y); jj<=-1; ++jj)
                  at(i, j, k) -= grid()->getInterpolateV2S()->getC(Y, j, jj)*at(i, j+jj, k)/grid()->getInterpolateV2S()->getC(Y, j, 0);
              } else if(BC::Neumann == grid()->bcu(Y)) {
                for(OT jj=SW::DL(Y); jj<=-1; ++jj)
                  at(i, j, k) -= c_(j, jj)*at(i, j+jj, k)/c_(j, 0);
              }
            }
        }
        break;
      }
      case(F::W) : {
        using StencD = Stencil<ST, OT, 0, SW::DL(0), SW::DU(0) >;

        StencD c_(grid()->nLoc(Z));

        if(grid()->bcl(Z) == BC::Neumann || grid()->bcu(Z) == BC::Neumann)
          FD_getDiffCoeff(
            1,
            grid()->nLoc(Z),
            grid()->bl(Z),
            grid()->bu(Z),
            grid()->dl(Z),
            grid()->du(Z),
            grid()->getBCLocal()->getBCL(Z),
            grid()->getBCLocal()->getBCU(Z),
            grid()->getShift(Z),
            3,
            1,
            1,
            0,
            false, // mapping
            grid()->getStencilWidths()->getDimNcbD(Z),
            grid()->getStencilWidths()->getNcbD(Z),
            grid()->getCoordinatesLocal()->getX(F::W, Z),
            grid()->getCoordinatesLocal()->getX(F::S, Z),
            c_.get());

        if(grid()->bcl(Z) > 0) {
          OT k = grid()->si(fType_, Z, B::Y);
          for(OT j=grid()->si(fType_, Y, B::Y); j<=grid()->ei(fType_, Y, B::Y); ++j)
            for(OT i=grid()->si(fType_, X, B::Y); i<=grid()->ei(fType_, X, B::Y); ++i) {
              at(i, j, k) = 0.;
              if(BC::Dirichlet == grid()->bcl(Z)) {
                for(OT kk=0; kk<=SW::DU(Z); ++kk)
                  at(i, j, k) -= grid()->getInterpolateV2S()->getC(Z, 1, kk)*at(i, j, 1+kk)/grid()->getInterpolateV2S()->getC(Z, 1, -1);
              } else if(BC::Neumann == grid()->bcl(Z)) {
                for(OT kk=0; kk<=SW::DU(Z); ++kk)
                  at(i, j, k) -= c_(1, kk)*at(i, j, 1+kk)/c_(1, -1);
              }
            }
        }
        if(grid()->bcu(Z) > 0) {
          OT k = grid()->ei(fType_, Z, B::Y);
          for(OT j=grid()->si(fType_, Y, B::Y); j<=grid()->ei(fType_, Y, B::Y); ++j)
            for(OT i=grid()->si(fType_, X, B::Y); i<=grid()->ei(fType_, X, B::Y); ++i) {
              at(i, j, k) = 0.;
              if(BC::Dirichlet == grid()->bcu(Z)) {
                for(OT kk=SW::DL(Z); kk<=-1; ++kk)
                  at(i, j, k) -= grid()->getInterpolateV2S()->getC(Z, k, kk)*at(i, j, k+kk)/grid()->getInterpolateV2S()->getC(Z, k, 0);
              } else if(BC::Neumann == grid()->bcu(Z)) {
                for(OT kk=SW::DL(Z); kk<=-1; ++kk)
                  at(i, j, k) -= c_(k, kk)*at(i, j, k+kk)/c_(k, 0);
              }
            }
        }
        break;
      }
      case(F::S) :
        break;
        //case(F::end) : break;
      }
      break;
    }
    case Belos::TRANS : {

      switch(fType_) {
      case(F::U) : {
        if(grid()->bcl(X) > 0) {
          OT i = grid()->si(fType_, X, B::Y);
          for(OT k=grid()->si(fType_, Z, B::Y); k<=grid()->ei(fType_, Z, B::Y); ++k)
            for(OT j=grid()->si(fType_, Y, B::Y); j<=grid()->ei(fType_, Y, B::Y); ++j) {
              for(OT ii=0; ii<=SW::DU(X); ++ii)
                at(i+ii+1, j, k) -= at(i, j, k)*grid()->getInterpolateV2S()->getC(X, 1, ii)/grid()->getInterpolateV2S()->getC(X, 1, -1);
              at(i, j, k) = 0.;
            }
        }
        if(grid()->bcu(X) > 0) {
          OT i = grid()->ei(fType_, X, B::Y);
          for(OT k=grid()->si(fType_, Z, B::Y); k<=grid()->ei(fType_, Z, B::Y); ++k)
            for(OT j=grid()->si(fType_, Y, B::Y); j<=grid()->ei(fType_, Y, B::Y); ++j) {
              for(OT ii=SW::DL(X); ii<=-1; ++ii)
                at(i+ii, j, k) -= grid()->getInterpolateV2S()->getC(X, i, ii)*at(i, j, k)/grid()->getInterpolateV2S()->getC(X, i, 0);
              at(i, j, k) = 0.;
            }
        }
        break;
      }
      case(F::V) : {
        if(grid()->bcl(Y) > 0) {
          OT j = grid()->si(fType_, Y, B::Y);
          for(OT k=grid()->si(fType_, Z, B::Y); k<=grid()->ei(fType_, Z, B::Y); ++k)
            for(OT i=grid()->si(fType_, X, B::Y); i<=grid()->ei(fType_, X, B::Y); ++i) {
              for(OT jj=0; jj<=SW::DU(Y); ++jj)
                at(i, j+jj+1, k) -= at(i, j, k)*grid()->getInterpolateV2S()->getC(Y, 1, jj)/grid()->getInterpolateV2S()->getC(Y, 1, -1);
              at(i, j, k) = 0.;
            }
        }
        if(grid()->bcu(Y) > 0) {
          OT j = grid()->ei(fType_, Y, B::Y);
          for(OT k=grid()->si(fType_, Z, B::Y); k<=grid()->ei(fType_, Z, B::Y); ++k)
            for(OT i=grid()->si(fType_, X, B::Y); i<=grid()->ei(fType_, X, B::Y); ++i) {
              for(OT jj=SW::DL(Y); jj<=-1; ++jj)
                at(i, j+jj, k) -= grid()->getInterpolateV2S()->getC(Y, j, jj)*at(i, j, k)/grid()->getInterpolateV2S()->getC(Y, j, 0);
              at(i, j, k) = 0.;
            }
        }
        break;
      }
      case(F::W) : {
        if(grid()->bcl(Z) > 0) {
          OT k = grid()->si(fType_, Z, B::Y);
          for(OT j=grid()->si(fType_, Y, B::Y); j<=grid()->ei(fType_, Y, B::Y); ++j)
            for(OT i=grid()->si(fType_, X, B::Y); i<=grid()->ei(fType_, X, B::Y); ++i) {
              at(i, j, k) /= grid()->getInterpolateV2S()->getC(Z, 1, -1);
              for(OT kk=0; kk<=SW::DU(Z); ++kk)
                at(i, j, k+kk+1) -= grid()->getInterpolateV2S()->getC(Z, 1, kk)*at(i, j, k);
              at(i, j, k) = 0.;
            }
        }
        if(grid()->bcu(Z) > 0) {
          OT k = grid()->ei(fType_, Z, B::Y);
          for(OT j=grid()->si(fType_, Y, B::Y); j<=grid()->ei(fType_, Y, B::Y); ++j)
            for(OT i=grid()->si(fType_, X, B::Y); i<=grid()->ei(fType_, X, B::Y); ++i) {
              at(i, j, k) /= grid()->getInterpolateV2S()->getC(Z, k, 0);
              for(OT kk=SW::DL(Z); kk<=-1; ++kk)
                at(i, j, k+kk) -= grid()->getInterpolateV2S()->getC(Z, k, kk)*at(i, j, k);
              at(i, j, k) = 0.;
            }
        }
        break;
      }
      case(F::S) :
        break;
      }
    }
    case Belos::CONJTRANS :
      break;
    }
  }


  /// \brief levels field if scalar field
  void level() const {

    if(F::S == fType_) {

      ST pre0 = Teuchos::ScalarTraits<ST>::zero();

      for(OT k=grid()->si(fType_, Z); k<=grid()->ei(fType_, Z); ++k)
        for(OT j=grid()->si(fType_, Y); j<=grid()->ei(fType_, Y); ++j)
          for(OT i=grid()->si(fType_, X); i<=grid()->ei(fType_, X); ++i)
            pre0 += at(i, j, k);

      pre0 = this->reduce(grid()->comm(), pre0);
      pre0 /= static_cast<ST>(getLength());

      for(OT k=grid()->si(fType_, Z); k<=grid()->ei(fType_, Z); ++k)
        for(OT j=grid()->si(fType_, Y); j<=grid()->ei(fType_, Y); ++j)
          for(OT i=grid()->si(fType_, X); i<=grid()->ei(fType_, X); ++i)
            const_cast<ScalarField*>(this)->at(i, j, k) -= pre0;
    }
  }

  /// \}

  /// Print the vector.  To be used for debugging only.
  void print(std::ostream& out=std::cout)  const {

    out << "--- FieldType: " << fType_ << "--- \n";
    out << "--- StorageSize: " << getStorageSize() << "---\n";
    out << "--- owning: " << owning_ << "---\n";
    out << "--- exchangedState: " << exchangedState_ << "--\n\n";
    out << "i, \tj, \tk, \tphi(i, j, k)\n";

    Teuchos::Tuple<OT, 3> cw;
    for(int i=0; i<3; ++i)
      cw[i] = grid()->nLoc(i) + SW::BU(i) - SW::BL(i) + 1;

    for(OT k=grid()->si(fType_, Z, B::Y); k<=grid()->ei(fType_, Z, B::Y); ++k)
      for(OT j=grid()->si(fType_, Y, B::Y); j<=grid()->ei(fType_, Y, B::Y); ++j)
        for(OT i=grid()->si(fType_, X, B::Y); i<=grid()->ei(fType_, X, B::Y); ++i)
          out << i << "\t" << j << "\t" << k << "\t" << at(i, j, k) << "\n";
  }


  /// Write the ScalarField to an hdf5 file, the velocities are interpolated to the pressure points
  void write(const int count=0 , const bool restart=false) const {

    if(0 == grid()->rankS())
      switch(fType_) {
        case F::U:
          std::cout << "writing velocity field x(" << count << ") ...\n";
          break;
        case F::V:
          std::cout << "writing velocity field y(" << count << ") ...\n";
          break;
        case F::W:
          std::cout << "writing velocity field z(" << count << ") ...\n";
          break;
        case F::S:
          std::cout << "writing pressure field  (" << count << ") ...\n";
          Teuchos::Tuple<OT, 3> N;
          for(int i=0; i<3; ++i) {
            N[i] = grid()->nGlo(i);
            if(grid()->getBCGlobal()->getBCL(i) == Pimpact::BC::Periodic)
              N[i] = N[i]-1;
          }
          if(!restart) {
            std::ofstream xfile;
            std::ostringstream ss;
            ss << std::setw(5) << std::setfill('0') << count;
            std::string fname = "pre_"+ss.str();
            xfile.open(fname+".xmf", std::ofstream::out);
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
            xfile << "\t\t\t\t\t" << fname << ".h5:/VectorX\n";
            xfile << "\t\t\t\t</DataItem>\n";
            xfile << "\t\t\t\t<DataItem ItemType=\"Uniform\"\n";
            xfile << "\t\t\t\t\tDimensions=\""<< N[1] << "\"\n";
            xfile << "\t\t\t\t\tNumberType=\"Float\"\n";
            xfile << "\t\t\t\t\tPrecision=\"8\"\n";
            xfile << "\t\t\t\t\tFormat=\"HDF\">\n";
            xfile << "\t\t\t\t\t" << fname << ".h5:/VectorY\n";
            xfile << "\t\t\t\t</DataItem>\n";
            xfile << "\t\t\t\t<DataItem ItemType=\"Uniform\"\n";
            xfile << "\t\t\t\t\tDimensions=\""<< N[2] << "\"\n";
            xfile << "\t\t\t\t\tNumberType=\"Float\"\n";
            xfile << "\t\t\t\t\tPrecision=\"8\"\n";
            xfile << "\t\t\t\t\tFormat=\"HDF\">\n";
            xfile << "\t\t\t\t\t" << fname << ".h5:/VectorZ\n";
            xfile << "\t\t\t\t</DataItem>\n";
            xfile << "\t\t\t</Geometry>\n";
            xfile << "\t\t\t<Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Node\">\n";
            xfile << "\t\t\t\t<DataItem Dimensions=\""<< N[2] << " " << N[1] << " " << N[0] << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
            xfile << "\t\t\t\t\t" << fname << ".h5:/pre\n";
            xfile << "\t\t\t\t</DataItem>\n";
            xfile << "\t\t\t</Attribute>\n";
            xfile << "\t\t</Grid>\n";
            xfile << "\t</Domain>\n";
            xfile << "</Xdmf>\n";
            xfile.close();
          }
          break;
    }

    if(!restart) {
      Teuchos::RCP<ScalarField<GridT> > temp;

      if(F::S != fType_) {
        temp = Teuchos::rcp(new ScalarField<GridT>(grid(), Owning::Y, F::S));
        grid()->getInterpolateV2S()->apply(*this, *temp);
      }

      if(2 == GridT::sdim) {

        write_hdf5_2D(
          grid()->rankS(),
          MPI_Comm_c2f(grid()->comm()),
          grid()->nGlo(),
          grid()->getBCGlobal()->getBCL(),
          grid()->getBCGlobal()->getBCU(),
          grid()->nLoc(),
          grid()->bl(),
          grid()->bu(),
          grid()->sInd(F::S),
          grid()->eInd(F::S),
          grid()->getStencilWidths()->getLS(),
          grid()->np(),
          grid()->ib(),
          grid()->getShift(),
          (int)fType_,
          count,
          (F::S == fType_)?9:10,
          (F::S == fType_)?s_:temp->s_,
          grid()->getCoordinatesGlobal()->getX(F::S, 0),
          grid()->getCoordinatesGlobal()->getX(F::S, 1),
          grid()->getDomainSize()->getRe(),
          grid()->getDomainSize()->getAlpha2());
      } else if(3 == GridT::sdim) {

        int stride[3] = {1, 1, 1};

        write_hdf_3D(
            false,
            grid()->rankS(),
            MPI_Comm_c2f(grid()->comm()),
            grid()->nGlo(),
            grid()->getBCGlobal()->getBCL(),
            grid()->getBCGlobal()->getBCU(),
            grid()->nLoc(),
            grid()->bl(),
            grid()->bu(),
            grid()->sInd(F::S),
            grid()->eInd(F::S),
            grid()->getStencilWidths()->getLS(),
            grid()->np(),
            grid()->ib(),
            grid()->getShift(),
            (int)fType_+1,
            (int)F::S+1,
            count,
            (F::S == fType_)?9:10,
            stride,
            (F::S == fType_)?s_:temp->s_,
            grid()->getCoordinatesGlobal()->getX(F::S, 0),
            grid()->getCoordinatesGlobal()->getX(F::S, 1),
            grid()->getCoordinatesGlobal()->getX(F::S, 2),
            grid()->getCoordinatesGlobal()->getX(F::U, 0),
            grid()->getCoordinatesGlobal()->getX(F::V, 1),
            grid()->getCoordinatesGlobal()->getX(F::W, 2),
            grid()->getDomainSize()->getRe()/*, */
            /*grid()->getDomainSize()->getAlpha2()*/);

      }
    } else {

      int stride[3] = {1, 1, 1};

      write_hdf_3D(
          true,
          grid()->rankS(),
          MPI_Comm_c2f(grid()->comm()),
          grid()->nGlo(),
          grid()->getBCGlobal()->getBCL(),
          grid()->getBCGlobal()->getBCU(),
          grid()->nLoc(),
          grid()->bl(),
          grid()->bu(),
          grid()->sIndB(fType_),
          grid()->eIndB(fType_),
          grid()->getStencilWidths()->getLS(),
          grid()->np(),
          grid()->ib(),
          grid()->getShift(),
          (int)fType_+1,
          (int)fType_+1,
          count,
          ((F::S == fType_)?9:10)+7,
          stride,
          s_,
          grid()->getCoordinatesGlobal()->getX(F::S, 0),
          grid()->getCoordinatesGlobal()->getX(F::S, 1),
          grid()->getCoordinatesGlobal()->getX(F::S, 2),
          grid()->getCoordinatesGlobal()->getX(F::U, 0),
          grid()->getCoordinatesGlobal()->getX(F::V, 1),
          grid()->getCoordinatesGlobal()->getX(F::W, 2),
          grid()->getDomainSize()->getRe()/*, */
          /*grid()->getDomainSize()->getAlpha2() */);
    }
  }

  void read(const int count=0) {

    int vel_dir = static_cast<int>(fType_) + 1;
    read_hdf(
        grid()->rankS(),
        MPI_Comm_c2f(grid()->comm()),
        grid()->getBCGlobal()->getBCL(),
        grid()->getBCGlobal()->getBCU(),
        grid()->nLoc(),
        grid()->bl(),
        grid()->bu(),
        grid()->sIndB(fType_),
        grid()->eIndB(fType_),
        grid()->getStencilWidths()->getLS(),
        grid()->ib(),
        grid()->getShift(),
        vel_dir,
        count,
        ((F::S == fType_)?9:10)+7,
        s_);

    changed();
  }


public:

  constexpr F getType() {
    return fType_;
  }

  /// \name storage methods.
  /// \brief highly dependent on underlying storage should only be used by
  /// Operator or on top field implementer.
  /// @{

  OT getStorageSize() const {

    OT n = 1;
    for(int i=0; i<3; ++i)
      n *= grid()->nLoc(i)+SW::BU(i)-SW::BL(i)+1; // seems wrong: there a one was added for AMG, but it is not neede error seem to be in Impact there it should be (B1L+1:N1+B1U) probably has to be changed aganin for 3D

    return n;
  }

  void setStoragePtr(ST*  array) {
    s_ = array;
  }

  constexpr ScalarArray getRawPtr() {
    return s_;
  }

  constexpr const ST* getConstRawPtr() {
    return s_;
  }


  /// @}

  constexpr const Teuchos::RCP<const GridT>& grid() {
    return AbstractField<GridT>::grid_;
  }

  constexpr const MPI_Comm& comm() {
    return grid()->comm();
  }

  /// \name comunication methods.
  /// \brief highly dependent on underlying storage should only be used by
  /// Operator or on top field implementer.
  /// \{

  void changed(const int dir) const {
    exchangedState_[dir] = false;
  }
  void changed() const {
    for(int dir=0; dir<GridT::sdim; ++dir)
      changed(dir);
  }


  bool is_exchanged(const int dir) const {
    return exchangedState_[dir];
  }
  bool is_exchanged() const {
    bool all_exchanged = true;
    for(int dir=0; dir<GridT::sdim; ++dir)
      all_exchanged = all_exchanged && is_exchanged(dir);
    return all_exchanged;
  }

  /// \brief updates ghost layers
  void exchange(const int dir) const {
    int ones[3] = {0, 0, 0};
    if(!exchangedState_[dir]) {
      F_exchange(
        static_cast<int>(GridT::sdim),
        MPI_Comm_c2f(grid()->getProcGrid()->getCommWorld()),
        grid()->getProcGrid()->getRankL(),
        grid()->getProcGrid()->getRankU(),
        grid()->nLoc(),
        grid()->bl(),
        grid()->bu(),
        grid()->getBCLocal()->getBCL(),
        grid()->getBCLocal()->getBCU(),
        grid()->sInd(F::S),
        grid()->eInd(F::S),
//				grid()->sIndB(fType_), // should it work
//        grid()->eIndB(fType_),
        ones,
        grid()->nLoc(),
        1+dir,
        1+(int)fType_,
        s_);
      exchangedState_[dir] = true;
    }
  }
  void exchange() const {
    for(int dir=0; dir<GridT::sdim; ++dir)
      exchange(dir);
  }

  void setExchanged(const int dir) const {
    exchangedState_[dir] = true;
  }
  void setExchanged() const {
    for(int dir=0; dir<GridT::sdim; ++dir)
      changed(dir);
  }

  /// \}


protected:

  /// \name indexing
  /// @{

  /// \brief stride in Y direction
  //constexpr OT stride1() {
    //return grid()->nLoc(0)+SW::BU(0)-SW::BL(0)+1;
  //}

  /// \brief stride in Z direction
  //constexpr OT stride2() {
    //return
          //(grid()->nLoc(0)+SW::BU(0)-SW::BL(0)+1)*(
            //grid()->nLoc(1)+SW::BU(1)-SW::BL(1)+1);
  //}



  /// \brief comput index
  ///
  /// \param i index in x-direction
  /// \param j index in y-direction
  /// \param k index in z-direction
  constexpr OT index(const OT* const i) {
    return (i[0]-SW::BL(0)) +
      (i[1]-SW::BL(1))*stride1_ +
      (i[2]-SW::BL(2))*stride2_;
  }


  /// \brief comput index
  ///
  /// \param i index in x-direction
  /// \param j index in y-direction
  /// \param k index in z-direction
  constexpr OT index(const OT i, const OT j, const OT k) {
    return (i-SW::BL(0)) +
      (j-SW::BL(1))*stride1_ +
      (k-SW::BL(2))*stride2_;
  }


  /// \brief field access
  ///
  /// \param i index in x-direction
  /// \param j index in y-direction
  /// \param k index in z-direction
  ///
  /// \return const reference
  constexpr ST at(const OT i, const OT j, const OT k) {
    return s_[ index(i, j, k) ];
  }

  /// \brief field access
  ///
  /// \param i index in x-direction
  /// \param j index in y-direction
  /// \param k index in z-direction
  ///
  /// \return reference
  ST& at(const OT i, const OT j, const OT k)  {
    return s_[ index(i, j, k) ];
  }

  /// \brief field access
  ///
  /// \param i index coordinate
  constexpr ST at(const OT* const i) {
    return s_[ index(i) ];
  }
  /// \brief field access
  ///
  /// \param i index coordinate
  ST& at(const OT* const i) {
    return s_[ index(i) ];
  }

  /// \brief field access
  ///
  /// \param i index coordinate
  ST& at(const Teuchos::Tuple<const OT, 3>& i) {
    return s_[ index(i[0], i[1], i[2]) ];
  }

  /// \brief field access
  ///
  /// \param i index coordinate
  constexpr ST at(const Teuchos::Tuple<const OT, 3>& i) {
    return s_[ index(i[0], i[1], i[2]) ];
  }

public:

  /// \brief field access
  ///
  /// \param i index in x-direction
  /// \param j index in y-direction
  /// \param k index in z-direction
  ///
  /// \return const reference
  constexpr ST operator()(const OT i, const OT j, const OT k) {
    return at(i, j, k);
  }

  /// \brief field access
  ///
  /// \param i index in x-direction
  /// \param j index in y-direction
  /// \param k index in z-direction
  ///
  /// \return reference
  ST& operator()(const OT i, const OT j, const OT k)  {
    return at(i, j, k);
  }

  /// \brief field access
  ///
  /// \param i index coordinate
  constexpr ST operator()(const OT* const i) {
    return at(i);
  }
  /// \brief field access
  ///
  /// \param i index coordinate
  ST& operator()(const OT* const i) {
    return at(i);
  }

  /// \brief field access
  ///
  /// \param i index coordinate
  ST& operator()(const Teuchos::Tuple<const OT, 3>& i) {
    return at(i);
  }
  /// \brief field access
  ///
  /// \param i index coordinate
  constexpr ST operator()(const Teuchos::Tuple<const OT, 3>& i) {
    return at(i);
  }

}; // end of class ScalarField




/// \brief creates a scalar field(vector) belonging to a grid
///
/// \param grid scalar Vector Grid to which returned vector belongs
/// \param fType
/// \return scalar vector
/// \relates ScalarField
template<class GridT>
Teuchos::RCP<ScalarField<GridT> >
createScalarField(const Teuchos::RCP<const GridT >& grid, F fType=F::S) {

  return Teuchos::rcp(new ScalarField<GridT>(grid, true, fType));
}

///  @}

} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_SCALARFIELD_HPP
