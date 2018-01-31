#pragma once
#ifndef PIMPACT_COORDINATESGLOBAL_HPP
#define PIMPACT_COORDINATESGLOBAL_HPP


#include <cmath>
#include <ostream>

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_DomainSize.hpp"
#include "Pimpact_GridSizeGlobal.hpp"
#include "Pimpact_Stencil.hpp"
#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \brief global grid coordinates
///
/// Coordinate | index   | domain
/// -----------| --------| -------
/// xS         | 1..nGlo | 0..L
/// xV         | 0..nGlo | 0.5..L+0.5
///
/// \tparam ScalarT scalar type
/// \tparam OrdinalT index type
/// \tparam dim computational dimension as soon as Time is own class -> sdim
///
/// \note: - local processor-block coordinates and grid spacings are
///          automatically derived from global grid
///        - dy3p = dy3w = 1. for 2D (may simplify programming)
///        - ensure that for all i
///             y1p(i) <y1p(i+1)
///             y1u(i) <y1u(i+1)
///             y1p(i) <y1u(i) <y1p(i+1)
///                etc.
///        - code is tested only for
///             y1p(1) = 0.
///             y1p(M1) = L1
///                etc.
///
/// \todo make nice interface for getter
/// \relates CoordinatesGlobal
/// \ingroup SpaceObject
template<class ScalarT, class OrdinalT, int dim>
class CoordinatesGlobal {

  template<class ST, class OT, int sdT, int dT>
  friend Teuchos::RCP<const CoordinatesGlobal<ST, OT, dT> > createCoordinatesGlobal(
    const Teuchos::RCP<const GridSizeGlobal<OT, sdT> >& gridSize,
    const Teuchos::RCP<const DomainSize<ST, sdT> >& domainSize,
    const Teuchos::Tuple<Teuchos::RCP<Teuchos::ParameterList>, 3 >& gridStretching);

  template<class ST, class OT, int sdT, int dT>
  friend Teuchos::RCP<const CoordinatesGlobal<ST, OT, dT> > createCoordinatesGlobal(
    const Teuchos::RCP<const GridSizeGlobal<OT, sdT> >& gridSize,
    const Teuchos::RCP<const CoordinatesGlobal<ST, OT, dT> >& coordinates);

protected:

  using AS = Array<ScalarT, OrdinalT, 1>;
  using AV = Array<ScalarT, OrdinalT, 0>;

  using TAS = const Teuchos::Tuple<AS, dim >;
  using TAV = const Teuchos::Tuple<AV, dim >;

  TAS xS_;
  TAV xV_;

  Teuchos::Tuple<Teuchos::RCP<Teuchos::ParameterList>, 3 > stretchPara_;

  /// \name coordinate stretchings
  /// @{

  /// \brief equidistant grid
  ///
  /// \f[\mathrm{ x[i] = \frac{iL}{M-1} - x0}\f]
  ///
  /// \param[in] i
  /// \param[in] L
  /// \param[in] M
  /// \param[in] x0
  /// \param[out] x
  void coord_equi(const ScalarT i, const ScalarT L, const ScalarT M, const ScalarT x0, ScalarT& x) {
    x  = i*L/(M-1.) - x0;
  }


  /// \brief coordinate stretching for parabolas
  ///
  /// \f[ \mathrm{ x[i] = L\left(\frac{ i^2 }{ (M-1)^2 } + 2\alpha \frac{i}{M-1}  \right)\frac{1}{1+2\alpha)} - x0 } \f]
  ///
  /// \param[in] i index
  /// \param[in] L length of domain
  /// \param[in] M number of global grid points
  /// \param[in] x0 origin
  /// \param[in] alpha parameter for parabola alpha=0 very parabolic alpha>>0 equidistant
  /// \param[out] x coordinate
  void coord_parab(const ScalarT i, const ScalarT L, const ScalarT M, const ScalarT x0, const ScalarT alpha, ScalarT& x) {
    x  = L*(std::pow(i, 2)/std::pow(M-1., 2) + 2.*alpha*i/(M-1.))/(1.+2.*alpha) - x0;
  }


  /// \brief coordinate stretching for parabolas at the end
  ///
  /// \f[ \mathrm{ x[i] = L\left(\frac{ i^2 }{ (M-1)^2 } + 2\alpha \frac{i}{M-1}  \right)\frac{1}{1+2\alpha)} - x0 } \f]
  ///
  /// \param[in] i index
  /// \param[in] L length of domain
  /// \param[in] M number of global grid points
  /// \param[in] x0 origin
  /// \param[in] alpha parameter for parabola alpha=0 very parabolic alpha>>0 equidistant
  /// \param[out] x coordinate
  void coord_parab_end(const ScalarT i, const ScalarT L, const ScalarT M, const ScalarT x0, const ScalarT alpha, ScalarT& x) {
    x  = L*(-std::pow(i-M+1, 2)/std::pow(M-1., 2) + 1. + 2.*alpha*i/(M-1.))/(1.+2.*alpha) - x0;
  }



  /// \brief cos stretching
  ///
  /// \param[in] i index
  /// \param[in] L length of domain
  /// \param[in] M number of global grid points
  /// \param[in] x0 origin
  /// \param[in] iML ???
  /// \param[in] iMU ???
  /// \param[in] i0L ???
  /// \param[in] i0U ???
  /// \param[out] x coordinate
  ///
  /// \f[\mathrm{ wL = \frac\pi{2(i0L + iML-1)} }\f]
  /// \f[\mathrm{ xL = \frac{\cos(wL *i0L)}{wL} }\f]
  /// \f[\mathrm{ wU = \frac\pi{2(M + i0U - iML)} }\f]
  /// \f[\mathrm{ xU = \frac{\cos(wU * i0U)}{wU} }\f]
  /// \f[\mathrm{ x[i] = i*L/(M-1) - x0 }\f]
  ///
  /// \note
  ///	- i0L >= 0., i0U >= 0., iML >= 1 and iMU <= M is already tested.
  /// - so following is satisfied wL, wU >= 0..
  /// - identical to coord_tan except of std::cos functions.
  void coord_cos(
    const ScalarT i,
    const ScalarT L,
    const ScalarT M,
    const ScalarT x0,
    const ScalarT iML,
    const ScalarT iMU,
    const ScalarT i0L,
    const ScalarT i0U,
    ScalarT& x) {

    ScalarT wL;
    ScalarT wU;
    ScalarT xL;
    ScalarT xU;

    ScalarT pi = 4.*std::atan(1.);

    //--- parameters for grid stretching
    if(iML<=1.) {
      wL = 0.;
      xL = 0.;
    } else {
      wL = pi/(2.*(iML + i0L - 1.));
      xL = std::cos(wL*i0L)/wL;
    }

    if(iMU>=M) {
      wU = 0.;
      xU = 0.;
    } else {
      wU = pi/(2.*(M - iMU + i0U));
      xU = std::cos(wU*i0U)/wU;
    }

    //--- coordinates in the physical space
    if(i<iML && wL!= 0.) {
      if((i + i0L - 1.) <0.)
        // mirroring of the function
        x = - xL + std::cos(wL*(i + i0L - 1.))/wL;
      else
        x =   xL - std::cos(wL*(i + i0L - 1.))/wL;
    } else if(i>iMU && wU!=0.) {
      if((M - i + i0U) <0.)
        // mirroring of the function
        x = (2. - std::cos(wU*(i - i0U - M)))/wU;
      else
        x =  std::cos(wU*(i - i0U - M))/wU;
      x = x + xL - iML + iMU;
    } else
      x = i + xL - iML;


    //--- Normalization
    x *= L/(xL + xU - iML + iMU);
    x -= x0;
  }

  ///  @}



  /// \brief helper function getting number for switch statement
  /// from name
  /// \param[in] name input name
  /// \return according int number
  int string2int(const std::string& name) {
    std::string lcName = name;
    std::transform(lcName.begin(), lcName.end(), lcName.begin(), ::tolower);
    if("none" == lcName) return 0;
    else if("parabola" == lcName) return 1;
    else if("parab" == lcName) return 1;
    else if("para" == lcName) return 1;
    else if("para end" == lcName) return 3;
    else if("cos" == lcName) return 2;
    else {
      const bool& Stertch_Type_not_known = true;
      TEUCHOS_TEST_FOR_EXCEPT(Stertch_Type_not_known);
    }
    return 0;
  }



  /// \brief
  ///
  /// \param[in] gridSize
  /// \param[in] domainSize
  /// \param[in] stretchPara
  template<int sd>
  CoordinatesGlobal(
    const Teuchos::RCP<const GridSizeGlobal<OrdinalT, sd> >& gridSize,
    const Teuchos::RCP<const DomainSize<ScalarT, sd> >& domainSize,
    const Teuchos::Tuple<Teuchos::RCP<Teuchos::ParameterList>, 3 >& stretchPara):
    stretchPara_(stretchPara) {

    for(int dir=0; dir<dim; ++dir) {

      OrdinalT M = std::max(gridSize->get(dir), 1); // max need in the case of nf=0
      ScalarT Ms = std::max(M, 1);

      xS_ [dir] = AS(M);
      xV_ [dir] = AV(M);

      if(dir<3) {

        ScalarT L  = domainSize->getSize(dir);
        ScalarT x0 = domainSize->getOrigin(dir);

        int stretchType = string2int(stretchPara_[dir]->get<std::string>("Stretch Type", "none"));
        for(OrdinalT i=1; i<=M; ++i) {
          ScalarT is = i-1;

          switch(stretchType) {
          case 0:
            coord_equi(is, L, Ms, x0, xS_[dir][i]);
            break;
          case 1:
            coord_parab(is, L, Ms, x0, stretchPara_[dir]->get<ScalarT>("alpha", 0.5), xS_[dir][i]);
            break;
          case 2:
            coord_cos(
              is+1.,
              L,
              Ms,
              x0,
              stretchPara_[dir]->get<ScalarT>("N metr L", 1.),
              stretchPara_[dir]->get<ScalarT>("N metr U", M ),
              stretchPara_[dir]->get<ScalarT>("x0 L", 0.),
              stretchPara_[dir]->get<ScalarT>("x0 U", 0.),
              xS_[dir][i]);
            break;
          case 3:
            coord_parab_end(is, L, Ms, x0, stretchPara_[dir]->get<ScalarT>("alpha", 0.5), xS_[dir][i]);
            break;
          default:
            coord_equi(is, L, Ms, x0, xS_[dir][i]);
            break;
          }
        }
        for(OrdinalT i=0; i<=M; ++i) {
          ScalarT is = static_cast<ScalarT>(i) - 0.5;
          switch(stretchType) {
          case 0:
            coord_equi(is, L, Ms, x0, xV_[dir][i]);
            break;
          case 1:
            coord_parab(is, L, Ms, x0, stretchPara_[dir]->get<ScalarT>("alpha", 0.5), xV_[dir][i]);
            break;
          case 2:
            coord_cos(
              is+1.,
              L,
              Ms,
              x0,
              stretchPara_[dir]->get<ScalarT>("N metr L", 1.),
              stretchPara_[dir]->get<ScalarT>("N metr U", M ),
              stretchPara_[dir]->get<ScalarT>("x0 L", 0.),
              stretchPara_[dir]->get<ScalarT>("x0 U", 0. ),
              xV_[dir][i]);
            break;
          case 3:
            coord_parab_end(is, L, Ms, x0, stretchPara_[dir]->get<ScalarT>("alpha", 0.5), xV_[dir][i]);
            break;
          default:
            coord_equi(is, L, Ms, x0, xV_[dir][i]);
            break;
          }
        }
      } else if(3==dir) {
        // in time direction no stretching is considered as long as
        // time-periodic problems are considered equidistant should be best

        ScalarT pi2 = 8.*std::atan(1.);

        for(OrdinalT i=1; i<=M; ++i) {
          ScalarT is = i-1;
          xS_ [dir][i] = is*pi2/(Ms);
        }
        for(OrdinalT i=0; i<=M; ++i) {
          ScalarT is = static_cast<ScalarT>(i) - 0.5;
          xV_ [dir][i] = is*pi2/(Ms-1.);
        }
      }
    }
  }


  /// \brief constructor from fine grid
  ///
  /// \param[in] gridSizeC
  /// \param[in] coordinatesF
  template<int sd>
  CoordinatesGlobal(
    const Teuchos::RCP<const GridSizeGlobal<OrdinalT, sd> >& gridSizeC,
    const Teuchos::RCP<const CoordinatesGlobal<ScalarT, OrdinalT, dim> >& coordinatesF) {

    for(int dir=0; dir<dim; ++dir) {

      OrdinalT Mc = std::max(gridSizeC->get(dir), 1);
      OrdinalT Mf = coordinatesF->xS_[dir].NN();

      if(Mc==Mf) {
        xS_[dir] = coordinatesF->xS_[dir];
        xV_[dir] = coordinatesF->xV_[dir];
      } else {

        xS_[dir] = AS(Mc);
        xV_[dir] = AV(Mc);

        OrdinalT d = 1;

        if(dir<3) {
          if(Mc>1) // shouldn't be necessary when strategy makes its job correct. maybe throw exception
            d = (Mf - 1)/(Mc - 1);
        } else {
          if(Mc>0) // shouldn't be necessary when strategy makes its job correct. maybe throw exception
            d = Mf / Mc;
        }

        for(OrdinalT j=1; j<=Mc; ++j)
          xS_[dir][j] = coordinatesF->xS_[dir][(j-1)*d+1];

        for(OrdinalT j=1; j<Mc; ++j)
          xV_[dir][j] = coordinatesF->xS_[dir][j*d];

        // otherwise restriction of the boundary is wrong
        xV_[dir][0 ] = coordinatesF->xV_[dir][0];
        xV_[dir][Mc] = coordinatesF->xV_[dir][Mf];

        // for many grids this seems better for less grids the above seems to be better ~4
        // makes restriction on boundaries more complicated but Operators on coarse grid smoother
        //const double& alpha = 0.5; // 0.5 increase convergence for unstreched grids beyond maxGirds>2 but for stretch even worse for maxGrids=2
        //xV_[dir][0 ] = (1.+alpha)*xS_[dir][1] -alpha*xV_[dir][1];
        //xV_[dir][Mc] = (1.+alpha)*xS_[dir][Mc]-alpha*xV_[dir][Mc-1];
        {
          //ScalarT is = static_cast<ScalarT>(i) - 0.5;
          //switch(stretchType) {
          //case 0:
          //coord_equi(-0.5,   L, Mc, x0, xV_[dir][0]);
          //coord_equi(Mc-0.5, L, Mc, x0, xV_[dir][Mc]);
          //break;
          //case 1:
          //coord_parab(-0.5,   L, Mc, x0, stretchPara_[dir]->get<ScalarT>("alpha", 0.5), xV_[dir][0]);
          //coord_parab(Mc-0.5, L, Mc, x0, stretchPara_[dir]->get<ScalarT>("alpha", 0.5), xV_[dir][Mc]);
          //break;
          //case 2:
          //coord_cos(
          //0.5,
          //L,
          //Mc,
          //x0,
          //stretchPara_[dir]->get<ScalarT>("N metr L", 1.),
          //stretchPara_[dir]->get<ScalarT>("N metr U", M ),
          //stretchPara_[dir]->get<ScalarT>("x0 L", 0.),
          //stretchPara_[dir]->get<ScalarT>("x0 U", 0. ),
          //xV_[dir][0]);
          //coord_cos(
          //Mc+0.5,
          //L,
          //Mc,
          //x0,
          //stretchPara_[dir]->get<ScalarT>("N metr L", 1.),
          //stretchPara_[dir]->get<ScalarT>("N metr U", M ),
          //stretchPara_[dir]->get<ScalarT>("x0 L", 0.),
          //stretchPara_[dir]->get<ScalarT>("x0 U", 0. ),
          //xV_[dir][0]);
          //break;
          //default:
          //coord_equi(-0, 5,   L, Mc, x0, xV_[dir][0]);
          //coord_equi(MC-0.5, L, Mc, x0, xV_[dir][Mc]);
          //break;
          //}
        }
      }
    }
  }

public:

  /// \name getter
  /// @{

  constexpr const ScalarT* getX(const F ftype, const int dir) {
    return (F::S==ftype || dir!=ftype)?
      xS_[dir].get() : xV_[dir].get();
  }

  constexpr const ScalarT operator()(const F ftype, const int dir, const OrdinalT i) {
    return (F::S==ftype || dir!=ftype)?
      xS_[dir][i] : xV_[dir][i];
  }

  constexpr const ScalarT getX(const F ftype, const int dir, const OrdinalT i) {
    return (F::S==ftype || dir!=ftype)?
      xS_[dir][i] : xV_[dir][i];
  }

  constexpr const Teuchos::Tuple<Teuchos::RCP<Teuchos::ParameterList> , 3>& getStretchParameter() const {
    return stretchPara_;
  }

  ///  @}

  void print(std::ostream& out=std::cout) const {

    for(int i=0; i<dim; ++i) {
      out <<"Global coordinates of scalars in dir: " <<static_cast<ECoord>(i) <<"\n";
      xS_[i].print(out);
    }

    for(int i=0; i<dim; ++i) {
      out <<"Global coordinates of velocity in dir: " <<static_cast<ECoord>(i) <<"\n";
      xV_[i].print(out);
    }
  };


}; // end of class CoordinatesGlobal



/// \brief create Grid coordinates Global
/// \relates CoordinatesGlobal
///
/// \tparam ST
/// \tparam OT
/// \tparam d
/// \param gridSize
/// \param domainSize
/// \param gridStretching
///
/// \return
template<class ST, class OT, int sd, int d>
Teuchos::RCP<const CoordinatesGlobal<ST, OT, d> >
createCoordinatesGlobal(
  const Teuchos::RCP<const GridSizeGlobal<OT, sd> >& gridSize,
  const Teuchos::RCP<const DomainSize<ST, sd> >& domainSize,
  const Teuchos::Tuple<Teuchos::RCP<Teuchos::ParameterList> , 3>& gridStretching) {

  return Teuchos::rcp(
      new CoordinatesGlobal<ST, OT, d>(
        gridSize,
        domainSize,
        gridStretching));
}



/// \brief creates coarse coordinates
///
/// \tparam ST
/// \tparam OT
/// \tparam d
/// \param gridSize
/// \param coordinates
///
/// \return
template<class ST, class OT, int sd, int d>
Teuchos::RCP<const CoordinatesGlobal<ST, OT, d> >
createCoordinatesGlobal(
    const Teuchos::RCP<const GridSizeGlobal<OT, sd> >& gridSize,
    const Teuchos::RCP<const CoordinatesGlobal<ST, OT, d> >& coordinates) {

  return Teuchos::rcp(new CoordinatesGlobal<ST, OT, d>(gridSize, coordinates));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COORDINATESGLOBAL_HPP
