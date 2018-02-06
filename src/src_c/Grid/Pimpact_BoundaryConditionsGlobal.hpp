#pragma once
#ifndef PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP
#define PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \brief global boundary conditions
///
/// \tparam dim as soon time is another class dim ->sdim
/// \ingroup GridObject
/// \todo make int -> BC
template<int dim>
class BoundaryConditionsGlobal {

public:

  using TBC3 = const Teuchos::Tuple<BC, dim>;
  using Ti3 = const Teuchos::Tuple<int, dim>;

  template<int d>
  friend Teuchos::RCP<const BoundaryConditionsGlobal<d> > createBoudaryConditionsGlobal();

  template<int d>
  friend Teuchos::RCP<const BoundaryConditionsGlobal<d> >
  createBoudaryConditionsGlobal(const Teuchos::RCP<Teuchos::ParameterList>& pl);


protected:

  Ti3 BCL_int_;
  Ti3 BCU_int_;

  BoundaryConditionsGlobal(
    int BC1L, int BC1U,
    int BC2L, int BC2U,
    int BC3L, int BC3U) {

    BCL_int_[0] = BC1L;
    BCU_int_[0] = BC1U;
    BCL_int_[1] = BC2L;
    BCU_int_[1] = BC2U;
    BCL_int_[2] = BC3L;
    BCU_int_[2] = BC3U;

    if(4==dim) {
      BCL_int_[3] = static_cast<int>(BC::Periodic);
      BCU_int_[3] = static_cast<int>(BC::Periodic);
    }
  }

  BoundaryConditionsGlobal(TBC3 BCL_global, TBC3 BCU_global) {

    for(int i=0; i<3; ++i) {
      BCL_int_[i] = static_cast<int>(BCL_global[i]);
      BCU_int_[i] = static_cast<int>(BCU_global[i]);
    };
    if(4==dim) {
      BCL_int_[3] = static_cast<int>(BC::Periodic);
      BCU_int_[3] = static_cast<int>(BC::Periodic);
    }
  }

public:

  /// \name getter
  /// @{

  constexpr const int getBCL(const int dir) const {
    return BCL_int_[dir];
  }
  constexpr const int getBCU(const int dir) const {
    return BCU_int_[dir];
  }

  constexpr const int* getBCL() const {
    return BCL_int_.getRawPtr();
  }
  constexpr const int* getBCU() const {
    return BCU_int_.getRawPtr();
  }

  constexpr const Teuchos::Tuple<int, dim> periodic() const {

    Teuchos::Tuple<int, dim> periodic;
    for(int i=0; i<3; ++i) {
      if(getBCL(i)==BC::Periodic)
        periodic[i] = 1;
      else
        periodic[i] = 0;
    }
    if(4==dim) periodic[3] = 1;
    return periodic;
  }

  ///  @}

  /// \brief prints BC tuples
  ///
  /// \param out output stream
  void print(std::ostream& out=std::cout) const {
    out << "---BoundaryConditionsGlobal: ---\n";
    out << " BCL_global: " << BCL_int_ << "\n";
    out << " BCU_global: " << BCU_int_ << "\n";
  }

}; // end of class BoundaryConditionsGlobal




template<int dim=3>
Teuchos::RCP<const BoundaryConditionsGlobal<dim> >
createBoudaryConditionsGlobal(const Teuchos::RCP<Teuchos::ParameterList>& pl) {

  return Teuchos::rcp(new BoundaryConditionsGlobal<dim>(
        pl->get<int>("lower X", 1),
        pl->get<int>("upper X", 1),
        pl->get<int>("lower Y", 1),
        pl->get<int>("upper Y", 1),
        pl->get<int>("lower Z", 1),
        pl->get<int>("upper Z", 1)));
}




} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP
