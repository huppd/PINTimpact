#pragma once
#ifndef PIMPACT_GRIDSIZEGLOBAL_HPP
#define PIMPACT_GRIDSIZEGLOBAL_HPP


#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Utils.hpp"




namespace Pimpact {


/// \brief global grid size(independent of FieldType)
///
/// \todo remove tpara dimension and rm getters
/// \ingroup SpaceObject
template<class OrdinalT, int sdim >
class GridSizeGlobal : public Teuchos::Tuple<OrdinalT, 4> {

  template<class OT, int sd>
  friend Teuchos::RCP<const GridSizeGlobal<OT, sd> > createGridSizeGlobal(OT n1, OT n2, OT n3, OT nt);

  template<class OT, int sd>
  friend Teuchos::RCP<const GridSizeGlobal<OT, sd> > createGridSizeGlobal(const Teuchos::Tuple<OT, 4>& tuple);

protected:

  GridSizeGlobal(const Teuchos::Tuple<OrdinalT, 4>& gridSize):
    Teuchos::Tuple<OrdinalT, 4>(gridSize) {

    for(int i=0; i<sdim; ++i)
      assert(((*this)[i]-1)%2 == 0);
  };

public:

  constexpr const OrdinalT& get(const int i) const  {
    return (*this)[i];
  }

  void print(std::ostream& out=std::cout) const {
    out <<" --- GridSizeGlobal: " <<*this <<" ---\n";
  };


}; // end of class GridSizeGlobal



/// \brief create GridSize Global
/// \relates GridSizeGlobal
template<class OT, int sd>
Teuchos::RCP<const GridSizeGlobal<OT, sd> > createGridSizeGlobal(OT n1, OT n2, OT n3, OT nt=1) {

  Teuchos::Tuple<OT, 4> temp;

  temp[0] = n1;
  temp[1] = n2;
  temp[2] = n3;
  temp[3] = nt;

  return Teuchos::rcp(new GridSizeGlobal<OT, sd>(temp));
}



/// \brief create GridSize Global
/// \relates GridSizeGlobal
template<class OT, int sd>
Teuchos::RCP<const GridSizeGlobal<OT, sd> > createGridSizeGlobal(const Teuchos::Tuple<OT, 4>& to) {

  return Teuchos::rcp(new GridSizeGlobal<OT, sd>(to));
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDSIZEGLOBAL_HPP
