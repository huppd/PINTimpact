#pragma once
#ifndef PIMPACT_H5WRITEROP_HPP
#define PIMPACT_H5WRITEROP_HPP

#include "Pimpact_Types.hpp"

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"

namespace Pimpact{


extern "C" {

void write_hdf5_2D(
        const int& COMM_CART,
        const int* const M,
        const int* const BC_L_global,
        const int* const BC_U_global,
        const int* const N,
        const int* const bL,
        const int* const bU,
        const int* const SS,
        const int* const NN,
        const int* const ls,
        const int* const NB,
        const int* const iB,
        const int* const iShift,
        const int& ftype,
        const int& filecount,
        const int& filenamelen,
        const double* const phi,
        const double* const y1p,
        const double* const y2p,
//        const double* const y3p,
        const double& Re,
        const double& alpha2 );

}


/// \brief output operator.
/// \ingroup BaseOperator
template<class Scalar,class Ordinal,int dimension=3>
class H5WriterOp {

protected:

  Teuchos::RCP< const Space<Scalar,Ordinal,dimension> > space_;

public:

  typedef ScalarField<Scalar,Ordinal,dimension>  SF;

  H5WriterOp( const Teuchos::RCP<const Space<Scalar,Ordinal,dimension> >& space ):
    space_(space) {

  };


  void write(const SF& x, int filecount=0 ) const {


  }


};


/// \relates H5WriterOp
template< class S, class O, int d=3 >
Teuchos::RCP< H5WriterOp<S,O,d> > createH5WriterOp(
    const Teuchos::RCP<const Space<S,O,d> >& space ) {
  return( Teuchos::rcp( new H5WriterOp<S,O,d>( space ) ) );
}

} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_H5WRITEROP_HPP
