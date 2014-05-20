#pragma once
#define PIMPACT_MLHELMHOLTDOP_HPP
#ifndef PIMPACT_MLHELMHOLTZOP_HPP



#include "ml_include.h"

#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"
//#include "Pimpact_MultiField.hpp"

//#include "Pimpact_HelmholtzOp.hpp"

//#include "Pimpact_LinearProblem.hpp"

#include "Pimpact_Space.hpp"



namespace Pimpact{


extern "C" {

void VF_extract_dof(
    const int& dimens,
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const sU,
    const int* const nU,
    const int* const sV,
    const int* const nV,
    const int* const sW,
    const int* const nW,
    const double* const phiiU,
    const double* const phiiV,
    const double* const phiiW,
    double* const phioU,
    double* const phioV,
    double* const phioW );
}


/// \ingroup BaseOperator
template< class Scalar, class Ordinal >
class MLHelmholtzOp {

public:

  typedef VectorField<Scalar,Ordinal>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal>  RangeFieldT;

  typedef Scalar* ScalarArray;

private:

  Teuchos::Tuple<ML*,3> mlObject_;
  Teuchos::Tuple<ML_Aggregate*,3> agg_object_;

  Teuchos::RCP<const Space<Ordinal> > space_;

  const int& dim() const { return( space_->filedS_->dim_   ); }

  Teuchos::Tuple<ScalarArray,3> x_;
  Teuchos::Tuple<ScalarArray,3> rhs_;

  Teuchos::Tuple<int,3> nloc_;

  int nLevels_;


public:

  MLHelmholtzOp(
      const Teuchos::RCP<const Space>& space,
      int nGrids=20 ):
    space_(space) {

    // comput nlo_
    int nTot = 0;
    for( int field=0; field<3; ++field ) {
      nloc_[field] = 1;
       for( int i=0; i<3; ++i )
         nloc_[field] *= space_->eInd(field)[i]-space_->sInd(field)[i];
       nTot += nloc_[field];
    }

    // init x
    x_[0] = new Scalar[ nTot ];
    x_[1] = x_[0] + nloc_[0];
    x_[2] = x_[1] + nloc_[1];

    // init rhs
    rhs_[0] = new Scalar[ nTot ];
    rhs_[1] = rhs_[0] + nloc_[0];
    rhs_[2] = rhs_[1] + nloc_[1];



    for( int i=0; i<dim(); ++i ) {
      ML_Create( &mlObject_[i], nGrids );

      ML_Init_Amatrix( mlObject_[i], 0,  nloc_[i], nloc_[i], NULL );
      ML_Set_Amatrix_Getrow( mlObject_[i], 0, getrow, NULL, 129 );
      ML_Set_Amatrix_Matvec( mlObject_[i], 0, matvec );


      ML_Aggregate_Create( &agg_object_[i]);
      ML_Aggregate_Set_MaxCoarseSize( agg_object_[i], 1 );

      nLevels_ = ML_Gen_MGHierarchy_UsingAggregation( mlObject_[i], 0,
                                                   ML_INCREASING, agg_object_[i] );

      ML_Gen_Smoother_GaussSeidel( mlObject_[i], ML_ALL_LEVELS, ML_PRESMOOTHER, 1, ML_DEFAULT );

      ML_Gen_Solver( mlObject_[i], ML_MGV, 0, nLevels_-1 );
    }

  };

  ~MLHelmholtzOp() {

    for( int i=0; i<3; ++i ) {
      ML_Aggregate_Destroy( &agg_object_[i] );
      ML_Destroy( &mlObject_[i] );
    }

  }

  int getrow(
      ML_Operator *Amat,
      int N_requested_rows,
      int requested_rows[],
      int allocated_space,
      int columns[],
      double values[],
      int row_lengths[] );

  int matvec(
      ML_Operator *Amat,
      int in_length,
      double p[],
      int out_length,
      double ap[] ) ;


  void apply(const DomainFieldT& x, RangeFieldT& y) const {

    VF_extract_dof(
        space_->dim(),
        space_->nLoc(),
        space_->bl(), space_->bu(),
        space_->sInd(0), space_->eInd(0),
        space_->sInd(1), space_->eInd(1),
        space_->sInd(2), space_->eInd(2),
        x.vec_[0], x.vec_[1], x.vec_[2],
        x_    [0], x_    [1], x_    [2] );

    VF_extract_dof(
        space_->dim(),
        space_->nLoc(),
        space_->bl(), space_->bu(),
        space_->sInd(0), space_->eInd(0),
        space_->sInd(1), space_->eInd(1),
        space_->sInd(2), space_->eInd(2),
        y.vec_[0], y.vec_[1], y.vec_[2],
        rhs_  [0], rhs_  [1], rhs_  [2] );

    for( int i=0; i<space_->dim(); ++i )
      ML_Iterate( mlObject_[i], *x_[i], *rhs_[i] );

  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class MLHelmholtzOp



/// \relates MLHelmholtzOp
template< class Scalar, class Ordinal>
Teuchos::RCP< MLHelmholtzOp<Scalar,Ordinal> > createMLHelmholtzOp(
    const Teuchos::RCP<LinearProblem<MultiField<VectorField<Scalar,Ordinal> > > > lap_prob ) {

  return( Teuchos::rcp( new MLHelmholtzOp<Scalar,Ordinal>( lap_prob ) ) );
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MLHELMHOLTDOP_HPP
