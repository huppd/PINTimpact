#pragma once
#ifndef mlguide_hpp
#define mlguide_hpp

#include <math.h>
#include <iostream>
#include "ml_include.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"
#ifdef ML_MPI
#include "mpi.h"
#endif


//extern int Poisson_getrow(ML_Operator *A_data, int N_requested_rows, int requested_rows[],
//    int allocated_space, int columns[], double values[], int row_lengths[]);
//
//extern int Poisson_matvec(ML_Operator *A_data, int in_length, double p[], int out_length,
//    double ap[]);

namespace bla{

template<class S,class O>
class MyML {

  Teuchos::Tuple<ML*,3> ml_object;
  Teuchos::Tuple<ML_Aggregate*,3> agg_object;

public:

  MyML( int N_grids=20 ) {

    for( int dim=0; dim<2; ++dim ) {
      ML_Create         (&ml_object[dim], N_grids);

      ML_Init_Amatrix      (ml_object[dim], 0,  129, 129, &dim);
      ML_Set_Amatrix_Getrow(ml_object[dim], 0,  Poisson_getrow, NULL, 129);
      ML_Set_Amatrix_Matvec(ml_object[dim], 0,  Poisson_matvec);
      ML_Set_PrintLevel(10);

      ML_Aggregate_Create(&agg_object[dim]);
      ML_Aggregate_Set_MaxCoarseSize(agg_object[dim],1);
      int N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_object[dim], 0,
          ML_INCREASING, agg_object[dim]);
      /******** Begin code to set a Jacobi smoother ******/

      //   ML_Gen_Smoother_Jacobi(ml_object, ML_ALL_LEVELS, ML_PRESMOOTHER, 10, ML_DEFAULT);
      ML_Gen_Smoother_GaussSeidel(ml_object[dim], ML_ALL_LEVELS, ML_PRESMOOTHER, 10, ML_DEFAULT);
      //   ML_Gen_Smoother_BlockGaussSeidel(ml_object, ML_ALL_LEVELS, ML_PRESMOOTHER, 1, ML_DEFAULT, 3 );

      /******** End code to set a Jacobi smoother ******/

      ML_Gen_Solver (ml_object[dim], ML_MGV, 0, N_levels-1);
    }
  }
  void apply( double* sol, double* rhs ) {

    for(int dim=0; dim<2; ++dim) {
      ML_Iterate(ml_object[dim], sol, rhs);
    }
  }
  ~MyML() {
    for(int dim=0; dim<2; ++dim) {
      ML_Aggregate_Destroy(&agg_object[dim]);
      ML_Destroy(&ml_object[dim]);
    }
  }

private:
  static int Poisson_getrow(ML_Operator *A_data, int N_requested_rows, int requested_rows[],
      O allocated_space, int columns[], double values[], int row_lengths[]) {

    //    std::cout << "getrow\n";
    int count = 0, i, start, row;

    for( i= 0; i < N_requested_rows; i++) {
      if (allocated_space < count+3) return(0);
      start = count;
      row = requested_rows[i];
      if ( (row >= 0) || (row <= (129-1)) ) {
        columns[count] = row; values[count++] = 2.;
        if (row != 0) { columns[count] = row-1; values[count++] = -1.; }
        if (row != (129-1)) { columns[count] = row+1; values[count++] = -1.; }
      }
      row_lengths[i] = count - start;
    }
    return(1);
  }
  static int Poisson_matvec(ML_Operator *A_data, int in_length, double p[], int out_length,
      double ap[]) {

    std::cout << "blalblamatvec\n";
    for( O i = 0; i < 129; i++ ) {
      ap[i] = 2*p[i];
      if (i != 0) ap[i] -= p[i-1];
      if (i != (129-1)) ap[i] -= p[i+1];
    }
    return( 0 );
  }
};

template<class S,class O>
Teuchos::RCP<MyML<S,O> > createMyML( int N_grids ) {
  return(
      Teuchos::rcp(new bla::MyML<S,O>(N_grids) ) );
}

}




#endif
