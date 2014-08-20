#pragma once
#define PIMPACT_MLHELMHOLTZOP_HPP
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

void VF_extract_dof_reverse(
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

void OP_innerhelmholtz(
    const int& m,
    const int* const ss,
    const int* const nn,
    const double& mulI,
    const double& multL,
    double* const x,
    double* const rhs );

void OP_HelmholtzGetRowEntries(
    const int& m,
    const int* const ss,
    const int* const nn,
    const double& mulI,
    const double& multL,
    const int& i,
    const int& j,
    const int& k,
    const int& maxRowEntries,
    int* const ic,
    int* const jc,
    int* const kc,
    double* const val,
    int& rowEntries );
}




/// \ingroup BaseOperator
/// \todo make Adat pretty
/// \todo make mulL better
template< class Scalar, class Ordinal >
class MLHelmholtzOp {

public:

  typedef VectorField<Scalar,Ordinal>  DomainFieldT;
  typedef VectorField<Scalar,Ordinal>  RangeFieldT;

  typedef Scalar* ScalarArray;

private:

  struct Adat {
    Teuchos::RCP< const Pimpact::Space<Scalar,Ordinal> > space_;
    int field_;
    double mulI_;
    double mulL_;
  };

  Teuchos::Tuple<ML*,3> mlObject_;
  Teuchos::Tuple<ML_Aggregate*,3> agg_object_;

  //  ML* mlObject_[3]; ML_Aggregate* agg_object_[3];

  Teuchos::Tuple<Adat,3> adat_;

  const int& dim() const { return( adat_[0].space_->dim() ); }
  Teuchos::RCP<const Space<Scalar,Ordinal> > space() const {return( adat_[0].space_ );}

  Teuchos::Tuple<int,3> nloc_;

  int nLevels_;

  Teuchos::Tuple<ScalarArray,3> sol_;
  Teuchos::Tuple<ScalarArray,3> rhs_;

public:

  MLHelmholtzOp(
      const Teuchos::RCP< const Space<Scalar,Ordinal> >& space=Teuchos::null,
      int nGrids=8,
      Scalar mulI=1.,
      Scalar mulL=1.,
      double tol = 1.e-3 ):
        nloc_(Teuchos::tuple(1,1,1)) {

    //    std::cout << "construct me baby\n";
    for( int field=0; field<space->dim(); ++field) {
      adat_[field].space_=space;
      adat_[field].field_ = field;
      adat_[field].mulI_ = mulI;
      adat_[field].mulL_ = mulL;
      //      adat_[field].space_->print();
    }
    //    space->print();

    //    MLHelmholtzOp::space_=space;
    // comput nlo_
    int nTot = 0;
    for( int field=0; field<3; ++field ) {
      for( int i=0; i<dim(); ++i ) {
        //        space()->print();
        //        std::cout << "\nnloc: " << nloc_[i]<< "\n";
        nloc_[field] *= space()->eInd(field)[i]-space()->sInd(field)[i]+1;
      }
      nTot += nloc_[field];
    }

    // init x
    sol_[0] = new Scalar[ nTot ];
    sol_[1] = sol_[0] + nloc_[0];
    sol_[2] = sol_[1] + nloc_[1];

    // init rhs
    rhs_[0] = new Scalar[ nTot ];
    rhs_[1] = rhs_[0] + nloc_[0];
    rhs_[2] = rhs_[1] + nloc_[1];



    for( int field=0; field<dim(); ++field ) {
      //      std::cout << "\nfield: " << field << "\tnloc: " << nloc_[field] << "\n";
      ML_Create( &mlObject_[field], nGrids );

      ML_Init_Amatrix( mlObject_[field], 0,  nloc_[field], nloc_[field], &adat_[field] );
      ML_Set_Amatrix_Getrow( mlObject_[field], 0, pimp_getrow, NULL, nloc_[field] );
      ML_Set_Amatrix_Matvec( mlObject_[field], 0, pimp_matvec );


      ML_Set_Tolerance( mlObject_[field], tol );

      //      ML_Set_OutputLevel( mlObject_[field], 10 );
//      ML_Set_ResidualOutputFrequency( mlObject_[field], -1 );
//      ML_Set_PrintLevel(10);

      ML_Aggregate_Create( &agg_object_[field]);
//      ML_Aggregate_Set_MaxCoarseSize( agg_object_[field], 4 );

      nLevels_ = ML_Gen_MGHierarchy_UsingAggregation( mlObject_[field], 0,
          ML_INCREASING, agg_object_[field] );

      //      std::cout << "\nnLeves: " << nLevels_ << "\n";

      //      ML_Gen_Smoother_Jacobi( mlObject_[field], ML_ALL_LEVELS, ML_PRESMOOTHER, 4, ML_DEFAULT);
      //      ML_Gen_Smoother_Jacobi( mlObject_[field], ML_ALL_LEVELS, ML_BOTH, 4, ML_DEFAULT);
      ML_Gen_Smoother_Cheby( mlObject_[field], ML_ALL_LEVELS, ML_BOTH, 0., 4);
      //                  ML_Gen_Smoother_GaussSeidel( mlObject_[field], ML_ALL_LEVELS, ML_PRESMOOTHER, 1, ML_DEFAULT);
      //            ML_Gen_Smoother_GaussSeidel( mlObject_[field], ML_ALL_LEVELS, ML_BOTH, 1, ML_DEFAULT);
      //      ML_Gen_Smoother_SymGaussSeidel( mlObject_[field], ML_ALL_LEVELS, ML_PRESMOOTHER, 1, ML_DEFAULT);
      //            ML_Gen_Smoother_SymGaussSeidel( mlObject_[field], ML_ALL_LEVELS, ML_BOTH, 1, ML_DEFAULT);
      ML_Gen_CoarseSolverSuperLU( mlObject_[field], nLevels_ );

      ML_Gen_Solver( mlObject_[field], ML_MGW, 0, nLevels_-1 );
    }

  };

  ~MLHelmholtzOp() {

    delete[] sol_[0];
    delete[] rhs_[0];
    for( int i=0; i<dim(); ++i ) {
      ML_Aggregate_Destroy( &agg_object_[i] );
      ML_Destroy( &mlObject_[i] );
    }

  }

private:

  static int pimp_getrow(
      ML_Operator *Amat,
      int nRequestedRows,
      int requestedRows[],
      int allocatedSpace,
      int columns[],
      double values[],
      int rowLengths[] ) {



    Adat* itemp;
    itemp  = (Adat *) ML_Get_MyGetrowData( Amat );
    int field = itemp->field_;
    double mulI = itemp->mulI_;
    double mulL = itemp->mulL_;
    Teuchos::RCP<const Pimpact::Space<Scalar,Ordinal> > space = itemp->space_;
    //    std::cout << "\ngetrow: # requested rows: "<< nRequestedRows << " # allocated space: " << allocatedSpace << "\n";
    //    space->print();

    int count = 0;
    int start, row;
    int i,j,k;

    int maxRowEntries = 0;
    for( int l=0; l<space->dim(); ++l )
      maxRowEntries += ( -space->bl()[l]+space->bu()[l]+1 ) ;
    maxRowEntries -= space->dim()-1;

    int rowEntries;

    int* ic = new int[maxRowEntries];
    int* jc = new int[maxRowEntries];
    int* kc = new int[maxRowEntries];
    Scalar* val = new Scalar[maxRowEntries];

    for( int l=0; l<nRequestedRows; l++ ) {
      if( allocatedSpace < count+maxRowEntries )
        return(0);
      start = count;
      row = requestedRows[l];
      //       if ( (row >= 0) || (row <= (129-1)) ) {
      //         columns[count] = row;
      //         values[count++] = 2.;
      //         if (row != 0) {
      //           columns[count] = row-1;
      //           values[count++] = -1.; }
      //         if (row != (129-1)) {
      //           columns[count] = row+1;
      //           values[count++] = -1.;
      //         }
      //       }
      k = row
          /(space->eInd(field)[0]-space->sInd(field)[0]+1)
          /(space->eInd(field)[1]-space->sInd(field)[1]+1)
          + space->sInd(field)[2];
      j = ( row
          -(space->eInd(field)[0]-space->sInd(field)[0]+1)
          *(space->eInd(field)[1]-space->sInd(field)[1]+1)
          *(k-space->sInd(field)[2]) )
                                            /(space->eInd(field)[0]-space->sInd(field)[0]+1)
                                            +space->sInd(field)[1];
      i = row
          -(space->eInd(field)[0]-space->sInd(field)[0]+1)
          *(space->eInd(field)[1]-space->sInd(field)[1]+1)
          *(k-space->sInd(field)[2])
          -(space->eInd(field)[0]-space->sInd(field)[0]+1)
          *(j-space->sInd(field)[1])
          +space->sInd(field)[0];
      //      std::cout << "\nrow: " << row <<" i: " << i << " j: " << j << " k: " << k << "\n";
      //      std::cout << "\nmax row entries: " << maxRowEntries << "\n";

      OP_HelmholtzGetRowEntries(
          field+1,
          space->sInd(field),
          space->eInd(field),
          mulI,
          mulL,
          i,
          j,
          k,
          maxRowEntries,
          ic,
          jc,
          kc,
          val,
          rowEntries );

      for( int m=0; m<rowEntries-1; ++m) {
        values[count] = val[m];
        columns[count++] = ic[m]-space->sInd(field)[0]
                                                    + (jc[m]-space->sInd(field)[1])*(space->eInd(field)[0]-space->sInd(field)[0]+1)
                                                    + (kc[m]-space->sInd(field)[2])*(space->eInd(field)[0]-space->sInd(field)[0]+1)*(space->eInd(field)[1]-space->sInd(field)[1]+1);
        //        std::cout << "\t\tic: " << ic[m] << " jc: " << jc[m] << " kc: " << kc[m] ;
        //        std::cout << "\tcol: "<< columns[count-1] << "\t";
        //        std::cout << "\tval: "<< val[m] << "\n";
      }
      rowLengths[l] = count - start;
      //      std::cout << "\n impact rowEntries: " << rowEntries;
      //      std::cout << "\npimpact rowEntries: " << rowLengths[l];
    }

    //    if( allocatedSpace<nRequestedRows ) return( 0 );
    //    for( int i=0; i<nRequestedRows; ++i) {
    //      values[i] = 1.;
    //      rowLengths[i] = 1;
    //      columns[i] = requestedRows[i];
    //    }
    return( 1 );

  }

  static int pimp_matvec(
      ML_Operator *Amat,
      int in_length,
      double* p,
      int out_length,
      double* ap ) {

    //    typedef double Scalar;
    //    typedef int Ordinal;

    //    std::cout << "\nmatvec: in length: " << in_length << " outlength: " << out_length << " ";
    //    std::cout << "\nmatvec\n";
    Adat* itemp;
    itemp  = (Adat *) ML_Get_MyGetrowData( Amat );
    int field = itemp->field_;
    double mulI = itemp->mulI_;
    double mulL = itemp->mulL_;
    Teuchos::RCP<const Pimpact::Space<double,int> > space = itemp->space_;
    //    space->print();
    //    int nloc = 1;
    //    for( int i=0; i<2;++i) {
    //      nloc*=-space->sInd(field)[i]+space->eInd(field)[i]+1;
    //    }
    //    std::cout << "nloc: " << nloc << "\n";

    OP_innerhelmholtz(
        field+1,
        space->sInd(field),
        space->eInd(field),
        mulI,
        mulL,
        p,
        ap );

    //    for( int i=0; i<out_length; ++i) {
    //      ap[i] = 1*p[i];
    //    }
    return( 0 );
  }

public:

  void apply(const DomainFieldT& x, RangeFieldT& y) const {
    x.exchange();

    VF_extract_dof(
        space()->dim(),
        space()->nLoc(),
        space()->bl(), space()->bu(),
        space()->sInd(0), space()->eInd(0),
        space()->sInd(1), space()->eInd(1),
        space()->sInd(2), space()->eInd(2),
        x.vec_[0], x.vec_[1], x.vec_[2],
        rhs_  [0], rhs_  [1], rhs_  [2] );

    VF_extract_dof(
        space()->dim(),
        space()->nLoc(),
        space()->bl(), space()->bu(),
        space()->sInd(0), space()->eInd(0),
        space()->sInd(1), space()->eInd(1),
        space()->sInd(2), space()->eInd(2),
        y.vec_[0], y.vec_[1], y.vec_[2],
        sol_  [0], sol_  [1], sol_  [2] );

    for( int i=0; i<dim(); ++i )
      ML_Iterate( mlObject_[i], sol_[i], rhs_[i] );
//            ML_Solve( mlObject_[i], nloc_[i], sol_[i], nloc_[i], rhs_[i]);
      //      ML_Solve_MGFull( mlObject_[i], sol_[i], rhs_[i] );
//      ML_Solve_MGV( mlObject_[i], sol_[i], rhs_[i] ); // not working?
    //          ML_Solve_AMGV( mlObject_[i], sol_[i], rhs_[i] );

    VF_extract_dof_reverse(
        space()->dim(),
        space()->nLoc(),
        space()->bl(), space()->bu(),
        space()->sInd(0), space()->eInd(0),
        space()->sInd(1), space()->eInd(1),
        space()->sInd(2), space()->eInd(2),
        sol_  [0], sol_  [1], sol_  [2],
        y.vec_[0], y.vec_[1], y.vec_[2] );

    y.changed();

  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

}; // end of class MLHelmholtzOp


/// \relates MLHelmholtzOp
template< class Scalar, class Ordinal>
Teuchos::RCP< MLHelmholtzOp<Scalar,Ordinal> > createMLHelmholtzOp(
    const Teuchos::RCP< const Space<Scalar,Ordinal> >& space,
    int nGrids=20,
    Scalar mulI=1.,
    Scalar mulL=1.,
    double tol=1.e-6 ) {

  return(
      Teuchos::rcp(
          new MLHelmholtzOp<Scalar,Ordinal>( space, nGrids, mulI, mulL, tol ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MLHELMHOLTDOP_HPP
