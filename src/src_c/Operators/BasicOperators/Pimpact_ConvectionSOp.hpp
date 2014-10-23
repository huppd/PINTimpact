#pragma once
#ifndef PIMPACT_CONVECTIONSOP_HPP
#define PIMPACT_CONVECTIONSOP_HPP


#include "Pimpact_Types.hpp"

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_ScalarField.hpp"



namespace Pimpact {


extern "C" {

void OP_convection(
    const int& dimens,
    const int* const N,
    const int* const bL,
    const int* const bU,
    const int* const nL,
    const int* const nU,
    const int* const SS,
    const int* const NN,
    const double* const c1D,
    const double* const c2D,
    const double* const c3D,
    const double* const c1U,
    const double* const c2U,
    const double* const c3U,
    const double* const phiU,
    const double* const phiV,
    const double* const phiW,
    const double* const phi,
    const double* nlu,
    const double& mul );

//void OP_nonlinear(
//    double* const phi1U, double* const phi1V, double* const phi1W,
//    double* const phi2U, double* const phi2V, double* const phi2W,
//    double* const nl1,   double* const nl2,   double* const nl3,
//    const double& mul );

}


/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup BaseOperator
template<class Scalar,class Ordinal, int dimension=3>
class ConvectionSOp {

public:

  typedef Teuchos::Tuple< Teuchos::RCP< ScalarField<Scalar,Ordinal,dimension> >, 3 > FluxFieldT;
  typedef ScalarField<Scalar,Ordinal,dimension>  DomainFieldT;
  typedef ScalarField<Scalar,Ordinal,dimension>  RangeFieldT;

protected:

  typedef const Teuchos::Tuple<Scalar*,3> TO;

  Teuchos::RCP<const Space<Scalar,Ordinal,dimension> > space_;

  TO cSD_;
  TO cVD_;

  TO cSU_;
  TO cVU_;

public:

  ConvectionSOp( const Teuchos::RCP<const Space<Scalar,Ordinal,dimension> >& space  ):
    space_(space) {

    for( int i=0; i<3; ++i ) {
      Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->nu(i) - space_->nl(i) + 1);

      cSD_[i] = new Scalar[ nTemp ];
      if( i<space_->dim() )
        FD_getDiffCoeff(
            space_->rankST(),
            space_->nLoc(i),
            space_->bl(i),
            space_->bu(i),
            space_->nl(i),
            space_->nu(i),
            space_->getDomain()->getBCLocal()->getBCL(i),
            space_->getDomain()->getBCLocal()->getBCU(i),
            space_->getShift(i),
            int(EField::S)+1,
            i+1,
            1,
            -1,
            true,
            space_->getFieldSpace()->getDimNcbC(i),
            space_->getFieldSpace()->getNcbC(i),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            cSD_[i] );

      cSU_[i] = new Scalar[ nTemp ];
      if( i<space_->dim() )
        FD_getDiffCoeff(
            space_->rankST(),
            space_->nLoc(i),
            space_->bl(i),
            space_->bu(i),
            space_->nl(i),
            space_->nu(i),
            space_->getDomain()->getBCLocal()->getBCL(i),
            space_->getDomain()->getBCLocal()->getBCU(i),
            space_->getShift(i),
            int(EField::S)+1,
            i+1,
            1,
            +1,
            true,
            space_->getFieldSpace()->getDimNcbC(i),
            space_->getFieldSpace()->getNcbC(i),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            cSU_[i] );

      cVD_[i] = new Scalar[ nTemp ];
      if( i<space_->dim() )
        FD_getDiffCoeff(
            space_->rankST(),
            space_->nLoc(i),
            space_->bl(i),
            space_->bu(i),
            space_->nl(i),
            space_->nu(i),
            space_->getDomain()->getBCLocal()->getBCL(i),
            space_->getDomain()->getBCLocal()->getBCU(i),
            space_->getShift(i),
            1,
            i+1,
            1,
            -1,
            true,
            space_->getFieldSpace()->getDimNcbC(i),
            space_->getFieldSpace()->getNcbC(i),
            space_->getCoordinatesLocal()->getX( i, i ),
            space_->getCoordinatesLocal()->getX( i, i ),
            cVD_[i] );

      cVU_[i] = new Scalar[ nTemp ];
      if( i<space_->dim() )
        FD_getDiffCoeff(
            space_->rankST(),
            space_->nLoc(i),
            space_->bl(i),
            space_->bu(i),
            space_->nl(i),
            space_->nu(i),
            space_->getDomain()->getBCLocal()->getBCL(i),
            space_->getDomain()->getBCLocal()->getBCU(i),
            space_->getShift(i),
            1,
            i+1,
            1,
            +1,
            true,
            space_->getFieldSpace()->getDimNcbC(i),
            space_->getFieldSpace()->getNcbC(i),
            space_->getCoordinatesLocal()->getX( i, i ),
            space_->getCoordinatesLocal()->getX( i, i ),
            cVU_[i] );
    }

  };

  void assignField( const RangeFieldT& mv ) {};



  void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, Scalar mul=0. ) const {

    int m = (int)z.fType_;

    TEUCHOS_TEST_FOR_EXCEPTION(
         z.fType_ != y.fType_,
         std::logic_error,
         "Pimpact::ConvectionSOP can only be applied to same fieldType !!!\n");




    for( int i =0; i<space_->dim(); ++i ) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          x[i]->fType_ != y.fType_,
          std::logic_error,
          "Pimpact::ConvectionSOP can only be applied to same fieldType !!!\n");
    }

    for( int vel_dir=0; vel_dir<space_->dim(); ++vel_dir )
      x[vel_dir]->exchange();

    y.exchange();

    // why not use default parameter 1? because one has to init z equal 0.
    if( std::abs(mul) < 1.e-16 ) {
      z.init( 0. );
      mul = 1.;
    }

    OP_convection(
        space_->dim(),
        space_->nLoc(),
        space_->bl(),
        space_->bu(),
        space_->nl(),
        space_->nu(),
        space_->sInd(m),
        space_->eInd(m),
        (m==0)?cVD_[0]:cSD_[0],
        (m==1)?cVD_[1]:cSD_[1],
        (m==2)?cVD_[2]:cSD_[2],
        (m==0)?cVU_[0]:cSU_[0],
        (m==1)?cVU_[1]:cSU_[1],
        (m==2)?cVU_[2]:cSU_[2],
        x[0]->s_,
        x[1]->s_,
        x[2]->s_,
        y.s_,
        z.s_,
        mul );

    z.changed();

  }

  void print( std::ostream& out=std::cout ) const {
     out << " --- scalar stencil: ---";
     for( int i=0; i<3; ++i ) {
       out << "\ni: " << i << "\n( ";
       Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->nu(i) - space_->nl(i) + 1);
       for( int j=0; j<nTemp; ++j ) {
         out << cSD_[i][j] <<", ";
         out << cVU_[i][j] <<"\t";
       }
       out << ")\n";
     }
     out << " --- velocity stencil: ---";
     for( int i=0; i<3; ++i ) {
       out << "\ni: " << i << "\n( ";
       Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->nu(i) - space_->nl(i) + 1);
       for( int j=0; j<nTemp; ++j ) {
         out << cVD_[i][j] <<", ";
         out << cVU_[i][j] <<"\t";
       }
       out << ")\n";
     }
   }


  bool hasApplyTranspose() const { return( false ); }


}; // end of class ConvectionSOp



/// \relates ConvectionSOp
template< class S=double, class O=int, int d=3 >
Teuchos::RCP<ConvectionSOp<S,O,d> > createConvectionSOp(
    const Teuchos::RCP<const Space<S,O,d> >& space ) {
  return( Teuchos::rcp( new ConvectionSOp<S,O,d>( space ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_CONVECTIONSOP_HPP
