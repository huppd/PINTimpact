#pragma once
#ifndef PIMPACT_GRADOP_HPP
#define PIMPACT_GRADOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{



extern "C" {

void OP_grad(
    const int& dimens,
    const int* const N,
    const int* const bl,
    const int* const bu,
    const int* const gl,
    const int* const gu,
    const int* const su,
    const int* const nu,
    const int* const sv,
    const int* const nv,
    const int* const sw,
    const int* const nw,
    const double* const c1,
    const double* const c2,
    const double* const c3,
    const double* const phi,
    double* const grad );

void OP_SetBCZero(
    const int* const N,
    const int* const bl,
    const int* const bu,
    const int* const BCL,
    const int* const BCU,
    const int* const ss,
    const int* const nn,
    const double* phi );

void OP_extrapolateBC(
		const int& m,         
    const int* const N,         
    const int* const bL,
		const int* const bU,     
    const int& dL,
		const int& dU,     
		const int* const BC_L,
		const int* const BC_U, 
		const int* const SB,
		const int* const NB,
		const double* const c,    
		const double*       phi );

}



/// \ingroup BaseOperator
template<class ST>
class GradOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

protected:

  typedef const Teuchos::Tuple<Scalar*,3> TO;

  Teuchos::RCP<const SpaceT> space_;

  TO c_;

public:

  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef VectorField<SpaceT>  RangeFieldT;

  GradOp( const Teuchos::RCP< const SpaceT>& space):
    space_(space) {

    for( int i=0; i<3; ++i ) {
      Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->gu(i) - space_->gl(i) + 1);

      c_[i] = new Scalar[ nTemp ];

      if( i<space_->dim() )
        FD_getDiffCoeff(
//            space_->rankST(),
            space_->nLoc(i),
            space_->bl(i),
            space_->bu(i),
            space_->gl(i),
            space_->gu(i),
            space_->getBCLocal()->getBCL(i),
            space_->getBCLocal()->getBCU(i),
            space_->getShift(i),
            2,
            i+1,
            1,
            0,
            true,
            space_->getStencilWidths()->getDimNcbG(i),
            space_->getStencilWidths()->getNcbG(i),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            space_->getCoordinatesLocal()->getX( i, i ),
            c_[i] );
    }

  };


  ~GradOp() {
    for( int i=0; i<3; ++i ) {
      delete[] c_[i];
    }
  }



  void apply(const DomainFieldT& x, RangeFieldT& y) const {

		int dim = space_->dim();
    for( int i=0; i<dim; ++i) {
//			x.level();
      x.exchange(i);
    }

		OP_grad(
				space_->dim(),
				space_->nLoc(),
				space_->bl(),
				space_->bu(),
				space_->gl(),
				space_->gu(),
				space_->sInd(U),
				space_->eInd(U),
				space_->sInd(V),
				space_->eInd(V),
				space_->sInd(W),
				space_->eInd(W),
        getC(X),
        getC(Y),
        getC(Z),
				x.getConstRawPtr(),
				y.getRawPtr() );

		// necessary?
		for( int i=0; i<space()->dim(); ++i ) {
			OP_SetBCZero(
					space_->nLoc(),
					space_->bl(),
					space_->bu(),
					space_->getBCLocal()->getBCL(),
					space_->getBCLocal()->getBCU(),
					space_->sIndB(i),
					space_->eIndB(i),
					y.getRawPtr(i) );
			OP_extrapolateBC(
					i+1,
					space_->nLoc(),
					space_->bl(),
					space_->bu(),
					space_->dl(i),
					space_->du(i),
					space_->getBCLocal()->getBCL(),
					space_->getBCLocal()->getBCU(),
					space_->sIndB(i),
					space_->eIndB(i),
					space_->getInterpolateV2S()->getC(
						static_cast<ECoord>(i) ),
					y.getRawPtr(i) );
		}

    y.changed();
  }

  void assignField( const RangeFieldT& mv ) {};
  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(space_); };

  const Scalar* getC( const ECoord& dir ) const {
      return( c_[dir] );
  }

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    out << " --- stencil: ---";
    for( int i=0; i<3; ++i ) {
      out << "\ndir: " << i << "\n";
      Ordinal nTemp1 = ( space_->nLoc(i) + 1 );
      Ordinal nTemp2 = ( space_->gu(i) - space_->gl(i) + 1 );
      for( int j=0; j<nTemp1; ++j ) {
        out << "\ni: " << j << "\t(";
        for( int k=0; k<nTemp2; ++k ) {
          out << c_[i][k+nTemp2*j] <<", ";
        }
        out << ")\n";
      }
      out << "\n";
    }
  }

	const std::string getLabel() const { return( "Grad" ); };

}; // end of class GradOp


} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::GradOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::GradOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::GradOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::GradOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_GRADOP_HPP
