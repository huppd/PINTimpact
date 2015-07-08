#pragma once
#ifndef PIMPACT_HELMHOLTZOP_HPP
#define PIMPACT_HELMHOLTDOP_HPP

#include "Pimpact_extern_FDCoeff.hpp"

#include "Pimpact_Types.hpp"

#include "Pimpact_VectorField.hpp"



namespace Pimpact{


extern "C" {
void OP_helmholtz(
    const int& dimens,
    const int* const N,
    const int* const bl,
    const int* const bu,
    const int* const ss,
    const int* const nn,
    const double* const c11,
    const double* const c22,
    const double* const c33,
    const double& mulI,
    const double& multL,
    const double* const phi,
    double* const lap );
}



/// \brief Helmholtz operator
/// \ingroup BaseOperator
///
/// computes \f$ y = ( mulI_ I - mulL_ \Delta) x \f$
template<class ST>
class HelmholtzOp {


public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef const Teuchos::Tuple<Scalar*,3> TO;

protected:

  Teuchos::RCP<const SpaceT> space_;

	Scalar mulI_;
	Scalar mulL_;
	
  TO cS_;
  TO cV_;

public:

  typedef VectorField<SpaceT>  DomainFieldT;
  typedef VectorField<SpaceT>  RangeFieldT;

	HelmholtzOp(
			const Teuchos::RCP<const SpaceT>& space ):
		space_(space),
		mulI_( (Scalar)0. ),
		mulL_( 1./space_->getDomain()->getDomainSize()->getRe() )	{

    for( int i=0; i<3; ++i ) {
      Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->bu(i) - space_->bl(i) + 1);

      cS_[i] = new Scalar[ nTemp ];
      if( i<space_->dim() )
        FD_getDiffCoeff(
            space_->rankST(),
            space_->nLoc(i),
            space_->bl(i),
            space_->bu(i),
            space_->bl(i),
            space_->bu(i),
            space_->getDomain()->getBCLocal()->getBCL(i),
            space_->getDomain()->getBCLocal()->getBCU(i),
            space_->getShift(i),
            int(EField::S)+1,
            i+1,
            2,
            0,
            true,
            space_->getStencilWidths()->getDimNcbC(i),
            space_->getStencilWidths()->getNcbC(i),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            cS_[i] );

      cV_[i] = new Scalar[ nTemp ];
      if( i<space_->dim() )
        FD_getDiffCoeff(
            space_->rankST(),
            space_->nLoc(i),
            space_->bl(i),
            space_->bu(i),
            space_->bl(i),
            space_->bu(i),
            space_->getDomain()->getBCLocal()->getBCL(i),
            space_->getDomain()->getBCLocal()->getBCU(i),
            space_->getShift(i),
            1,
            i+1,
            2,
            0,
            true,
            space_->getStencilWidths()->getDimNcbC(i),
            space_->getStencilWidths()->getNcbC(i),
            space_->getCoordinatesLocal()->getX( i, i ),
            space_->getCoordinatesLocal()->getX( i, i ),
            cV_[i] );
    }

  };

  ~HelmholtzOp() {
    for( int i=0; i<3; ++i ) {
      delete[] cS_[i];
      delete[] cV_[i];
    }
  }


  void apply(const DomainFieldT& x, RangeFieldT& y ) const {

    x.exchange();

    for( int i=0; i<space_->dim(); ++i ) {
      EField fType = (EField)i;
      OP_helmholtz(
          space_->dim(),
          space_->nLoc(),
          space_->bl(),
          space_->bu(),
          space_->sInd(fType),
          space_->eInd(fType),
          getC(X,fType),
          getC(Y,fType),
          getC(Z,fType),
          mulI_,
          mulL_,
          x.getConstRawPtr(fType),
          y.getRawPtr(fType) );
    }

    y.changed();

  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

  void print( std::ostream& out=std::cout ) const {

    out << "--- " << getLabel() << " ---\n";
    out << " --- scalar stencil: ---";

    for( int i=0; i<3; ++i ) {

      out << "\ndir: " << i << "\n";

      Ordinal nTemp1 = ( space_->nLoc(i) + 1 );
      Ordinal nTemp2 = ( space_->bu(i) - space_->bl(i) + 1 );

      for( int j=0; j<nTemp1; ++j ) {
        out << "\ni: " << j << "\t(";
        for( int k=0; k<nTemp2; ++k ) {
          out << cS_[i][k+nTemp2*j] <<", ";
        }
        out << ")\n";
      }
      out << "\n";
    }
    out << " --- velocity stencil: ---";
    for( int i=0; i<3; ++i ) {
      out << "\ndir: " << i << "\n";
      Ordinal nTemp1 = ( space_->nLoc(i) + 1 );
      Ordinal nTemp2 = ( space_->bu(i) - space_->bl(i) + 1 );
      for( int j=0; j<nTemp1; ++j ) {
        out << "\ni: " << j << "\t(";
        for( int k=0; k<nTemp2; ++k ) {
          out << cV_[i][k+nTemp2*j] <<", ";
        }
        out << ")\n";
      }
      out << "\n";
    }
  }


	Teuchos::RCP<const SpaceT> space() const { return(space_); };

  const Scalar* getC( const ECoord& dir, const EField& ftype ) const {
      return( ((int)dir==(int)ftype)?cV_[dir]:cS_[dir] );
  }

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		mulI_ = para->get<Scalar>( "mulI", 0. );
		mulL_ = para->get<Scalar>( "mulL", 1./space_->getDomain()->getDomainSize()->getRe() );
	}

	const std::string getLabel() const { return( "Helmholtz" ); };

}; // end of class HelmholtzOp



} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::HelmholtzOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::HelmholtzOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::HelmholtzOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::HelmholtzOp< Pimpact::Space<double,int,4,4> >;
#endif
	

#endif // end of #ifndef PIMPACT_HELMHOLTZOP_HPP
