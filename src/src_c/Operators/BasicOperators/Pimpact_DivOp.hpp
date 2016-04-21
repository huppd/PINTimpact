#pragma once
#ifndef PIMPACT_DIVOP_HPP
#define PIMPACT_DIVOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{



extern "C"
void OP_div(
		const int& dimens,
		const int* const N,
		const int* const bl,
		const int* const bu,
		const int* const dl,
		const int* const du,
		const int* const ss,
		const int* const nn,
		const double* const c1,
		const double* const c2,
		const double* const c3,
		const double* const phiU,
		double* const phi );



/// \brief Divergence operator.
/// \ingroup BaseOperator
template<class ST>
class DivOp {

public:

  using SpaceT = ST;

  using DomainFieldT = VectorField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using TO = const Teuchos::Tuple<Scalar*,3>;

  Teuchos::RCP<const SpaceT> space_;

  TO c_;

public:

  DivOp( const Teuchos::RCP<const SpaceT>& space ):
    space_(space) {

    for( int i=0; i<3; ++i ) {
      Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->du(i) - space_->dl(i) + 1);

      c_[i] = new Scalar[ nTemp ];
      if( i<space_->dim() )
        FD_getDiffCoeff(
//            space_->rankST(),
            space_->nLoc(i),
            space_->bl(i),
            space_->bu(i),
            space_->dl(i),
            space_->du(i),
            space_->getBCLocal()->getBCL(i),
            space_->getBCLocal()->getBCU(i),
            space_->getShift(i),
            3,
            i+1,
            1,
            0,
						//true, // mapping
						false,
            space_->getStencilWidths()->getDimNcbD(i),
            space_->getStencilWidths()->getNcbD(i),
            space_->getCoordinatesLocal()->getX( i, i ),
            space_->getCoordinatesLocal()->getX( i, EField::S ),
            c_[i] );
    }
  };


	~DivOp() {
		for( int i=0; i<3; ++i )
			delete[] c_[i];
	}


  void apply( const DomainFieldT& x, RangeFieldT& y,
      Belos::ETrans trans=Belos::NOTRANS ) const {

    for( int dir=0; dir<space_->dim(); ++dir )
      x.exchange( dir, dir );

    OP_div(
        space_->dim(),
        space_->nLoc(),
        space_->bl(),
        space_->bu(),
        space_->dl(),
        space_->du(),
        space_->sInd(S),
        space_->eInd(S),
        getC(X),
        getC(Y),
        getC(Z),
        x.getConstRawPtr(),
        y.getRawPtr() );

    y.changed();
  }

  void assignField( const RangeFieldT& mv ) const {};
  void assignField( const DomainFieldT& mv ) const {};

  bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(space_); };

  inline constexpr const Scalar* getC( const ECoord& dir ) const {
      return( c_[dir] );
  }

	//inline constexpr const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const  {
		//return( c_[dir][ off - space_->dl(i) + 1 + (i-1)*( space_->du(i) - space_->dl(i) + 1) ] );
  //}

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    out << " --- stencil: ---";
    for( int i=0; i<3; ++i ) {
      out << "\ndir: " << i << "\n";
      Ordinal nTemp1 = ( space_->nLoc(i) + 1 );
      Ordinal nTemp2 = ( space_->du(i) - space_->dl(i) + 1);
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

	const std::string getLabel() const { return( "Div" ); };

};


} // end of namespace Pimpact



#ifdef COMPILE_ETI
extern template class Pimpact::DivOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::DivOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::DivOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::DivOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_DIVOP_HPP
