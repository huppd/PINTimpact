#pragma once
#ifndef PIMPACT_CONVECTIONSOP_HPP
#define PIMPACT_CONVECTIONSOP_HPP


#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact {



extern "C"
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
    double* const nlu,
    const double& mul,
    const double& mulI,
    const double& mulC );



/// \brief convection operator, that takes the free interpolated velocity components and advects accordingly
/// \ingroup NonliearOperator
template<class ST>
class ConvectionSOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef Teuchos::Tuple< Teuchos::RCP< ScalarField<SpaceT> >, 3 > FluxFieldT;
  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;

protected:

  typedef const Teuchos::Tuple<Scalar*,3> TO;

  Teuchos::RCP<const SpaceT> space_;

  TO cSD_;
  TO cVD_;

  TO cSU_;
  TO cVU_;

public:

	ConvectionSOp( const Teuchos::RCP<const SpaceT>& space  ):
		space_(space) {

			for( int i=0; i<3; ++i ) {

				Ordinal nTemp = ( space_->nLoc(i) + 1 )*( space_->nu(i) - space_->nl(i) + 1);
				cSD_[i] = new Scalar[ nTemp ];

				if( i<space_->dim() )
					FD_getDiffCoeff(
							space_->nLoc(i),
							space_->bl(i),
							space_->bu(i),
							space_->nl(i),
							space_->nu(i),
							space_->getBCLocal()->getBCL(i),
							space_->getBCLocal()->getBCU(i),
							space_->getShift(i),
							int(EField::S)+1,
							i+1,
							1,
							-1,
							true,
							space_->getStencilWidths()->getDimNcbC(i),
							space_->getStencilWidths()->getNcbC(i),
							space_->getCoordinatesLocal()->getX( i, EField::S ),
							space_->getCoordinatesLocal()->getX( i, EField::S ),
							cSD_[i] );

				cSU_[i] = new Scalar[ nTemp ];

				if( i<space_->dim() )
					FD_getDiffCoeff(
							space_->nLoc(i),
							space_->bl(i),
							space_->bu(i),
							space_->nl(i),
							space_->nu(i),
							space_->getBCLocal()->getBCL(i),
							space_->getBCLocal()->getBCU(i),
							space_->getShift(i),
							int(EField::S)+1,
							i+1,
							1,
							+1,
							true,
							space_->getStencilWidths()->getDimNcbC(i),
							space_->getStencilWidths()->getNcbC(i),
							space_->getCoordinatesLocal()->getX( i, EField::S ),
							space_->getCoordinatesLocal()->getX( i, EField::S ),
							cSU_[i] );

				cVD_[i] = new Scalar[ nTemp ];

				if( i<space_->dim() )
					FD_getDiffCoeff(
							space_->nLoc(i),
							space_->bl(i),
							space_->bu(i),
							space_->nl(i),
							space_->nu(i),
							space_->getBCLocal()->getBCL(i),
							space_->getBCLocal()->getBCU(i),
							space_->getShift(i),
							1,
							i+1,
							1,
							-1,
							true,
							space_->getStencilWidths()->getDimNcbC(i),
							space_->getStencilWidths()->getNcbC(i),
							space_->getCoordinatesLocal()->getX( i, i ),
							space_->getCoordinatesLocal()->getX( i, i ),
							cVD_[i] );

				cVU_[i] = new Scalar[ nTemp ];
				if( i<space_->dim() )
					FD_getDiffCoeff(
							space_->nLoc(i),
							space_->bl(i),
							space_->bu(i),
							space_->nl(i),
							space_->nu(i),
							space_->getBCLocal()->getBCL(i),
							space_->getBCLocal()->getBCU(i),
							space_->getShift(i),
							1,
							i+1,
							1,
							+1,
							true,
							space_->getStencilWidths()->getDimNcbC(i),
							space_->getStencilWidths()->getNcbC(i),
							space_->getCoordinatesLocal()->getX( i, i ),
							space_->getCoordinatesLocal()->getX( i, i ),
							cVU_[i] );

				for( Ordinal j=0; j<nTemp; ++j ) {
					if( std::isnan( cSD_[i][j] ) || std::isinf( cSD_[i][j] ) ) cSD_[i][j] = 0.;
					if( std::isnan( cSU_[i][j] ) || std::isinf( cSU_[i][j] ) ) cSU_[i][j] = 0.;
					if( std::isnan( cVD_[i][j] ) || std::isinf( cVD_[i][j] ) ) cVD_[i][j] = 0.;
					if( std::isnan( cVU_[i][j] ) || std::isinf( cVU_[i][j] ) ) cVU_[i][j] = 0.;
				}

			}

	};


  void assignField( const RangeFieldT& mv ) {};


	void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z,
			Scalar mul, Scalar mulI, Scalar mulC, Scalar mulL ) const {
		apply( x, y, z, mul, mulC );
	}


	void apply( const FluxFieldT& x, const DomainFieldT& y, RangeFieldT& z, Scalar mul=0., Scalar mulC=1. ) const {

		int m = (int)z.getType();

		TEUCHOS_TEST_FOR_EXCEPTION(
				z.getType() != y.getType(),
				std::logic_error,
				"Pimpact::ConvectionSOP can only be applied to same fieldType !!!\n");


		for( int i =0; i<space_->dim(); ++i ) {
			TEUCHOS_TEST_FOR_EXCEPTION(
					x[i]->getType() != y.getType(),
					std::logic_error,
					"Pimpact::ConvectionSOP can only be applied to same fieldType !!!\n");
		}

		for( int vel_dir=0; vel_dir<space_->dim(); ++vel_dir )
			x[vel_dir]->exchange();

		y.exchange();

		OP_convection(
				space_->dim(),
				space_->nLoc(),
				space_->bl(),
				space_->bu(),
				space_->nl(),
				space_->nu(),
				space_->sInd(m),
				space_->eInd(m),
				getCD(X,z.getType()),
				getCD(Y,z.getType()),
				getCD(Z,z.getType()),
				getCU(X,z.getType()),
				getCU(Y,z.getType()),
				getCU(Z,z.getType()),
				x[0]->getConstRawPtr(),
				x[1]->getConstRawPtr(),
				x[2]->getConstRawPtr(),
				y.getConstRawPtr(),
				z.getRawPtr(),
				mul,
				0.,
				mulC );

		z.changed();

	}


  void print( std::ostream& out=std::cout ) const {
     out << " --- ConvectioSOp ---\n";
     for( int i=0; i<3; ++i ) {
       out << "dir: " << i << "\n ";
       out << "i\n";// cSD,\t\t cVD,\t\t cSU,\t\t cVU\n ";
       for( Ordinal j=0; j<( space_->nLoc(i) + 1 ); ++j ) {
         out << j << "\n";
         Ordinal km = ( space_->nu(i) - space_->nl(i) + 1);
         out << "cSD:\t";
         for( Ordinal k=0; k<km; ++k )
           out << getCD((ECoord)i,S)[k+j*km] <<",\t";
         out << "\ncVD:\t";
         for( Ordinal k=0; k<km; ++k )
           out << getCD((ECoord)i,(EField)i)[k+j*km] <<",\t";
         out << "\ncSU:\t";
         for( Ordinal k=0; k<km; ++k )
           out << getCU((ECoord)i,S)[k+j*km] <<",\t";
         out << "\ncVU:\t";
         for( Ordinal k=0; k<km; ++k )
           out << getCU((ECoord)i,(EField)i)[k+j*km] <<",\t";
         out << "\n";
       }
//       out << "i\t cSD,\t\t cVD,\t\t cSU,\t\t cVU\n ";
//       for( Ordinal j=0; j<( space_->nLoc(i) + 1 ); ++j ) {
//         out << j << "\t ";
//         Ordinal km = ( space_->nu(i) - space_->nl(i) + 1);
//         for( Ordinal k=0; k<km; ++k )
//           out << getCD((ECoord)i,S)[k+j*km] <<", ";
//         out << "\t";
//         for( Ordinal k=0; k<km; ++k )
//           out << getCD((ECoord)i,(EField)i)[k+j*km] <<", ";
//         out << "\t";
//         for( Ordinal k=0; k<km; ++k )
//           out << getCU((ECoord)i,S)[k+j*km] <<", ";
//         out << "\t";
//         for( Ordinal k=0; k<km; ++k )
//           out << getCU((ECoord)i,(EField)i)[k+j*km] <<", ";
//         out << "\n";
//       }
     }
  }


  bool hasApplyTranspose() const { return( false ); }

  Teuchos::RCP<const SpaceT>  space() const { return( space_ ); }

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}


  const Scalar* getCU( const ECoord& dir, const EField& ftype ) const  {
    return( ( ((int)dir)==((int)ftype) )?cVU_[dir]:cSU_[dir] );
//      return( cVU_[dir] );
  }

  const Scalar* getCD( const ECoord& dir, const EField& ftype ) const  {
    return( ( ((int)dir)==((int)ftype) )?cVD_[dir]:cSD_[dir] );
//      return( cVD_[dir] );
  }

	const std::string getLabel() const { return( "Convection" ); };


}; // end of class ConvectionSOp


} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::ConvectionSOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::ConvectionSOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::ConvectionSOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::ConvectionSOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_CONVECTIONSOP_HPP
