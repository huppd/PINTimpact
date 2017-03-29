#pragma once
#ifndef PIMPACT_GRADOP_HPP
#define PIMPACT_GRADOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Stencil.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{



/// \ingroup BaseOperator
template<class ST>
class GradOp {

public:

  using SpaceT = ST;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

	static const int dimNC = ST::dimNC;
	static const int dim = ST::dimension;

	using StencD = Stencil< Scalar, Ordinal, 0, ST::SW::DL(0), ST::SW::DU(0) >;
	using StencG = Stencil< Scalar, Ordinal, 0, ST::SW::GL(0), ST::SW::GU(0) >;

  using TD = const Teuchos::Tuple<StencD,ST::sdim>; 
  using TG = const Teuchos::Tuple<StencG,ST::sdim>; 

  Teuchos::RCP<const SpaceT> space_;

  TG c_;
  TD cT_;

public:

  using DomainFieldT =  ScalarField<SpaceT>;
  using RangeFieldT =  VectorField<SpaceT>;

	/// \todo for the transposed stencil it would be better to compute locally all global stencils instead of MPI_allreduce them 
	/// \todo make it MG readey (participant, all reduce)
	GradOp( const Teuchos::RCP< const SpaceT>& space):
		space_(space) {

			//const bool mapping=true; // order ~2
			const bool mapping=false; // order ~6

			for( int dir=0; dir<ST::sdim; ++dir ) {
				// Gradient stencil
				
				F fdir = static_cast<F>( dir );

				c_[dir] = StencG( space_->nLoc(dir) );

				FD_getDiffCoeff(
						0,
						space_->nLoc(dir),
						space_->bl(dir),
						space_->bu(dir),
						space_->gl(dir),
						space_->gu(dir),
						space_->getBCLocal()->getBCL(dir),
						space_->getBCLocal()->getBCU(dir),
						space_->getShift(dir),
						2,
						dir+1,
						1,
						0,
						mapping,
						space_->getStencilWidths()->getDimNcbG(dir),
						space_->getStencilWidths()->getNcbG(dir),
						space_->getCoordinatesLocal()->getX( F::S, dir ),
						space_->getCoordinatesLocal()->getX( fdir, dir ),
						c_[dir].get() );

				// transposed Gradient stencil
				cT_[dir] = StencD( space_->nLoc(dir) );

				Ordinal nTempG = ( space_->nGlo(dir) + space_->bu(dir) - space_->bl(dir) + 1 )
					*( space_->bu(dir) - space_->bl(dir) + 1);


				Stencil< Scalar, Ordinal, ST::SW::BL(0), ST::SW::BL(0), ST::SW::BU(0) >
				cG1( space_->nGlo(dir) + space_->bu(dir) );
				Stencil< Scalar, Ordinal, ST::SW::BL(0), ST::SW::BL(0), ST::SW::BU(0) >
				cG2( space_->nGlo(dir) + space_->bu(dir) );

				for( Ordinal i = space_->si(fdir,dir,B::Y); i<=space_->ei(fdir,dir,B::Y); ++i )
					for( Ordinal ii = space_->gl(dir); ii<=space_->gu(dir); ++ii )
						cG1( i+space_->getShift(dir), ii ) = getC( static_cast<ECoord>(dir), i, ii );

				MPI_Allreduce(
						cG1.get(),	                              // const void *sendbuf,
						cG2.get(),                                // void *recvbuf,
						nTempG,			                              // int count,
						MPI_REAL8,	                              // MPI_Datatype datatype,
						MPI_SUM,		                              // MPI_Op op,
						space_->getProcGrid()->getCommBar(dir) ); // MPI_Comm comm )

				if( -1==space_->getBCGlobal()->getBCL(dir) ) {

					Ordinal ls1 = space_->getStencilWidths()->getLS(dir);
					Ordinal M1 = space_->nGlo(dir);

					for( Ordinal i=space->bl(dir); i<=-1; ++i )
						for( Ordinal ii=space->bl(dir); ii<=space->bu(dir); ++ii )
							cG2(2+ls1+i,ii) = cG2(M1+1+ls1+i,ii);

					for( Ordinal i=1; i<=space->bu(dir); ++i )
						for( Ordinal ii=space->bl(dir); ii<=space->bu(dir); ++ii )
							cG2(M1+ls1+i,ii) = cG2(1+ls1+i,ii);
				}

				for( Ordinal
						i =space_->si(F::S,static_cast<ECoord>(dir),B::Y);
						i<=space_->ei(F::S,static_cast<ECoord>(dir),B::Y);
						++i ) 
					for( Ordinal ii=space->dl(dir); ii<=space->du(dir); ++ii ) {
						cT_[dir](i,ii) = cG2( i+ii+space_->getShift(dir), -ii );
					}
			}
		};


  void apply( const DomainFieldT& x, RangeFieldT& y, const Add& add=Add::N ) const {

		applyG( x, y, add );

		if( Add::N==add )
			applyJ( y );
		//else
			//y.setBZero();
  }


  void applyG( const DomainFieldT& x, RangeFieldT& y, const Add& add=Add::N ) const {

		const B& b = ( (Add::N==add) ? B::Y : B::N );

		x.exchange(X);
		for( Ordinal k=space()->si(F::U,Z,b); k<=space()->ei(F::U,Z,b); ++k )
			for( Ordinal j=space()->si(F::U,Y,b); j<=space()->ei(F::U,Y,b); ++j )
				for( Ordinal i=space()->si(F::U,X,b); i<=space()->ei(F::U,X,b); ++i ) {
					if( Add::N==add ) y(F::U)(i,j,k) = 0.;
					y(F::U)(i,j,k) += innerStencU( x, i, j, k );
				}

		x.exchange(Y);
		for( Ordinal k=space()->si(F::V,Z,b); k<=space()->ei(F::V,Z,b); ++k )
			for( Ordinal j=space()->si(F::V,Y,b); j<=space()->ei(F::V,Y,b); ++j )
				for( Ordinal i=space()->si(F::V,X,b); i<=space()->ei(F::V,X,b); ++i ) {
					if( Add::N==add ) y(F::V)(i,j,k) = 0.;
					y(F::V)(i,j,k) += innerStencV( x, i, j, k );
				}

		if( 3==SpaceT::sdim )  {

			x.exchange(Z);
			for( Ordinal k=space()->si(F::W,Z,b); k<=space()->ei(F::W,Z,b); ++k )
				for( Ordinal j=space()->si(F::W,Y,b); j<=space()->ei(F::W,Y,b); ++j )
					for( Ordinal i=space()->si(F::W,X,b); i<=space()->ei(F::W,X,b); ++i ) {
						if( Add::N==add ) y(F::W)(i,j,k) = 0.;
						y(F::W)(i,j,k) += innerStencW( x, i, j, k );
					}
		}
  }


  void applyJ( RangeFieldT& y ) const {

		// BC scaling 
		const Scalar& eps = 0.1;
		
		for( F dir=F::U; dir<SpaceT::sdim; ++dir ) {
			B bc2 = B::Y;
			if( F::U!=dir ) {
				if( space()->getBCLocal()->getBCL(X) > 0 ) {
					Ordinal i = space()->si(dir,X,B::Y);
					for( Ordinal k=space()->si(dir,Z, bc2); k<=space()->ei(dir,Z,bc2); ++k )
						for( Ordinal j=space()->si(dir,Y,bc2); j<=space()->ei(dir,Y,bc2); ++j )
							y(dir)(i,j,k) *= eps;  
				}
				if( space()->getBCLocal()->getBCU(X) > 0 ) {
					Ordinal i = space()->ei(dir,X,B::Y);
					for( Ordinal k=space()->si(dir,Z,bc2); k<=space()->ei(dir,Z,bc2); ++k )
						for( Ordinal j=space()->si(dir,Y,bc2); j<=space()->ei(dir,Y,bc2); ++j )
							y(dir)(i,j,k) *= eps;  
				}
				bc2 = B::N;
			}

			if( F::V!=dir ) {
				if( space()->getBCLocal()->getBCL(Y) > 0 ) {
					Ordinal j = space()->si(dir,Y,B::Y);
					for( Ordinal k=space()->si(dir,Z,bc2); k<=space()->ei(dir,Z,bc2); ++k )
						for( Ordinal i=space()->si(dir,X,bc2); i<=space()->ei(dir,X,bc2); ++i ) 
							y(dir)(i,j,k) *= eps;  
				}
				if( space()->getBCLocal()->getBCU(Y) > 0 ) {
					Ordinal j = space()->ei(dir,Y,B::Y);
					for( Ordinal k=space()->si(dir,Z,bc2); k<=space()->ei(dir,Z,bc2); ++k )
						for( Ordinal i=space()->si(dir,X,bc2); i<=space()->ei(dir,X,bc2); ++i )
							y(dir)(i,j,k) *= eps;  
				}
				bc2 = B::N;
			}

			if( F::W!=dir ) {
				if( space()->getBCLocal()->getBCL(Z) > 0 ) {
					Ordinal k = space()->si(dir,Z,B::Y);
					for( Ordinal j=space()->si(dir,Y,bc2); j<=space()->ei(dir,Y,bc2); ++j )
						for( Ordinal i=space()->si(dir,X,bc2); i<=space()->ei(dir,X,bc2); ++i )
							y(dir)(i,j,k) *= eps;  
				}
				if( space()->getBCLocal()->getBCU(Z) > 0 ) {
					Ordinal k = space()->ei(dir,Z,B::Y);
					for( Ordinal j=space()->si(dir,Y,bc2); j<=space()->ei(dir,Y,bc2); ++j )
						for( Ordinal i=space()->si(dir,X,bc2); i<=space()->ei(dir,X,bc2); ++i )
							y(dir)(i,j,k) *= eps;  
				}
				bc2 = B::N;
			}
		}

		y.extrapolateBC();
		
    y.changed();
  }


  void apply( const RangeFieldT& x, DomainFieldT& y, const Add& add=Add::N ) const {

		for( int dir=0; dir<SpaceT::sdim; ++dir )
			x.exchange( dir, dir );

		if( 3==SpaceT::sdim )  {

			for( Ordinal k=space()->si(F::S,Z); k<=space()->ei(F::S,Z); ++k )
				for( Ordinal j=space()->si(F::S,Y); j<=space()->ei(F::S,Y); ++j )
					for( Ordinal i=space()->si(F::S,X); i<=space()->ei(F::S,X); ++i ) {
						if( Add::N==add ) y(i,j,k) = 0.;
						y(i,j,k) += innerStenc3D( x, i, j, k );
					}
		}
		else{

			for( Ordinal k=space()->si(F::S,Z); k<=space()->ei(F::S,Z); ++k )
				for( Ordinal j=space()->si(F::S,Y); j<=space()->ei(F::S,Y); ++j )
					for( Ordinal i=space()->si(F::S,X); i<=space()->ei(F::S,X); ++i ) {
						if( Add::N==add ) y(i,j,k) = 0.;
						y(i,j,k) += innerStenc2D( x, i, j, k );
					}
		}

    y.changed();
	}



  void assignField( const RangeFieldT& mv ) {};
  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(space_); };

	constexpr const Scalar* getC( const ECoord& dir ) const {
		return( c_[dir].get() );
	}

	constexpr const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const {
		return( c_[dir](i,off) );
	}

	constexpr const Scalar& getCTrans( const ECoord& dir, Ordinal i, Ordinal off ) const {
		return( cT_[dir](i,off) );
	}

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "\n--- " << getLabel() << " ---\n";
    //out << " --- stencil: ---";
    for( int dir=0; dir<ST::sdim; ++dir ) {
			out << "\ndir: " << static_cast<ECoord>(dir) << "\n";

			c_[dir].print( out );
    }

		out << "--- " << getLabel() << "^T ---\n";
		out << " --- stencil: ---";
		for( int dir=0; dir<ST::sdim; ++dir ) {
			out << "\ndir: " << static_cast<ECoord>(dir) << "\n";

			cT_[dir].print( out );
		}
  }

	const std::string getLabel() const { return( "Grad" ); };

protected:

	constexpr Scalar innerStencU( const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar grad = 0.;

		for( int ii=ST::SW::GL(X); ii<=ST::SW::GU(X); ++ii ) 
			grad += getC(X,i,ii)*x(i+ii,j,k);

		return( grad );
	}

	constexpr Scalar innerStencV( const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar grad = 0.;

		for( int jj=ST::SW::GL(Y); jj<=ST::SW::GU(Y); ++jj ) 
			grad += getC(Y,j,jj)*x(i,j+jj,k);

		return( grad );
	}

	constexpr Scalar innerStencW( const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar grad = 0.;

		for( int kk=ST::SW::GL(Z); kk<=ST::SW::GU(Z); ++kk ) 
			grad += getC(Z,k,kk)*x(i,j,k+kk);

		return( grad );
	}

	constexpr Scalar innerStenc3D( const RangeFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar gradT = 0.;

		for( int ii=ST::SW::DL(X); ii<=ST::SW::DU(X); ++ii ) 
			gradT += getCTrans(X,i,ii)*x(F::U)(i+ii,j,k);

		for( int jj=ST::SW::DL(Y); jj<=ST::SW::DU(Y); ++jj ) 
			gradT += getCTrans(Y,j,jj)*x(F::V)(i,j+jj,k);

		for( int kk=ST::SW::DL(Z); kk<=ST::SW::DU(Z); ++kk ) 
			gradT += getCTrans(Z,k,kk)*x(F::W)(i,j,k+kk);

		return( gradT );
	}

	constexpr Scalar innerStenc2D( const RangeFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar gradT = 0.;

		for( int ii=ST::SW::DL(X); ii<=ST::SW::DU(X); ++ii ) 
			gradT += getCTrans(X,i,ii)*x(F::U)(i+ii,j,k);

		for( int jj=ST::SW::DL(Y); jj<=ST::SW::DU(Y); ++jj ) 
			gradT += getCTrans(Y,j,jj)*x(F::V)(i,j+jj,k);

		return( gradT );
	}


}; // end of class GradOp


} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::GradOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::GradOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::GradOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::GradOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_GRADOP_HPP
