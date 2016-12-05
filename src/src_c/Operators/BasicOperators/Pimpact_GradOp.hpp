#pragma once
#ifndef PIMPACT_GRADOP_HPP
#define PIMPACT_GRADOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Stencil.hpp"
#include "Pimpact_Types.hpp"
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

	using SW = StencilWidths<dim,dimNC>;

	using StencD = Stencil< Scalar, Ordinal, 0, SW::DL(0), SW::DU(0) >;
	using StencG = Stencil< Scalar, Ordinal, 0, SW::GL(0), SW::GU(0) >;

  using TD = const Teuchos::Tuple<StencD*,ST::sdim>; 
  using TG = const Teuchos::Tuple<StencG*,ST::sdim>; 

  Teuchos::RCP<const SpaceT> space_;

  TG c_;
  TD cT_;

public:

  using DomainFieldT =  ScalarField<SpaceT>;
  using RangeFieldT =  VectorField<SpaceT>;

	/// \todo for the transposed stencil it would be better to compute locally all global stencils instead of MPI_allreduce them 
	GradOp( const Teuchos::RCP< const SpaceT>& space):
		space_(space) {

			for( int dir=0; dir<ST::sdim; ++dir ) {
				// Gradient stencil
				
				c_[dir] = new StencG( space_->nLoc(dir) );

				FD_getDiffCoeff(
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
						true, // mapping
						//false,
						space_->getStencilWidths()->getDimNcbG(dir),
						space_->getStencilWidths()->getNcbG(dir),
						space_->getCoordinatesLocal()->getX( dir, EField::S ),
						space_->getCoordinatesLocal()->getX( dir, dir ),
						c_[dir]->get() );

				// transposed Gradient stencil
				cT_[dir] = new StencD( space_->nLoc(dir) );

				Ordinal nTempG = ( space_->nGlo(dir) + space_->bu(dir) - space_->bl(dir) + 1 )
					*( space_->bu(dir) - space_->bl(dir) + 1);


				Stencil< Scalar, Ordinal, SW::BL(0), SW::BL(0), SW::BU(0) >
				cG1( space_->nGlo(dir) + space_->bu(dir) );
				Stencil< Scalar, Ordinal, SW::BL(0), SW::BL(0), SW::BU(0) >
				cG2( space_->nGlo(dir) + space_->bu(dir) );

				for( Ordinal i = space_->begin(dir,dir,With::B); i<=space_->end(dir,dir,With::B); ++i )
					for( Ordinal ii = space_->gl(dir); ii<=space_->gu(dir); ++ii )
						cG1( i+space_->getShift(dir), ii ) = getC( static_cast<ECoord>(dir), i, ii );

				MPI_Allreduce(
						cG1.get(),	                                // const void *sendbuf,
						cG2.get(),                                  // void *recvbuf,
						nTempG,			                                // int count,
						MPI_REAL8,	                                // MPI_Datatype datatype,
						MPI_SUM,		                                // MPI_Op op,
						space_->getProcGrid()->getCommSlice(dir) ); // MPI_Comm comm )

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
						i =space_->begin(S,static_cast<ECoord>(dir),With::B);
						i<=space_->end(S,static_cast<ECoord>(dir),With::B);
						++i ) 
					for( Ordinal ii=space->dl(dir); ii<=space->du(dir); ++ii ) {
						cT_[dir]->operator()(i,ii) = cG2( i+ii+space_->getShift(dir), -ii );
					}
			}
		};


  ~GradOp() {
    for( int i=0; i<ST::sdim; ++i ) {
      delete c_[i];
      delete cT_[i];
		}
  }


  void apply( const DomainFieldT& x, RangeFieldT& y, const With& withB=With::B ) const {

		applyG( x,y, withB );

		if( With::B==withB )
			applyJ( y );
		else
			y.setBZero();
  }


  void applyG( const DomainFieldT& x, RangeFieldT& y, const With& withB=With::B ) const {

		x.exchange(X);
		for( Ordinal k=space()->begin(U,Z,withB); k<=space()->end(U,Z,withB); ++k )
			for( Ordinal j=space()->begin(U,Y,withB); j<=space()->end(U,Y,withB); ++j )
				for( Ordinal i=space()->begin(U,X,withB); i<=space()->end(U,X,withB); ++i )
					y.getField(U).at(i,j,k) = innerStencU( x, i, j, k );

		x.exchange(Y);
		for( Ordinal k=space()->begin(V,Z,withB); k<=space()->end(V,Z,withB); ++k )
			for( Ordinal j=space()->begin(V,Y,withB); j<=space()->end(V,Y,withB); ++j )
				for( Ordinal i=space()->begin(V,X,withB); i<=space()->end(V,X,withB); ++i )
					y.getField(V).at(i,j,k) = innerStencV( x, i, j, k );

		if( 3==SpaceT::sdim )  {

			x.exchange(Z);
			for( Ordinal k=space()->begin(W,Z,withB); k<=space()->end(W,Z,withB); ++k )
				for( Ordinal j=space()->begin(W,Y,withB); j<=space()->end(W,Y,withB); ++j )
					for( Ordinal i=space()->begin(W,X,withB); i<=space()->end(W,X,withB); ++i )
						y.getField(W).at(i,j,k) = innerStencW( x, i, j, k );
		}
  }


  void applyJ( RangeFieldT& y ) const {

		// BC scaling 
		const Scalar& eps = 0.1;
		//const Scalar& eps = 1.;
		
		for( int dir=0; dir<SpaceT::sdim; ++dir ) {
			With bc2 = With::B;
			if( 0!=dir ) {
				if( space()->getBCLocal()->getBCL(X) > 0 ) {
					Ordinal i = space()->begin(dir,X,With::B);
					for( Ordinal k=space()->begin(dir,Z, bc2); k<=space()->end(dir,Z,bc2); ++k )
						for( Ordinal j=space()->begin(dir,Y,bc2); j<=space()->end(dir,Y,bc2); ++j )
							y.getField(dir).at(i,j,k) *= eps;  
				}
				if( space()->getBCLocal()->getBCU(X) > 0 ) {
					Ordinal i = space()->end(dir,X,With::B);
					for( Ordinal k=space()->begin(dir,Z,bc2); k<=space()->end(dir,Z,bc2); ++k )
						for( Ordinal j=space()->begin(dir,Y,bc2); j<=space()->end(dir,Y,bc2); ++j )
							y.getField(dir).at(i,j,k) *= eps;  
				}
				bc2 = With::noB;
			}

			if( 1!=dir ) {
				if( space()->getBCLocal()->getBCL(Y) > 0 ) {
					Ordinal j = space()->begin(dir,Y,With::B);
					for( Ordinal k=space()->begin(dir,Z,bc2); k<=space()->end(dir,Z,bc2); ++k )
						for( Ordinal i=space()->begin(dir,X,bc2); i<=space()->end(dir,X,bc2); ++i ) 
							y.getField(dir).at(i,j,k) *= eps;  
				}
				if( space()->getBCLocal()->getBCU(Y) > 0 ) {
					Ordinal j = space()->end(dir,Y,With::B);
					for( Ordinal k=space()->begin(dir,Z,bc2); k<=space()->end(dir,Z,bc2); ++k )
						for( Ordinal i=space()->begin(dir,X,bc2); i<=space()->end(dir,X,bc2); ++i )
							y.getField(dir).at(i,j,k) *= eps;  
				}
				bc2 = With::noB;
			}

			if( 2!=dir ) {
				if( space()->getBCLocal()->getBCL(Z) > 0 ) {
					Ordinal k = space()->begin(dir,Z,With::B);
					for( Ordinal j=space()->begin(dir,Y,bc2); j<=space()->end(dir,Y,bc2); ++j )
						for( Ordinal i=space()->begin(dir,X,bc2); i<=space()->end(dir,X,bc2); ++i )
							y.getField(dir).at(i,j,k) *= eps;  
				}
				if( space()->getBCLocal()->getBCU(Z) > 0 ) {
					Ordinal k = space()->end(dir,Z,With::B);
					for( Ordinal j=space()->begin(dir,Y,bc2); j<=space()->end(dir,Y,bc2); ++j )
						for( Ordinal i=space()->begin(dir,X,bc2); i<=space()->end(dir,X,bc2); ++i )
							y.getField(dir).at(i,j,k) *= eps;  
				}
				bc2 = With::noB;
			}
		}

		y.extrapolateBC();
		
    y.changed();
  }


  void apply(const RangeFieldT& x, DomainFieldT& y) const {

		for( int dir=0; dir<SpaceT::sdim; ++dir )
			x.exchange( dir, dir );

		if( 3==SpaceT::sdim )  {

			for( Ordinal k=space()->begin(S,Z); k<=space()->end(S,Z); ++k )
				for( Ordinal j=space()->begin(S,Y); j<=space()->end(S,Y); ++j )
					for( Ordinal i=space()->begin(S,X); i<=space()->end(S,X); ++i )
						y.at(i,j,k) = innerStenc3D( x, i, j, k );
		}
		else{

			for( Ordinal k=space()->begin(S,Z); k<=space()->end(S,Z); ++k )
				for( Ordinal j=space()->begin(S,Y); j<=space()->end(S,Y); ++j )
					for( Ordinal i=space()->begin(S,X); i<=space()->end(S,X); ++i )
						y.at(i,j,k) = innerStenc2D( x, i, j, k );
		}

    y.changed();
	}



  void assignField( const RangeFieldT& mv ) {};
  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(space_); };

	constexpr const Scalar* getC( const ECoord& dir ) const {
		return( c_[dir]->get() );
	}

	constexpr const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const {
		return( c_[dir]->at(i,off) );
	}

	constexpr const Scalar& getCTrans( const ECoord& dir, Ordinal i, Ordinal off ) const {
		return( cT_[dir]->at(i,off) );
	}

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "\n--- " << getLabel() << " ---\n";
    //out << " --- stencil: ---";
    for( int dir=0; dir<ST::sdim; ++dir ) {
			out << "\ndir: " << toString(static_cast<ECoord>(dir)) << "\n";

			c_[dir]->print( out );
    }

		out << "--- " << getLabel() << "^T ---\n";
		out << " --- stencil: ---";
		for( int dir=0; dir<ST::sdim; ++dir ) {
			out << "\ndir: " << toString(static_cast<ECoord>(dir)) << "\n";

			cT_[dir]->print( out );
		}
  }

	const std::string getLabel() const { return( "Grad" ); };

protected:

	inline constexpr Scalar innerStencU( const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar grad = 0.;

		for( int ii=space_->gl(X); ii<=space_->gu(X); ++ii ) 
			grad += getC(X,i,ii)*x.at(i+ii,j,k);

		return( grad );
	}

	inline constexpr Scalar innerStencV( const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar grad = 0.;

		for( int jj=space_->gl(Y); jj<=space_->gu(Y); ++jj ) 
			grad += getC(Y,j,jj)*x.at(i,j+jj,k);

		return( grad );
	}

	inline constexpr Scalar innerStencW( const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar grad = 0.;

		for( int kk=space_->gl(Z); kk<=space_->gu(Z); ++kk ) 
			grad += getC(Z,k,kk)*x.at(i,j,k+kk);

		return( grad );
	}

	inline constexpr Scalar innerStenc3D( const RangeFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar gradT = 0.;

		for( int ii=space_->dl(X); ii<=space_->du(X); ++ii ) 
			gradT += getCTrans(X,i,ii)*x.getField(U).at(i+ii,j,k);

		for( int jj=space_->dl(Y); jj<=space_->du(Y); ++jj ) 
			gradT += getCTrans(Y,j,jj)*x.getField(V).at(i,j+jj,k);

		for( int kk=space_->dl(Z); kk<=space_->du(Z); ++kk ) 
			gradT += getCTrans(Z,k,kk)*x.getField(W).at(i,j,k+kk);

		return( gradT );
	}

	inline constexpr Scalar innerStenc2D( const RangeFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar gradT = 0.;

		for( int ii=space_->dl(X); ii<=space_->du(X); ++ii ) 
			gradT += getCTrans(X,i,ii)*x.getField(U).at(i+ii,j,k);

		for( int jj=space_->dl(Y); jj<=space_->du(Y); ++jj ) 
			gradT += getCTrans(Y,j,jj)*x.getField(V).at(i,j+jj,k);

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
