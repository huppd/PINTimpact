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



/// \ingroup BaseOperator
template<class ST>
class GradOp {

public:

  using SpaceT = ST;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using TO = const Teuchos::Tuple<Scalar*,3>; 
  Teuchos::RCP<const SpaceT> space_;

  TO c_;

  TO cT_;

public:

  using DomainFieldT =  ScalarField<SpaceT>;
  using RangeFieldT =  VectorField<SpaceT>;

	/// \todo for the transposed stencil it would be better to compute locally all global stencils instead of MPI_allreduce them 
	GradOp( const Teuchos::RCP< const SpaceT>& space):
		space_(space) {

			for( int dir=0; dir<3; ++dir ) {
				// Gradient stencil
				Ordinal nTemp = ( space_->nLoc(dir) + 1 )*( space_->gu(dir) - space_->gl(dir) + 1);

				c_[dir] = new Scalar[ nTemp ];

				if( dir<space_->dim() )
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
							//true, // mapping
							false,
							space_->getStencilWidths()->getDimNcbG(dir),
							space_->getStencilWidths()->getNcbG(dir),
							space_->getCoordinatesLocal()->getX( dir, EField::S ),
							space_->getCoordinatesLocal()->getX( dir, dir ),
							c_[dir] );

				// transposed Gradient stencil
				Ordinal nTempT  = ( space_->nLoc(dir) + 1 )*( space_->du(dir) - space_->dl(dir) + 1);
				cT_[dir] = new Scalar[ nTempT ];
				for( Ordinal i = 0; i<nTempT; ++i )
					cT_[dir][i] = 0.;

				Ordinal nTempG = ( space_->nGlo(dir) + space_->bu(dir) - space_->bl(dir) + 1 )
					*( space_->bu(dir) - space_->bl(dir) + 1);

				Scalar* cG1 = new Scalar[ nTempG ];
				Scalar* cG2 = new Scalar[ nTempG ];

				for( Ordinal i = 0; i<nTempG; ++i ) {
					cG1[i] = 0.;
					cG2[i] = 0.;
				}

				for( Ordinal i = space_->begin(dir,dir,With::B); i<=space_->end(dir,dir,With::B); ++i )
					for( Ordinal ii = space_->gl(dir); ii<=space_->gu(dir); ++ii ) {
						Ordinal ind = ( ii - space_->bl(dir) ) + ( i+space_->getShift(dir)-space_->bl(dir) )*( space_->bu(dir) - space_->bl(dir) + 1 );
						cG1[ ind ]= getC( static_cast<ECoord>(dir), i, ii );
					}

				MPI_Allreduce(
						cG1,    		                                // const void *sendbuf,
						cG2,    		                                // void *recvbuf,
						nTempG,			                                // int count,
						MPI_REAL8,	                                // MPI_Datatype datatype,
						MPI_SUM,		                                // MPI_Op op,
						space_->getProcGrid()->getCommSlice(dir) ); // MPI_Comm comm )

				if( -1==space_->getBCGlobal()->getBCL(dir) ) {

					Ordinal ls1 = space_->getStencilWidths()->getLS(dir);
					Ordinal M1 = space_->nGlo(dir);

					for( Ordinal i=space->bl(dir); i<=-1; ++i )
						for( Ordinal ii=space->bl(dir); ii<=space->bu(dir); ++ii ) {
							Ordinal indT = ii-space_->bl(dir) + (2+ls1+i-space_->bl(dir)   )*( space_->bu(dir) - space_->bl(dir) + 1 );
							Ordinal indS = ii-space_->bl(dir) + (M1+1+ls1+i-space_->bl(dir))*( space_->bu(dir) - space_->bl(dir) + 1 );
							cG2[ indT ] = cG2[ indS ];
						}

					for( Ordinal i=1; i<=space->bu(dir); ++i )
						for( Ordinal ii=space->bl(dir); ii<=space->bu(dir); ++ii ) {
							Ordinal indT = ii-space_->bl(dir) + (M1+ls1+i-space_->bl(dir))*( space_->bu(dir) - space_->bl(dir) + 1 );
							Ordinal indS = ii-space_->bl(dir) + ( 1+ls1+i-space_->bl(dir))*( space_->bu(dir) - space_->bl(dir) + 1 );
							cG2[ indT ] = cG2[ indS ];
						}
				}

				for( Ordinal
						i =space_->begin(S,static_cast<ECoord>(dir),With::B);
						i<=space_->end(S,static_cast<ECoord>(dir),With::B);
						++i ) 
					for( Ordinal ii=space->dl(dir); ii<=space->du(dir); ++ii ) {
						Ordinal ind1 = ( ii - space_->dl(dir) ) + ( i )*( space_->du(dir) - space_->dl(dir) + 1 );
						Ordinal ind2 = ( -ii - space_->bl(dir) ) + ( i+ii+space_->getShift(dir)-space_->bl(dir) )*( space_->bu(dir) - space_->bl(dir) + 1 );
						cT_[dir][ind1] = cG2[ ind2 ];
					}

				delete[] cG1;
				delete[] cG2;
			}
		};


  ~GradOp() {
    for( int i=0; i<3; ++i )
      delete[] c_[i];
  }


  void apply( const DomainFieldT& x, RangeFieldT& y ) const {

		applyG( x,y );
		applyJ( y );
  }


  void applyG( const DomainFieldT& x, RangeFieldT& y ) const {

		With bcY_ = With::B;
		//With bcY_ = With::noB;

		x.exchange(X);
		for( Ordinal k=space()->begin(U,Z,bcY_); k<=space()->end(U,Z,bcY_); ++k )
			for( Ordinal j=space()->begin(U,Y,bcY_); j<=space()->end(U,Y,bcY_); ++j )
				for( Ordinal i=space()->begin(U,X,bcY_); i<=space()->end(U,X,bcY_); ++i )
					y.getField(U).at(i,j,k) = innerStencU( x, i, j, k );

		x.exchange(Y);
		for( Ordinal k=space()->begin(V,Z,bcY_); k<=space()->end(V,Z,bcY_); ++k )
			for( Ordinal j=space()->begin(V,Y,bcY_); j<=space()->end(V,Y,bcY_); ++j )
				for( Ordinal i=space()->begin(V,X,bcY_); i<=space()->end(V,X,bcY_); ++i )
					y.getField(V).at(i,j,k) = innerStencV( x, i, j, k );

		if( 3==space_->dim() )  {

			x.exchange(Z);
			for( Ordinal k=space()->begin(W,Z,bcY_); k<=space()->end(W,Z,bcY_); ++k )
				for( Ordinal j=space()->begin(W,Y,bcY_); j<=space()->end(W,Y,bcY_); ++j )
					for( Ordinal i=space()->begin(W,X,bcY_); i<=space()->end(W,X,bcY_); ++i )
						y.getField(W).at(i,j,k) = innerStencW( x, i, j, k );
		}
  }


  void applyJ( RangeFieldT& y ) const {

		// BC scaling 
		const Scalar& eps = 0.1;
		//const Scalar& eps = 1.;
		
		for( int dir=0; dir<3; ++dir ) {
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

		for( int dir=0; dir<space_->dim(); ++dir )
			x.exchange( dir, dir );

		if( 3==space_->dim() )  {

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
		return( c_[dir] );
	}

	constexpr const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const {
		return( c_[dir][ off - space_->gl(dir) + i*( space_->gu(dir) - space_->gl(dir) + 1) ] );
	}

	constexpr const Scalar& getCTrans( const ECoord& dir, Ordinal i, Ordinal off ) const {
		return( cT_[dir][ off - space_->dl(dir) + i*( space_->du(dir) - space_->dl(dir) + 1) ] );
	}

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "\n--- " << getLabel() << " ---\n";
    //out << " --- stencil: ---";
    for( int dir=0; dir<3; ++dir ) {
			out << "\ndir: " << toString(static_cast<ECoord>(dir)) << "\n";

      for( int i=0; i<=space_->nLoc(dir); ++i ) {
        out << "\ni: " << i << "\t(";
        for( int ii=space_->gl(dir); ii<=space_->gu(dir); ++ii ) {
          out << getC( static_cast<ECoord>(dir), i, ii ) <<", ";
        }
        out << ")\n";
      }
      out << "\n";
    }

		out << "--- " << getLabel() << "^T ---\n";
		out << " --- stencil: ---";
		for( int dir=0; dir<3; ++dir ) {
			out << "\ndir: " << toString(static_cast<ECoord>(dir)) << "\n";

			for( int i=0; i<=space_->nLoc(dir); ++i ) {
        out << "\ni: " << i << "\t(";
				for( int k=space_->dl(dir); k<=space_->du(dir); ++k ) {
					out << getCTrans(static_cast<ECoord>(dir),i,k ) << ", ";
				}
        out << ")\n";
			}
			out << "\n";
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
			gradT += getC(X,i,ii)*x.getField(U).at(i+ii,j,k);

		for( int jj=space_->dl(Y); jj<=space_->du(Y); ++jj ) 
			gradT += getC(Y,j,jj)*x.getField(V).at(i,j+jj,k);

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
