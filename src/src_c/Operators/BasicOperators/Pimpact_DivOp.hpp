#pragma once
#ifndef PIMPACT_DIVOP_HPP
#define PIMPACT_DIVOP_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_extern_FDCoeff.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Stencil.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{



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

	static const int dimNC = ST::dimNC;
	static const int dim = ST::dimension;

	using SW = StencilWidths<dim,dimNC>;

	using StencD = Stencil< Scalar, Ordinal, 0, SW::DL(0), SW::DU(0) >;
	using StencG = Stencil< Scalar, Ordinal, 0, SW::GL(0), SW::GU(0) >;

	using TD = const Teuchos::Tuple< StencD, ST::sdim >;
	using TG = const Teuchos::Tuple< StencG, ST::sdim >;

  Teuchos::RCP<const SpaceT> space_;

  TD c_;
  TG cT_;

public:

	DivOp( const Teuchos::RCP<const SpaceT>& space ):
		space_(space) {

			for( int dir=0; dir<ST::sdim; ++dir ) {
				// Divergence stencil
				c_[dir] = StencD( space_->nLoc(dir) );

				FD_getDiffCoeff(
						space_->nLoc(dir),
						space_->bl(dir),
						space_->bu(dir),
						space_->dl(dir),
						space_->du(dir),
						space_->getBCLocal()->getBCL(dir),
						space_->getBCLocal()->getBCU(dir),
						space_->getShift(dir),
						3,
						dir+1,
						1,
						0,
						//true, // mapping
						false,
						space_->getStencilWidths()->getDimNcbD(dir),
						space_->getStencilWidths()->getNcbD(dir),
						space_->getCoordinatesLocal()->getX( dir, dir ),
						space_->getCoordinatesLocal()->getX( dir, EField::S ),
						c_[dir].get() );

				// Divergence stencil transposed
				cT_[dir] = StencG( space_->nLoc(dir) );

				Ordinal nTempG = ( space_->nGlo(dir) + space_->bu(dir) - space_->bl(dir) + 1 )
					*( space_->bu(dir) - space_->bl(dir) + 1);

				Stencil< Scalar, Ordinal, SW::BL(0), SW::BL(0), SW::BU(0) >
					cG1( space_->nGlo(dir) + space_->bu(dir) );
				Stencil< Scalar, Ordinal, SW::BL(0), SW::BL(0), SW::BU(0) >
					cG2( space_->nGlo(dir) + space_->bu(dir) );

				for( Ordinal i = space_->begin(S,dir); i<=space_->end(S,dir); ++i )
					for( Ordinal ii = space_->dl(dir); ii<=space_->du(dir); ++ii )
						cG1( i+space_->getShift(dir), ii )= getC( static_cast<ECoord>(dir), i, ii );

				MPI_Allreduce(
						cG1.get(),    		                          // const void *sendbuf,
						cG2.get(),    		                          // void *recvbuf,
						nTempG,			                                // int count,
						MPI_REAL8,	                                // MPI_Datatype datatype,
						MPI_SUM,		                                // MPI_Op op,
						space_->getProcGrid()->getCommSlice(dir) ); // MPI_Comm comm )


				if( -1==space_->getBCGlobal()->getBCL(dir) ) {

					Ordinal ls1 = space_->getStencilWidths()->getLS(dir);
					Ordinal M1 = space_->nGlo(dir);

					for( Ordinal i=space->bl(dir); i<=-1; ++i )
						for( Ordinal ii=space->bl(dir); ii<=space->bu(dir); ++ii )
							cG2( 2+ls1+i, ii ) = cG2( M1+1+ls1+i, ii );

					for( Ordinal i=1; i<=space->bu(dir); ++i )
						for( Ordinal ii=space->bl(dir); ii<=space->bu(dir); ++ii )
							cG2( M1+ls1+i, ii ) = cG2( 1+ls1+i, ii );
				}

				for( Ordinal i = space_->begin(dir,dir,With::B);
						i<=space_->end(dir,dir,With::B); ++i ) 
					for( Ordinal ii=space->gl(dir); ii<=space->gu(dir); ++ii )
						cT_[dir]( i, ii ) = cG2( i+ii+space_->getShift(dir), -ii );
			}
  };



	void apply( const DomainFieldT& x, RangeFieldT& y, const Add& add=Add::No ) const {

		for( int dir=0; dir<ST::sdim; ++dir )
			x.exchange( dir, dir );

		if( 3==ST::sdim ) {

			for( Ordinal k=space()->begin(S,Z); k<=space()->end(S,Z); ++k )
				for( Ordinal j=space()->begin(S,Y); j<=space()->end(S,Y); ++j )
					for( Ordinal i=space()->begin(S,X); i<=space()->end(S,X); ++i ) {
						if( Add::No==add ) y(i,j,k) = 0.;
						y(i,j,k) += innerStenc3D( x, i, j, k );
					}
		}
		else {

			for( Ordinal k=space()->begin(S,Z); k<=space()->end(S,Z); ++k )
				for( Ordinal j=space()->begin(S,Y); j<=space()->end(S,Y); ++j )
					for( Ordinal i=space()->begin(S,X); i<=space()->end(S,X); ++i ) {
						if( Add::No==add ) y(i,j,k) = 0.;
						y(i,j,k) += innerStenc2D( x, i, j, k );
					}
		}

		y.changed();
	}


	void apply( const RangeFieldT& x, DomainFieldT& y, const Add& add=Add::No ) const {

		x.exchange(X);
		for( Ordinal k=space()->begin(U,Z,With::B); k<=space()->end(U,Z,With::B); ++k )
			for( Ordinal j=space()->begin(U,Y,With::B); j<=space()->end(U,Y,With::B); ++j )
				for( Ordinal i=space()->begin(U,X,With::B); i<=space()->end(U,X,With::B); ++i ) {
					if( Add::No==add ) y.getField(U)(i,j,k) = 0.;
					y.getField(U)(i,j,k) += innerStencU( x, i, j, k );
				}

		x.exchange(Y);
		for( Ordinal k=space()->begin(V,Z,With::B); k<=space()->end(V,Z,With::B); ++k )
			for( Ordinal j=space()->begin(V,Y,With::B); j<=space()->end(V,Y,With::B); ++j )
				for( Ordinal i=space()->begin(V,X,With::B); i<=space()->end(V,X,With::B); ++i ) {
					if( Add::No==add ) y.getField(V)(i,j,k) = 0.;
					y.getField(V)(i,j,k) += innerStencV( x, i, j, k );
				}

		if( 3==SpaceT::sdim )  {

			x.exchange(Z);
			for( Ordinal k=space()->begin(W,Z,With::B); k<=space()->end(W,Z,With::B); ++k )
				for( Ordinal j=space()->begin(W,Y,With::B); j<=space()->end(W,Y,With::B); ++j )
					for( Ordinal i=space()->begin(W,X,With::B); i<=space()->end(W,X,With::B); ++i ) {
						if( Add::No==add ) y.getField(W)(i,j,k) = 0.;
						y.getField(W)(i,j,k) += innerStencW( x, i, j, k );
					}
		}

		y.extrapolateBC( Belos::TRANS );

		// BC scaling 
		const Scalar& eps = 1.e-1;
		for( int dir=0; dir<SpaceT::sdim; ++dir ) {
			With bc2 = With::B;
			if( 0!=dir ) {
				if( space()->getBCLocal()->getBCL(X) > 0 ) {
					Ordinal i = space()->begin(dir,X,With::B);
					for( Ordinal k=space()->begin(dir,Z, bc2); k<=space()->end(dir,Z,bc2); ++k )
						for( Ordinal j=space()->begin(dir,Y,bc2); j<=space()->end(dir,Y,bc2); ++j )
							y.getField(dir)(i,j,k) *= eps;  
				}
				if( space()->getBCLocal()->getBCU(X) > 0 ) {
					Ordinal i = space()->end(dir,X,With::B);
					for( Ordinal k=space()->begin(dir,Z,bc2); k<=space()->end(dir,Z,bc2); ++k )
						for( Ordinal j=space()->begin(dir,Y,bc2); j<=space()->end(dir,Y,bc2); ++j )
							y.getField(dir)(i,j,k) *= eps;  
				}
				bc2 = With::noB;
			}

			if( 1!=dir ) {
				if( space()->getBCLocal()->getBCL(Y) > 0 ) {
					Ordinal j = space()->begin(dir,Y,With::B);
					for( Ordinal k=space()->begin(dir,Z,bc2); k<=space()->end(dir,Z,bc2); ++k )
						for( Ordinal i=space()->begin(dir,X,bc2); i<=space()->end(dir,X,bc2); ++i ) 
							y.getField(dir)(i,j,k) *= eps;  
				}
				if( space()->getBCLocal()->getBCU(Y) > 0 ) {
					Ordinal j = space()->end(dir,Y,With::B);
					for( Ordinal k=space()->begin(dir,Z,bc2); k<=space()->end(dir,Z,bc2); ++k )
						for( Ordinal i=space()->begin(dir,X,bc2); i<=space()->end(dir,X,bc2); ++i )
							y.getField(dir)(i,j,k) *= eps;  
				}
				bc2 = With::noB;
			}

			if( 2!=dir ) {
				if( space()->getBCLocal()->getBCL(Z) > 0 ) {
					Ordinal k = space()->begin(dir,Z,With::B);
					for( Ordinal j=space()->begin(dir,Y,bc2); j<=space()->end(dir,Y,bc2); ++j )
						for( Ordinal i=space()->begin(dir,X,bc2); i<=space()->end(dir,X,bc2); ++i )
							y.getField(dir)(i,j,k) *= eps;  
				}
				if( space()->getBCLocal()->getBCU(Z) > 0 ) {
					Ordinal k = space()->end(dir,Z,With::B);
					for( Ordinal j=space()->begin(dir,Y,bc2); j<=space()->end(dir,Y,bc2); ++j )
						for( Ordinal i=space()->begin(dir,X,bc2); i<=space()->end(dir,X,bc2); ++i )
							y.getField(dir)(i,j,k) *= eps;  
				}
				bc2 = With::noB;
			}
		}

		y.changed();
	}


  void assignField( const RangeFieldT& mv ) const {};
  void assignField( const DomainFieldT& mv ) const {};

  bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(space_); };

	constexpr const Scalar* getC( const ECoord& dir ) const {
		return( c_[dir].get() );
	}

	constexpr const Scalar& getC( const ECoord& dir, Ordinal i, Ordinal off ) const {
		return( c_[dir]( i, off ) );
	}

	constexpr const Scalar& getCTrans( const ECoord& dir, Ordinal i, Ordinal off ) const {
		return( cT_[dir]( i, off ) );
	}

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

	void print( std::ostream& out=std::cout ) const {
		out << "\n--- " << getLabel() << " ---\n";
		out << " --- stencil: ---";
		for( int dir=0; dir<ST::sdim; ++dir ) {
			out << "\ndir: " << toString( static_cast<ECoord>(dir) ) << "\n";
			c_[dir].print( out );
		}

		out << "--- " << getLabel() << "^T ---\n";
		out << " --- stencil: ---";
		for( int dir=0; dir<ST::sdim; ++dir ) {
			out << "\ndir: " << toString(static_cast<ECoord>(dir)) << "\n\n";
			cT_[dir].print( out );
		}
	}

	const std::string getLabel() const { return( "Div" ); };

protected:

	inline constexpr Scalar innerStenc3D( const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar div = 0.;

		for( int ii=space_->dl(X); ii<=space_->du(X); ++ii ) 
			div += getC(X,i,ii)*x.getField(U)(i+ii,j,k);

		for( int jj=space_->dl(Y); jj<=space_->du(Y); ++jj ) 
			div += getC(Y,j,jj)*x.getField(V)(i,j+jj,k);

		for( int kk=space_->dl(Z); kk<=space_->du(Z); ++kk ) 
			div += getC(Z,k,kk)*x.getField(W)(i,j,k+kk);

		return( div );
	}

	inline constexpr Scalar innerStenc2D( const DomainFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar div = 0.;

		for( int ii=space_->dl(X); ii<=space_->du(X); ++ii ) 
			div += getC(X,i,ii)*x.getField(U)(i+ii,j,k);

		for( int jj=space_->dl(Y); jj<=space_->du(Y); ++jj ) 
			div += getC(Y,j,jj)*x.getField(V)(i,j+jj,k);

		return( div );
	}

	inline constexpr Scalar innerStencU( const RangeFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar divT = 0.;

		for( int ii=space_->gl(X); ii<=space_->gu(X); ++ii ) 
			divT += getCTrans(X,i,ii)*x(i+ii,j,k);

		return( divT );
	}

	inline constexpr Scalar innerStencV( const RangeFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar divT = 0.;

		for( int jj=space_->gl(Y); jj<=space_->gu(Y); ++jj ) 
			divT += getCTrans(Y,j,jj)*x(i,j+jj,k);

		return( divT );
	}

	inline constexpr Scalar innerStencW( const RangeFieldT& x,
			const Ordinal& i, const Ordinal& j, const Ordinal& k ) const {

		Scalar divT = 0.;

		for( int kk=space_->gl(Z); kk<=space_->gu(Z); ++kk ) 
			divT += getCTrans(Z,k,kk)*x(i,j,k+kk);

		return( divT );
	}

}; // end of class DivOp


} // end of namespace Pimpact



#ifdef COMPILE_ETI
extern template class Pimpact::DivOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::DivOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::DivOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::DivOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_DIVOP_HPP
