#pragma once
#ifndef PIMPACT_DIVGRADOP_HPP
#define PIMPACT_DIVGRADOP_HPP


#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{



/// \brief "laplace" for pressure.
/// \ingroup BaseOperator
/// \note todo not workin properly?
/// \warning does not hold test.
template<class SpT>
class DivGradOp {

public:

  using SpaceT = SpT;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

	using ST = typename SpaceT::Scalar;
	using OT = typename SpaceT::Ordinal;

  Teuchos::RCP<DivOp<SpaceT> > div_;
  Teuchos::RCP<GradOp<SpaceT> > grad_;

public:

  DivGradOp( const Teuchos::RCP<const SpaceT>& space ):
    div_ ( create<DivOp>( space ) ),
    grad_( create<GradOp>( space ) ) {};

  DivGradOp(
      const Teuchos::RCP< DivOp<SpaceT> >& div,
      const Teuchos::RCP< GradOp<SpaceT> >& grad ):
    div_ (div),
    grad_(grad) {};

  void apply(const DomainFieldT& x, RangeFieldT& y,
      const Belos::ETrans& trans=Belos::NOTRANS, const Add& add=Add::N ) const {

		VectorField<SpaceT> temp( space() );

		switch( trans ) {
			case Belos::NOTRANS : {
				grad_->applyG( x, temp, add );
				grad_->applyJ( temp );
				div_->apply( temp, y, add );
				break;
			}
			case Belos::TRANS : {
				div_->apply( x, temp, add );
				grad_->apply( temp, y, add );
				break;
			}
			case Belos::CONJTRANS : {
				div_->apply( x, temp, add );
				grad_->apply( temp, y, add );
				break;
			}
		}
	}


	void computeResidual( const RangeFieldT& b, const DomainFieldT& x, RangeFieldT& res ) const {
		apply( x, res );
		res.add( 1., b, -1., res );
	}


	/// \todo fix for periodic BC
	void applyInvDiag( const DomainFieldT& x, RangeFieldT& y ) const {


		for( OT k=space()->si(F::S,Z); k<=space()->ei(F::S,Z); ++k )
			for( OT j=space()->si(F::S,Y); j<=space()->ei(F::S,Y); ++j )
				for( OT i=space()->si(F::S,X); i<=space()->ei(F::S,X); ++i ) {

					ST diag = ( (3==SpaceT::sdim)?innerDiag3D(i,j,k):innerDiag2D(i,j,k) );

					assert( 0!=diag );
					y(i,j,k) = x(i,j,k)/std::fabs( diag );
				}

		y.changed();
	}

  Teuchos::RCP<const DivOp<SpaceT> > getDivOp() const { return( div_ ); }
  Teuchos::RCP<const GradOp<SpaceT> > getGradOp() const { return( grad_ ); }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( true ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return( div_->space() ); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

	constexpr const std::string getLabel() const { return( "DivGrad" ); };

  void print( std::ostream& out=std::cout ) const {
		out << "---" << getLabel() << "---\n";
		div_->print( out );
		grad_->print( out );
		space()->getInterpolateV2S()->print();
  }


	ST innerDiag3D( const OT& i, const OT& j,
			const OT& k ) const {

		ST diag = 0.;
		//const ST& eps = 0.1;
		const ST& eps = 1./static_cast<ST>(GradOp<SpaceT>::epsI);

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->si(F::S,X) ) ||
			(               space()->getBCLocal()->getBCU(X) > 0 && i==space()->ei(F::S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->si(F::S,Y) ) ||
			(               space()->getBCLocal()->getBCU(Y) > 0 && j==space()->ei(F::S,Y) ) ;
		const bool bcZ = (space()->getBCLocal()->getBCL(Z) > 0 && k==space()->si(F::S,Z) ) ||
			(               space()->getBCLocal()->getBCU(Z) > 0 && k==space()->ei(F::S,Z) ) ;

		const ST epsX = ( (bcY||bcZ)?eps:1. );
		const ST epsY = ( (bcX||bcZ)?eps:1. );
		const ST epsZ = ( (bcX||bcY)?eps:1. );

		// X direction
		for( OT ii=space()->dl(X); ii<=space()->du(X); ++ii ) {
			if( 0<space()->getBCLocal()->getBCL(X) && i+ii==space()->si(F::U,X,B::Y) ) {
				for( OT iii=0; iii<=space()->du(X); ++iii ) 
					diag -= div_->getC( X, i, ii ) * epsX * grad_->getC( X, 1+iii, -iii-ii-1 )
						* space()->getInterpolateV2S()->getC( X, 1, iii ) /
						space()->getInterpolateV2S()->getC( X, 1, -1 );
			}
			else if( 0<space()->getBCLocal()->getBCU(X) && i+ii==space()->ei(F::U,X,B::Y) ) {
				for( OT iii=space()->dl(X); iii<=-1; ++iii )
					diag -= div_->getC( X, i, ii ) * epsX * grad_->getC( X, space()->ei(F::U,X,B::Y)+iii, -iii-ii )
						* space()->getInterpolateV2S()->getC( X, space()->ei(F::U,X,B::Y), iii ) /
						space()->getInterpolateV2S()->getC( X, space()->ei(F::U,X,B::Y), 0 );
			}
			else if( i+ii>=0 && i+ii<=space()->nLoc(X) )
				diag += div_->getC( X, i, ii ) * epsX * grad_->getC( X, i+ii, -ii );
		}

		// Y direction
		for( OT jj=space()->dl(Y); jj<=space()->du(Y); ++jj ) {
			if( 0<space()->getBCLocal()->getBCL(Y) && j+jj==space()->si(F::V,Y,B::Y) ) {
				for( OT jjj=0; jjj<=space()->du(Y); ++jjj )
					diag -= div_->getC( Y, j, jj ) * epsY * grad_->getC( Y, 1+jjj, -jjj+j-1 )
						* space()->getInterpolateV2S()->getC( Y, 1, jjj ) /
						space()->getInterpolateV2S()->getC( Y, 1, -1 );
			}
			else if( 0<space()->getBCLocal()->getBCU(Y) && j+jj==space()->ei(F::V,Y,B::Y) ) {
				for( OT jjj=space()->dl(Y); jjj<=-1; ++jjj )
					diag -= div_->getC( Y, j, jj ) * epsY * grad_->getC( Y, space()->ei(F::V,Y,B::Y)+jjj, -jjj-jj )
						* space()->getInterpolateV2S()->getC( Y, space()->ei(F::V,Y,B::Y), jjj ) /
						space()->getInterpolateV2S()->getC( Y, space()->ei(F::V,Y,B::Y), 0 );
			}
			else if( j+jj>=0 && j+jj<=space()->nLoc(Y) )
				diag += div_->getC( Y, j, jj )*epsY*grad_->getC( Y, j+jj, -jj );
		}

		// Z direction
		for( OT kk=space()->dl(Z); kk<=space()->du(Z); ++kk ) {
			if( 0<space()->getBCLocal()->getBCL(Z) && k+kk==space()->si(F::W,Z,B::Y) ) {
				for( OT kkk=0; kkk<=space()->du(Z); ++kkk )
					diag -= div_->getC( Z, k, kk ) * epsZ * grad_->getC( Z, 1+kkk, -kkk+k-1 )
						* space()->getInterpolateV2S()->getC( Z, 1, kkk ) /
						space()->getInterpolateV2S()->getC( Z, 1, -1 );
			}
			else if( 0<space()->getBCLocal()->getBCU(Z) && k+kk==space()->ei(F::W,Z,B::Y) ) {
				for( OT kkk=space()->dl(Z); kkk<=-1; ++kkk )
					diag -= div_->getC( Z, k, kk ) * epsZ * grad_->getC( Z, space()->ei(F::W,Z,B::Y)+kkk, -kkk-kk )
						* space()->getInterpolateV2S()->getC(Z,space()->ei(F::W,Z,B::Y),kkk) /
						space()->getInterpolateV2S()->getC(Z,space()->ei(F::W,Z,B::Y),0);
			}
			else if( k+kk>=0 && k+kk<=space()->nLoc(Z) )
				diag += div_->getC( Z, k, kk ) * epsZ * grad_->getC( Z, k+kk, -kk );
		}

		return( diag );
	}

	ST innerDiag2D( const OT& i, const OT& j, const OT& k ) const {

		ST diag = 0.;
		//const ST& eps = 0.1;
		const ST& eps = 1./static_cast<ST>(GradOp<SpaceT>::epsI);

		const bool bcX = (space()->getBCLocal()->getBCL(X) > 0 && i==space()->si(F::S,X) ) ||
			(               space()->getBCLocal()->getBCU(X) > 0 && i==space()->ei(F::S,X) ) ;
		const bool bcY = (space()->getBCLocal()->getBCL(Y) > 0 && j==space()->si(F::S,Y) ) ||
			(               space()->getBCLocal()->getBCU(Y) > 0 && j==space()->ei(F::S,Y) ) ;

		const ST epsX = ( bcY?eps:1. );
		const ST epsY = ( bcX?eps:1. );

		// X direction
		for( OT ii=space()->dl(X); ii<=space()->du(X); ++ii ) {
			if( 0<space()->getBCLocal()->getBCL(X) && i+ii==space()->si(F::U,X,B::Y) ) {
				for( OT iii=0; iii<=space()->du(X); ++iii ) 
					diag -= div_->getC( X, i, ii ) * epsX * grad_->getC( X, 1+iii, -iii-ii-1 )
						* space()->getInterpolateV2S()->getC( X, 1, iii ) /
						space()->getInterpolateV2S()->getC( X, 1, -1 );
			}
			else if( 0<space()->getBCLocal()->getBCU(X) && i+ii==space()->ei(F::U,X,B::Y) ) {
				for( OT iii=space()->dl(X); iii<=-1; ++iii )
					diag -= div_->getC( X, i, ii ) * epsX * grad_->getC( X, space()->ei(F::U,X,B::Y)+iii, -iii-ii )
						* space()->getInterpolateV2S()->getC( X, space()->ei(F::U,X,B::Y), iii ) /
						space()->getInterpolateV2S()->getC( X, space()->ei(F::U,X,B::Y), 0 );
			}
			else if( i+ii>=0 && i+ii<=space()->nLoc(X) )
				diag += div_->getC( X, i, ii ) * epsX * grad_->getC( X, i+ii, -ii );
		}

		// Y direction
		for( OT jj=space()->dl(Y); jj<=space()->du(Y); ++jj ) {
			if( 0<space()->getBCLocal()->getBCL(Y) && j+jj==space()->si(F::V,Y,B::Y) ) {
				for( OT jjj=0; jjj<=space()->du(Y); ++jjj )
					diag -= div_->getC( Y, j, jj ) * epsY * grad_->getC( Y, 1+jjj, -jjj+j-1 )
						* space()->getInterpolateV2S()->getC( Y, 1, jjj ) /
						space()->getInterpolateV2S()->getC( Y, 1, -1 );
			}
			else if( 0<space()->getBCLocal()->getBCU(Y) && j+jj==space()->ei(F::V,Y,B::Y) ) {
				for( OT jjj=space()->dl(Y); jjj<=-1; ++jjj )
					diag -= div_->getC( Y, j, jj ) * epsY * grad_->getC( Y, space()->ei(F::V,Y,B::Y)+jjj, -jjj-jj )
						* space()->getInterpolateV2S()->getC( Y, space()->ei(F::V,Y,B::Y), jjj ) /
						space()->getInterpolateV2S()->getC( Y, space()->ei(F::V,Y,B::Y), 0 );
			}
			else if( j+jj>=0 && j+jj<=space()->nLoc(Y) )
				diag += div_->getC( Y, j, jj )*epsY*grad_->getC( Y, j+jj, -jj );
		}

		return( diag );
	}

}; // end of DivGradOp



/// \relates DivGradOp
template<class SpaceT>
Teuchos::RCP< DivGradOp<SpaceT> > createDivGradOp(
		const Teuchos::RCP< DivOp<SpaceT> >& div,
		const Teuchos::RCP< GradOp<SpaceT> >& grad ) {

	return(
			Teuchos::rcp( new DivGradOp<SpaceT>( div, grad ) )
			);
}


} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::DivGradOp< Pimpact::Space<double,int,3,2> >;
extern template class Pimpact::DivGradOp< Pimpact::Space<double,int,3,4> >;
extern template class Pimpact::DivGradOp< Pimpact::Space<double,int,4,2> >;
extern template class Pimpact::DivGradOp< Pimpact::Space<double,int,4,4> >;
#endif


#endif // end of #ifndef PIMPACT_DIVGRADOP_HPP
