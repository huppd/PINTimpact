#pragma once
#ifndef PIMPACT_DIVGRADOP_HPP
#define PIMPACT_DIVGRADOP_HPP


#include "Pimpact_DivOp.hpp"
#include "Pimpact_GradOp.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_Types.hpp"
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
      Belos::ETrans trans=Belos::NOTRANS ) const {

		Teuchos::RCP< VectorField<SpaceT> > temp = create<VectorField>( space() );

		switch( trans ) {
			case Belos::NOTRANS : {
				grad_->apply( x, *temp );
				div_->apply( *temp, y );
				break;
			}
			case Belos::TRANS : {
				div_->apply( x, *temp );
				grad_->apply( *temp, y );
				break;
			}
			case Belos::CONJTRANS : {
				div_->apply( x, *temp );
				grad_->apply( *temp, y );
				break;
			}
		}
  }


	void computeResidual( const RangeFieldT& b, const DomainFieldT& x, RangeFieldT& res ) const {
		apply( x, res );
		res.add( 1., b, -1., res );
	}


	/// \test
	/// \todo add eps
	void applyInvDiag( const DomainFieldT& x, RangeFieldT& y ) const {

		//const Scalar& eps = 0.1;

		for( OT k=space()->begin(S,Z,false); k<=space()->end(S,Z,false); ++k )
			for( OT j=space()->begin(S,Y,false); j<=space()->end(S,Y,false); ++j )
				for( OT i=space()->begin(S,X,false); i<=space()->end(S,X,false); ++i ) {
					ST diag = 0.;

					// X direction
					for( OT ii=space()->dl(X); ii<=space()->du(X); ++ii ) {
						if( 0<space()->getBCLocal()->getBCL(X) && i+ii==0 ) {
							for( OT iii=0; iii<=space()->du(X); ++iii )
								diag -= div_->getC( X, i, ii ) * grad_->getC(X,1+iii,-iii)
									* space()->getInterpolateV2S()->getC(X,1,iii) /
									space()->getInterpolateV2S()->getC(X,1,-1);
						}
						else if( 0<space()->getBCLocal()->getBCU(X) && i+ii==space()->end(S,X,false) ) {
							for( OT iii=space()->dl(X); iii<=-1; ++iii )
								diag -= div_->getC( X, i, ii ) * grad_->getC( X, i+iii, -iii)
									* space()->getInterpolateV2S()->getC(X,i,iii) /
									space()->getInterpolateV2S()->getC(X,i,0);
						}
						else
							diag += div_->getC( X, i, ii )*grad_->getC( X, i+ii, -ii );
					}

					// Y direction
					for( OT jj=space()->dl(Y); jj<=space()->du(Y); ++jj ) {
						if( 0<space()->getBCLocal()->getBCL(Y) && j+jj==0 ) {
							for( OT jjj=0; jjj<=space()->du(Y); ++jjj )
								diag -= div_->getC( Y, j, jj ) * grad_->getC( Y, 1+jjj, -jjj)
									* space()->getInterpolateV2S()->getC(Y,1,jjj) /
									space()->getInterpolateV2S()->getC(Y,1,-1);
						}
						else if( 0<space()->getBCLocal()->getBCU(Y) && j+jj==space()->end(S,Y,false) ) {
							for( OT jjj=space()->dl(Y); jjj<=-1; ++jjj )
								diag -= div_->getC( Y, j, jj ) * grad_->getC( Y, j+jjj, -jjj)
									* space()->getInterpolateV2S()->getC(Y,j,jjj) /
									space()->getInterpolateV2S()->getC(Y,j,0);
						}
						else
							diag += div_->getC( Y, j, jj )*grad_->getC( Y, j+jj, -jj );
					}

					if( 3==space()->dim() ) {
						// Z direction
						for( OT kk=space()->dl(Z); kk<=space()->du(Z); ++kk ) {
							if( 0<space()->getBCLocal()->getBCL(Z) && k+kk==0 ) {
								for( OT kkk=0; kkk<=space()->du(Z); ++kkk )
									diag -= div_->getC( Z, k, kk ) * grad_->getC( Z, 1+kkk, -kkk)
										* space()->getInterpolateV2S()->getC(Z,1,kkk) /
										space()->getInterpolateV2S()->getC(Z,1,-1);
							}
							else if( 0<space()->getBCLocal()->getBCU(Z) && k+kk==space()->end(S,Z,false) ) {
								for( OT kkk=space()->dl(Z); kkk<=-1; ++kkk )
									diag -= div_->getC( Z, k, kk ) * grad_->getC( Z, k+kkk, -kkk)
										* space()->getInterpolateV2S()->getC(Z,k,kkk) /
										space()->getInterpolateV2S()->getC(Z,k,0);
							}
							else
								diag += div_->getC( Z, k, kk )*grad_->getC( Z, k+kk, -kk );
						}
					}

#ifndef NDEBUG
					TEUCHOS_TEST_FOR_EXCEPT( 0==diag );
#endif
					y.at(i,j,k) = x.at(i,j,k)/std::abs( diag );
				}

		y.changed();
	}

	//inline ST innerStenc3D( const VectorField<SpaceT> & x,
			//const OT& i, const OT& j, const OT& k ) const {

		//ST div = 0.;

		//if( i==space()->begin(S,X) && space()->getBCLocal()->getBCL(X)==DirichletBC )
			//for( int ii=space()->dl(X); ii<=space()->du(X); ++ii ) 
				//div -= space()->getInterpolateV2S()->getC(X,i,ii)*x.getField(U).at(i+ii,j,k);
		//else if( i==space()->end(S,X) && space()->getBCLocal()->getBCU(X)==DirichletBC )
			//for( int ii=space()->dl(X); ii<=space()->du(X); ++ii ) 
				//div += space()->getInterpolateV2S()->getC(X,i,ii)*x.getField(U).at(i+ii,j,k);
		//else
			//for( int ii=space()->dl(X); ii<=space()->du(X); ++ii ) 
				//div += div_->getC(X,i,ii)*x.getField(U).at(i+ii,j,k);

		//if( j==space()->begin(S,Y) && space()->getBCLocal()->getBCL(Y)==DirichletBC )
			//for( int jj=space()->dl(Y); jj<=space()->du(Y); ++jj ) 
				//div -= space()->getInterpolateV2S()->getC(Y,j,jj)*x.getField(V).at(i,j+jj,k);
		//else if( j==space()->end(S,Y) && space()->getBCLocal()->getBCU(Y)==DirichletBC )
			//for( int jj=space()->dl(Y); jj<=space()->du(Y); ++jj ) 
				//div += space()->getInterpolateV2S()->getC(Y,j,jj)*x.getField(V).at(i,j+jj,k);
		//else
			//for( int jj=space()->dl(Y); jj<=space()->du(Y); ++jj ) 
				//div += div_->getC(Y,j,jj)*x.getField(V).at(i,j+jj,k);

		//if( k==space()->begin(S,Z) && space()->getBCLocal()->getBCL(Z)==DirichletBC )
			//for( int kk=space()->dl(Z); kk<=space()->du(Z); ++kk ) 
				//div -= space()->getInterpolateV2S()->getC(Z,k,kk)*x.getField(W).at(i,j,k+kk);
		//else if( k==space()->end(S,Z) && space()->getBCLocal()->getBCU(Z)==DirichletBC )
			//for( int kk=space()->dl(Z); kk<=space()->du(Z); ++kk ) 
				//div += space()->getInterpolateV2S()->getC(Z,k,kk)*x.getField(W).at(i,j,k+kk);
		//else
			//for( int kk=space()->dl(Z); kk<=space()->du(Z); ++kk ) 
				//div += div_->getC(Z,k,kk)*x.getField(W).at(i,j,k+kk);

		//return( div );
	//}

	//inline ST innerStenc2D( const VectorField<SpaceT> & x,
			//const OT& i, const OT& j, const OT& k ) const {

		//ST div = 0.;

		//if( i==space()->begin(S,X) && space()->getBCLocal()->getBCL(X)==DirichletBC )
			//for( int ii=space()->dl(X); ii<=space()->du(X); ++ii ) 
				//div -= space()->getInterpolateV2S()->getC(X,i,ii)*x.getField(U).at(i+ii,j,k);
		//else if( i==space()->end(S,X) && space()->getBCLocal()->getBCU(X)==DirichletBC )
			//for( int ii=space()->dl(X); ii<=space()->du(X); ++ii ) 
				//div += space()->getInterpolateV2S()->getC(X,i,ii)*x.getField(U).at(i+ii,j,k);
		//else
			//for( int ii=space()->dl(X); ii<=space()->du(X); ++ii ) 
				//div += div_->getC(X,i,ii)*x.getField(U).at(i+ii,j,k);

		//if( j==space()->begin(S,Y) && space()->getBCLocal()->getBCL(Y)==DirichletBC )
			//for( int jj=space()->dl(Y); jj<=space()->du(Y); ++jj ) 
				//div -= space()->getInterpolateV2S()->getC(Y,j,jj)*x.getField(V).at(i,j+jj,k);
		//else if( j==space()->end(S,Y) && space()->getBCLocal()->getBCU(Y)==DirichletBC )
			//for( int jj=space()->dl(Y); jj<=space()->du(Y); ++jj ) 
				//div += space()->getInterpolateV2S()->getC(Y,j,jj)*x.getField(V).at(i,j+jj,k);
		//else
			//for( int jj=space()->dl(Y); jj<=space()->du(Y); ++jj ) 
				//div += div_->getC(Y,j,jj)*x.getField(V).at(i,j+jj,k);

		//return( div );
	//}

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( true ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return( div_->space() ); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

	constexpr const std::string getLabel() const { return( "DivGrad" ); };

  void print( std::ostream& out=std::cout ) const {
		out << "---" << getLabel() << "---\n";
		div_->print( out );
		grad_->print( out );
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
