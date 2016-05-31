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


	void applyInvDiag( const DomainFieldT& x, RangeFieldT& y ) const {

		//for( int dir=0; dir<space()->dim(); ++dir )
		//x.exchange();

		if( 3==space()->dim() ) {

			for( OT k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				for( OT j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j )
					for( OT i=space()->sInd(S,X); i<=space()->eInd(S,X); ++i ) {
						ST diag = 
							div_->getC( X, i, 0) * grad_->getC( X, i, 0) +
							div_->getC( Y, j, 0) * grad_->getC( Y, j, 0) +
							div_->getC( Z, k, 0) * grad_->getC( Z, k, 0);

						// BC in X
						if( space()->sInd(S,X)==i &&  0<space()->getBCLocal()->getBCL(X) ) {
							diag -= div_->getC(X,i,-1) * grad_->getC( X, i, 0) *
								space()->getInterpolateV2S()->getC(X,i,0) /
								space()->getInterpolateV2S()->getC(X,i,-1);
						}

						if( space()->eInd(S,X)==i && 0<space()->getBCLocal()->getBCU(X) ) {
							diag -= div_->getC(X,i,0) * grad_->getC( X, i-1, 0) * 
								space()->getInterpolateV2S()->getC(X,i,-1) /
								space()->getInterpolateV2S()->getC(X,i, 0);
						}

						// BC in Y
						if( space()->sInd(S,Y)==j &&  0<space()->getBCLocal()->getBCL(Y) ) {
							diag -= div_->getC(Y,j,-1) * grad_->getC( Y, j, 0) *
								space()->getInterpolateV2S()->getC(Y,j,0) /
								space()->getInterpolateV2S()->getC(Y,j,-1);
						}

						if( space()->eInd(S,Y)==j && 0<space()->getBCLocal()->getBCU(Y) ) {
							diag -= div_->getC(Y,j,0) * grad_->getC( Y, j-1, 0) * 
								space()->getInterpolateV2S()->getC(Y,j,-1) /
								space()->getInterpolateV2S()->getC(Y,j, 0);
						}

						// BC in Z
						if( space()->sInd(S,Z)==k &&  0<space()->getBCLocal()->getBCL(Z) ) {
							diag -= div_->getC(Z,k,-1) * grad_->getC( Z, k, 0) *
								space()->getInterpolateV2S()->getC(Z,k,0) /
								space()->getInterpolateV2S()->getC(Z,k,-1);
						}

						if( space()->eInd(S,Z)==k && 0<space()->getBCLocal()->getBCU(Z) ) {
							diag -= div_->getC(Z,k,0) * grad_->getC( Z, k-1, 0) * 
								space()->getInterpolateV2S()->getC(Z,k,-1) /
								space()->getInterpolateV2S()->getC(Z,k, 0);
						}

						if( diag!=0 )
							y.at(i,j,k) = x.at(i,j,k)/std::abs( diag );
						else{
							y.at(i,j,k) = 0.;
							std::cout << "i: " << i << " j: " << j << " k: " << k << "\n";
						}
						//std::cout << y.at(i,j,k) << "\n";
					}
		}
		else{

			for( OT k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				for( OT j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j )
					for( OT i=space()->sInd(S,X); i<=space()->eInd(S,X); ++i ) {
						ST diag = 
							div_->getC( X, i, 0) * grad_->getC( X, i, 0) +
							div_->getC( Y, j, 0) * grad_->getC( Y, j, 0);

						if( space()->sInd(S,X)==i &&  0<space()->getBCLocal()->getBCL(X) )
							diag -= div_->getC(X,i,-1) * grad_->getC( X, i, 0) *
								space()->getInterpolateV2S()->getC(X,i,0) /
								space()->getInterpolateV2S()->getC(X,i,-1);
						
						if( space()->eInd(S,X)==i && 0<space()->getBCLocal()->getBCU(X) )
							diag -= div_->getC(X,i,1) * grad_->getC( X, i, 0) * x.at(i,j,k) *
								space()->getInterpolateV2S()->getC(X,i,-1) /
								space()->getInterpolateV2S()->getC(X,i, 0);

						y.at(i,j,k) = x.at(i,j,k)/diag;
					}
		}

		// Dirichlet boundary conditions in X
			//OT i = space()->sInd(S,X);
			//for( OT k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				//for( OT j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j ) {
					//y.at(i,j,k) -=
						//div_->getC(X,i,-1) * grad_->getC( X, i, 0) * x.at(i,j,k) *
						//space()->getInterpolateV2S()->getC(X,i,0) /
						//space()->getInterpolateV2S()->getC(X,i,-1);
				//}
			//std::cout  << space()->getInterpolateV2S()->getC(X,1,-1) << "\t\t";
			//std::cout  << space()->getInterpolateV2S()->getC(X,1, 0) << "\n";

		//}
		//if( 0<space()->getBCLocal()->getBCU(X) ) {
			//OT i = space()->eInd(S,X);
			//for( OT k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				//for( OT j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j )
					//y.at(i,j,k) -=
						//div_->getC(X,i,1) * grad_->getC( X, i, 0) * x.at(i,j,k) *
						//space()->getInterpolateV2S()->getC(X,i,-1) /
						//space()->getInterpolateV2S()->getC(X,i, 0);
		//}

		//// Dirichlet boundary conditions in Y
		//if( 0<space()->getBCLocal()->getBCL(Y) ) {
			//OT j=space()->sInd(S,Y);
			//for( OT k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				//for( OT i=space()->sInd(S,X); i<=space()->eInd(S,X); ++i )
					//y.at(i,j,k) -=
						//div_->getC(Y,j,-1) * grad_->getC( Y, j, 0) * x.at(i,j,k) *
						//space()->getInterpolateV2S()->getC(Y,j,0) /
						//space()->getInterpolateV2S()->getC(Y,j,-1);
		//}
		//if( 0<space()->getBCLocal()->getBCU(Y) ) {
			//OT j=space()->eInd(S,Y);
			//for( OT k=space()->sInd(S,Z); k<=space()->eInd(S,Z); ++k )
				//for( OT i=space()->sInd(S,X); i<=space()->eInd(S,X); ++i )
					//y.at(i,j,k) -=
						//div_->getC(Y,j,1) * grad_->getC( Y, j, 0) * x.at(i,j,k) *
						//space()->getInterpolateV2S()->getC(Y,j,-1) /
						//space()->getInterpolateV2S()->getC(Y,j, 0);
		//}

		//// Dirichlet boundary conditions in Z
		//if( 0<space()->getBCLocal()->getBCL(Z) ) {
			//OT k=space()->sInd(S,Z);
			//for( OT j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j )
				//for( OT i=space()->sInd(S,X); i<=space()->eInd(S,X); ++i )
					//y.at(i,j,k) -=
						//div_->getC(Z,k,-1) * grad_->getC( Z, k, 0) * x.at(i,j,k) *
						//space()->getInterpolateV2S()->getC(Z,k,0) /
						//space()->getInterpolateV2S()->getC(Z,k,-1);
		//}
		//if( 0<space()->getBCLocal()->getBCU(Z) ) {
			//OT k=space()->eInd(S,Z);
			//for( OT j=space()->sInd(S,Y); j<=space()->eInd(S,Y); ++j )
				//for( OT i=space()->sInd(S,X); i<=space()->eInd(S,X); ++i )
					//y.at(i,j,k) -=
						//div_->getC(Z,j,1) * grad_->getC( Z, k, 0) * x.at(i,j,k) *
						//space()->getInterpolateV2S()->getC(Z,k,-1) /
						//space()->getInterpolateV2S()->getC(Z,k, 0);
		//}

		y.changed();
	}

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
