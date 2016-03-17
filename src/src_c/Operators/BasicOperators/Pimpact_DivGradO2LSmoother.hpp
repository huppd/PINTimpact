#pragma once
#ifndef PIMPACT_DIVGRADO2LSMOOTHER_HPP
#define PIMPACT_DIVGRADO2LSMOOTHER_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Pimpact_TeuchosTransfer.hpp"




namespace Pimpact{



/// \brief \f$\omega\f$-Line smoother for second Order DivGradOp.
///
///
/// \relates DivGradO2Op
/// \ingroup BaseOperator
/// \todo work in progress
template<class OperatorT>
class DivGradO2LSmoother {

public:

  using SpaceT = typename OperatorT::SpaceT;

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  using DomainFieldT = ScalarField<SpaceT>;
  using RangeFieldT = ScalarField<SpaceT>;

protected:

	using VectorT = Teuchos::SerialDenseVector<Ordinal,Scalar>;
	using OVectorT = Teuchos::SerialDenseVector<Ordinal,Ordinal>;
	using MatrixT = Teuchos::SerialDenseMatrix<Ordinal,Scalar>;

  int nIter_;
  Scalar omega_;

	Teuchos::Tuple<bool,3> lineDirection_;

	Teuchos::Tuple<Ordinal,3> n_;

	bool levelYes_;

  Teuchos::RCP<DomainFieldT> temp_;

  const Teuchos::RCP<const OperatorT> op_;

	Teuchos::Tuple< Teuchos::RCP<VectorT>, 3 > dl_;
	Teuchos::Tuple< Teuchos::RCP<VectorT>, 3 > du_;
	Teuchos::Tuple< Teuchos::RCP<VectorT>, 3 > d_;
	Teuchos::Tuple< Teuchos::RCP<VectorT>, 3 > du2_;

	Teuchos::Tuple< Teuchos::RCP<OVectorT> , 3 > ipiv_;

public:

	/// \brief constructor
	///
	/// \param[in] op pointer to operator that is smoothed
	/// \param[in] pl  Parameter list of options for the multi grid solver.
	///   These are the options accepted by the solver manager:
	///   - "omega" - a \c Scalar damping factor. Default: for 2D 0.8 for 3D 6./7.  /
	///   - "numIters" - a \c int number of smoothing steps. Default: 4  /
	///   - "BC smoothing" - a \c int type of BC smoothing 0, 0: Lacbobian, else: direct. Default: 0 /
	///   - "level" - a \c bool number of smoothing steps. Default: false  /
	DivGradO2LSmoother(
			const Teuchos::RCP<const OperatorT>& op,
      const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList() ):
    nIter_( pl->get<int>( "numIters", 1 ) ),
    omega_( pl->get<Scalar>( "omega", 0.75 ) ),
    levelYes_( pl->get<bool>( "level", false ) ),
    temp_( create<DomainFieldT>( op->space() ) ),
    op_(op) {
		
			lineDirection_[X] = pl->get<bool>( "X", true );
			lineDirection_[Y] = pl->get<bool>( "Y", true );
			lineDirection_[Z] = pl->get<bool>( "Z", true );
			
			TEUCHOS_TEST_FOR_EXCEPT( !lineDirection_[X] && !lineDirection_[Y] && !lineDirection_[Z] );

			for( int dir=0; dir<3; ++dir ) {
				n_[dir] = space()->eInd( S, dir ) - space()->sInd( S, dir ) + 1 ;
				if( true==lineDirection_[dir] ) {
					//std::cout << n_[dir] << "\n";

					dl_ [dir] = Teuchos::rcp( new VectorT( n_[dir]-1, true ) );
					d_  [dir] = Teuchos::rcp( new VectorT( n_[dir]  , true ) );
					du_ [dir] = Teuchos::rcp( new VectorT( n_[dir]-1, true ) );
					du2_[dir] = Teuchos::rcp( new VectorT( n_[dir]-2, true ) );
                                                        
					ipiv_[dir] = Teuchos::rcp( new OVectorT( n_[dir], true ) );

					for( Ordinal i=space()->sInd( S, dir ); i<=space()->eInd( S, dir ); ++i ) {
						for( int j=0; j<space()->dim(); ++j )
							(*d_[dir])[i-1] += op_->getC( j, i,  0 );
					}
					for( Ordinal i=space()->sInd( S, dir ); i<space()->eInd( S, dir ); ++i ) {
						(*du_[dir])[i-1] = op_->getC( dir, i, +1 );
						(*dl_[dir])[i-1] = op_->getC( dir, i+1, -1 );
					}

					//std::cout << "dl: " << *dl_[dir];
					//std::cout << "d: " << *d_[dir] ;
					//std::cout << "du: " << *du_[dir];

					Ordinal info;
					Teuchos::LAPACK<Ordinal,Scalar> lapack;
					lapack.GTTRF(
							n_[dir],
							dl_[dir]->values(),
							d_[dir]->values(),
							du_[dir]->values(),
							du2_[dir]->values(),
							ipiv_[dir]->values(),
							&info );
					//std::cout << "\ninfo: " << info << "\n";
				}
			}
		}


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}( x - A y_k ) \f]
	/// \todo solve multiple problems at once
	void apply(const DomainFieldT& b, RangeFieldT& y,
			Belos::ETrans trans=Belos::NOTRANS ) const {

		for( int i=0; i<nIter_; ++i) {

			//for( int dir=0; dir<3; ++dir ) {
			for( int dir=2; dir>=0; --dir ) { // why?

				if( true==lineDirection_[dir] ) {

					op_->computeResidual( b, y, *temp_ );

					Ordinal SS[3];
					Ordinal NN[3];
					Ordinal i[3];

					int d1 = ( dir + 1 )%3;
					int d2 = ( dir + 2 )%3;
					if( d2>d1 ) std::swap( d2, d1 );

					for( int j=0; j<3; ++j ) {
						SS[j] = op_->getSR(j);
						NN[j] = op_->getER(j);
					}

					Teuchos::RCP<VectorT> B = Teuchos::rcp( new VectorT( n_[dir], false ) );

					for( i[d1]=SS[d1]; i[d1]<=NN[d1]; ++i[d1] )
						for( i[d2]=SS[d2]; i[d2]<=NN[d2]; ++i[d2] ) {

							// transfer
							for( i[dir]=1; i[dir]<=n_[dir]; ++i[dir] )
								(*B)[i[dir]-1] = temp_->at( i[0], i[1], i[2] );
								//(*B)[i[dir]-1] = b.at( i[0], i[1], i[2] );
							//std::cout << *B << "\n";

							Ordinal info;
							Teuchos::LAPACK<Ordinal,Scalar> lapack;
							lapack.GTTRS(
									'N',
									n_[dir],
									1,
									dl_[dir]->values(),
									d_[dir]->values(),
									du_[dir]->values(),
									du2_[dir]->values(),
									ipiv_[dir]->values(),
									B->values(),
									B->stride(),
									&info );
							// transfer back
							for( i[dir]=1; i[dir]<=n_[dir]; ++i[dir] )
								y.at( i[0], i[1], i[2] ) = y.at( i[0], i[1], i[2] ) + omega_*(*B)[i[dir]-1];
						}
					y.changed();
				}
			}
		}

		if( levelYes_ )
			y.level();

	}

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    out << "\t numIter: " << nIter_ << "\n";
    op_->print( out );
  }

	const std::string getLabel() const { return( "DivGradO2LSmoother" ); };


}; // end of class DivGradO2LSmoother




} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::DivGradO2LSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::DivGradO2LSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::DivGradO2LSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::DivGradO2LSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,4> > >;
#endif


#endif // end of #ifndef PIMPACT_DIVGRADO2LSMOOTHER_HPP
