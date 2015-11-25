#pragma once
#ifndef PIMPACT_COORDINATESGLOBAL_HPP
#define PIMPACT_COORDINATESGLOBAL_HPP


#include <cmath>
#include <ostream>

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_DomainSize.hpp"
#include "Pimpact_GridSizeGlobal.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact{



/// \brief global grid coordinates
///
/// Coordinate | index   | domain
/// -----------| --------| -------
/// xS         | 1..nGlo | 0..L
/// xV         | 0..nGlo | 0.5..L+0.5
///
/// \tparam Scalar scalar type
/// \tparam Ordinal index type
/// \tparam dim computiational dimension
///
/// \note: - local processor-block coordinates and grid spacings are
///          automatically derived from global grid
///        - dy3p = dy3w = 1. for 2D (may simplify programming)
///        - ensure that for all i
///             y1p(i) < y1p(i+1)
///             y1u(i) < y1u(i+1)
///             y1p(i) < y1u(i) < y1p(i+1)
///                etc.
///        - code is tested only for
///             y1p(1 ) = 0.
///             y1p(M1) = L1
///                etc.
///
/// \relates CoordinatesGlobal
/// \ingroup SpaceObject
template<class Scalar, class Ordinal, int dim>
class CoordinatesGlobal {

	template<class ST,class OT,int dT>
	friend Teuchos::RCP<const CoordinatesGlobal<ST,OT,dT> > createCoordinatesGlobal(
			const Teuchos::RCP<const GridSizeGlobal<OT> >& gridSize,
			const Teuchos::RCP<const DomainSize<ST> >& domainSize,
			const Teuchos::Tuple<Teuchos::ParameterList,3>& gridStretching );

	template<class ST,class OT,int dT>
	friend Teuchos::RCP<const CoordinatesGlobal<ST,OT,dT> > createCoordinatesGlobal(
			const Teuchos::RCP<const GridSizeGlobal<OT> >& gridSize,
			const Teuchos::RCP<const CoordinatesGlobal<ST,OT,dT> >& coordinates );

protected:

	typedef const Teuchos::Tuple< Teuchos::ArrayRCP<Scalar>, dim > TO;

  TO xS_;
  TO xV_;

	Teuchos::Tuple<Teuchos::ParameterList,3> stretchPara_;

  //TO dxS_; // dh: doesnot where needed
  //TO dxV_;

	/// \name coordinate stretchings
	/// @{ 

	/// \brief 
	///
	/// \param i
	/// \param L
	/// \param M
	/// \param x0
	/// \param x
	void coord_equi( const Scalar& i, const Scalar& L, const Scalar& M, const Scalar& x0, Scalar& x/*, Scalar& dx*/ ) {
		x  = i*L/( M-1. ) - x0;
		//dx =   L/( M-1. );
	}
  
	/// \brief coordinate stretching for parabulas
	///
	/// \param[in] i index
	/// \param[in] L length of domain
	/// \param[in] M number of global grid points
	/// \param[in] x0 origin
	/// \param bla parameter for parabola bla=0 very parabolic bla>>0 equidistant
	/// \param x coordinate
	void coord_parab( const Scalar& i, const Scalar& L, const Scalar& M, const Scalar& x0, const Scalar& bla, Scalar& x/*, Scalar& dx*/ ) {
		x  = L*( std::pow(i,2)/std::pow(M-1.,2) + 2.*bla*i/(M-1.) )/(1.+2.*bla) - x0;
		//dx = L*(      2.*(i)  /std::pow(M-1.,2) + 2.*bla  /(M-1.) )/(1.+2.*bla);
	}



	/// \brief awesome stretching
	///
	/// \param[in] i index
	/// \param[in] L length of domain
	/// \param[in] M number of global grid points
	/// \param[in] x0 origin
	/// \param[in] iML ???
	/// \param[in] iMU ???
	/// \param[in] i0L ???
	/// \param[in] i0U ???
	/// \param[out] x coordinate
	/// 
	/// \note - i0L >= 0., i0U >= 0., iML >= 1 and iMU <= M is already tested.
	///       - so following is satisfied wL, wU >= 0..
	///       - identical to coord_tan except of std::cos functions.
	void coord_cos(
			const Scalar& i,
			const Scalar& L,
			const Scalar& M,
			const Scalar& x0,
			const Scalar& iML,
			const Scalar& iMU,
			const Scalar& i0L,
			const Scalar& i0U,
			Scalar& x/*, Scalar& dx*/ ) {

		Scalar wL;
		Scalar wU;
		Scalar xL;
		Scalar xU;

		Scalar pi = 4.*std::atan( 1. );

		//--- parameters for grid stretching 
		if( iML == 1. ) {
			wL = 0.;
			xL = 0.;
		}
		else {
			wL = pi/( 2.*( iML + i0L - 1. ) );
			xL = std::cos( wL*i0L )/wL;
		}

		if( iMU == M ) {
			wU = 0.;
			xU = 0.;
		}
		else {
			wU = pi/( 2.*( M - iMU + i0U ) );
			xU = std::cos( wU*i0U )/wU;
		}

		//--- coordinates in the physical space
		if( i<iML && wL!= 0. ) {
			if( (i + i0L - 1.) < 0. ) {
				// mirroring of the function
				x = - xL + std::cos( wL*( i + i0L - 1. ) )/wL;
				//dx =     - std::sin( wL*( i + i0L - 1.   ) );
			}
			else {
				x =   xL - std::cos( wL*( i + i0L - 1. ) )/wL;
				//dx =     + std::sin( wL*( i + i0L - 1. ) );
			}
		}
		else if( i>iMU && wU!=0. ) {
			if( (M - i + i0U) < 0. ) {
				// mirroring of the function
				x = ( 2. - std::cos( wU*( i - i0U - M ) ) )/wU;
				//dx =     + std::sin( wU*( i - i0U - M ) );
			}
			else {
				x =  std::cos( wU*( i - i0U - M ) )/wU;
				//dx = - std::sin(wU*(i - i0U - M));
			}
			x = x + xL - iML + iMU;
		}
		else {
			x = i + xL - iML;
			//dx = 1.;
		}


		//--- Normalization
		x *= L/( xL + xU - iML + iMU );
		x -= x0;
		//dx = dx * L / (xL + xU - iML + iMU);

	}


	///  @} 



	/// \brief helper function getting number for switch statement
	/// from name 
	/// \param name input name
	/// \return accroding int number
	int string2int( const std::string& name ) {
		std::string lcName = name;
		std::transform(lcName.begin(), lcName.end(), lcName.begin(), ::tolower);
		if( "none" == lcName ) return( 0 );
		else if( "parabola" == lcName ) return( 1 );
		else if( "parab" == lcName ) return( 1 );
		else if( "cos" == lcName ) return( 2 );
		else {
			const bool& Stertch_Type_not_known = true; 
			TEUCHOS_TEST_FOR_EXCEPT( Stertch_Type_not_known );
		}
		return( 0 );
	}



	/// \brief 
	///
	/// \param gridSize
	/// \param domainSize
	/// \param stretchPara
	CoordinatesGlobal(
			const Teuchos::RCP<const GridSizeGlobal<Ordinal> >& gridSize,
			const Teuchos::RCP<const DomainSize<Scalar> >& domainSize,
			const Teuchos::Tuple< Teuchos::ParameterList, 3 >& stretchPara ):
		stretchPara_(stretchPara) {

			for( int dir=0; dir<dim; ++dir ) {

				Ordinal M = gridSize->get(dir);
				Scalar Ms = M;

				xS_ [dir] = Teuchos::arcp<Scalar>( M   );
				//dxS_[dir] = Teuchos::arcp<Scalar>( M   );
				xV_ [dir] = Teuchos::arcp<Scalar>( M+1 );
				//dxV_[dir] = Teuchos::arcp<Scalar>( M+1 );

				if( dir<3 ) {

					Scalar L  = domainSize->getSize(dir);
					Scalar x0 = domainSize->getOrigin(dir);

					int stretchType = string2int( stretchPara_[dir].get<std::string>( "Stretch Type", "none" ) );
					for( Ordinal i=0; i<M; ++i ) {
						Scalar is = i;

						switch( stretchType ) {
							case 0:
								coord_equi( is, L, Ms, x0, xS_[dir][i] );
								break;
							case 1:
								coord_parab( is, L, Ms, x0, stretchPara_[dir].get<Scalar>( "alpha", 0.5 ), xS_[dir][i] );
								break;
							case 2:
								coord_cos(
										is,
										L,
										Ms,
										x0,
										stretchPara_[dir].get<Scalar>( "N metr L", 1. ),
										stretchPara_[dir].get<Scalar>( "N metr U", M  ),
										stretchPara_[dir].get<Scalar>( "x0 L", 0. ),
										stretchPara_[dir].get<Scalar>( "x0 U", 0. ),
										xS_[dir][i] );
								break;
							default:
								coord_equi( is, L, Ms, x0, xS_[dir][i] );
								break;
						}
					}
					for( Ordinal i=0; i<=M; ++i ) {
						Scalar is = static_cast<Scalar>(i) - 0.5;
						switch( stretchType ) {
							case 0:
								coord_equi( is, L, Ms, x0, xV_[dir][i] );
								break;
							case 1:
								coord_parab( is, L, Ms, x0, stretchPara_[dir].get<Scalar>( "alpha", 0.5 ), xV_[dir][i] );
								break;
							case 2:
								coord_cos(
										is,
										L,
										Ms,
										x0,
										stretchPara_[dir].get<Scalar>( "N metr L", 1. ),
										stretchPara_[dir].get<Scalar>( "N metr U", M  ),
										stretchPara_[dir].get<Scalar>( "x0 L", 0. ),
										stretchPara_[dir].get<Scalar>( "x0 U", 0.  ),
										xV_[dir][i] );
								break;
							default:
								coord_equi( is, L, Ms, x0, xV_[dir][i] );
								break;
						}
						//coord_equi( is-0.5, L, Ms, x0, xV_[dir][i] );
						//coord_parab( is-0.5, L, Ms, x0, 0.05, xV_[dir][i] );
					}

				}
				else if( 3==dir ) {
					// in time direction no stretching is considered as long as
					// time-periodic problems are considered equidistant should be best

					Scalar L = 4.*std::atan(1.);

					for( Ordinal i=0; i<M; ++i ) {
						Scalar is = i;
						xS_ [dir][i] = is*L/( Ms-1. );
						//dxS_[dir][i] =    L/( Ms-1. );
					}
					for( Ordinal i=0; i<=M; ++i ) {
						Scalar is = static_cast<Scalar>(i) - 0.5;
						xV_ [dir][i] = is*L/( Ms-1. );
						//dxV_[dir][i] =  L/( Ms-1. );
					}
				}
			}
		};


	CoordinatesGlobal(
			const Teuchos::RCP<const GridSizeGlobal<Ordinal> >& gridSizeC,
			const Teuchos::RCP<const CoordinatesGlobal<Scalar,Ordinal,dim> >& coordinatesF ) {

		for( int dir=0; dir<dim; ++dir ) {

			Ordinal Mc = gridSizeC->get(dir);
			Ordinal Mf = coordinatesF->xS_[dir].size();

			if( Mc==Mf ) {
				xS_ [dir] = coordinatesF->xS_ [dir];
				xV_ [dir] = coordinatesF->xV_ [dir];
			}
			else {

				xS_ [dir] = Teuchos::arcp<Scalar>( Mc   );
				xV_ [dir] = Teuchos::arcp<Scalar>( Mc+1 );

				Ordinal d = 1;

				if( dir<3 ) {
					if( Mc>1 ) // shouldn't be neccessary when strategy makes its job correct. maybe throw exception
						d = ( Mf - 1 )/( Mc - 1);
				}
				else {
					if( Mc>0 ) // shouldn't be neccessary when strategy makes its job correct. maybe throw exception
						d = Mf / Mc;
				}

				for( Ordinal j=0; j<Mc; ++j )
					xS_[dir][j] = coordinatesF->xS_[dir][j*d];

				for( Ordinal j=1; j<Mc; ++j )
					xV_[dir][j] = coordinatesF->xS_[dir][j*d-1];

				xV_[dir][0 ] =   xS_[dir][0   ]-xV_[dir][1   ];
				xV_[dir][Mc] = 2*xS_[dir][Mc-1]-xV_[dir][Mc-1];

			}
		}
	};


public:


	/// \name getter
	/// @{ 

	const Scalar* getX( ECoord dir, EField ftype ) const {
		if( EField::S==ftype )
			return( xS_[static_cast<int>(dir)].getRawPtr() );
		else if( static_cast<int>(dir)==static_cast<int>(ftype) )
			return( xV_[ static_cast<int>(dir) ].getRawPtr() );
		else
			return( xS_[ static_cast<int>(dir) ].getRawPtr() );
	}

	const Scalar* get( int dir, int ftype ) const {
		return( getX( static_cast<ECoord>(dir), static_cast<EField>(ftype) ) );
	}
	const Scalar* get( ECoord dir, int ftype ) const {
		return( getX( dir, static_cast<EField>(ftype) ) );
	}
	const Scalar* get( int dir, EField ftype ) const {
		return( getX( static_cast<ECoord>(dir), ftype ) );
	}


	const Teuchos::Tuple<Teuchos::ParameterList,3>& getStretchParameter() const {
		return( stretchPara_ );
	}
	///  @} 
	
	void print( std::ostream& out=std::cout ) const {

		for( int i=0; i<dim; ++i ) {
			out << "Global coordinates of scalars in dir: " << i << "\n";
			out << "i\txS\n";
			Ordinal j = 0;
			for( typename Teuchos::ArrayRCP<Scalar>::iterator jp=xS_[i].begin(); jp<xS_[i].end(); ++jp )
				out << ++j << "\t" << *jp << "\n";
		}

		for( int i=0; i<dim; ++i ) {
			out << "Global coordinates of vectors in dir: " << i << "\n";
			out << "i\txV\n";
			Ordinal j = 0;
			for( typename Teuchos::ArrayRCP<Scalar>::iterator jp=xV_[i].begin(); jp<xV_[i].end(); ++jp )
				out << j++ << "\t" << *jp << "\n";
		}

	};


}; // end of class CoordinatesGlobal



/// \brief create Grid coordinates Global
/// \relates CoordinatesGlobal
///
/// \tparam S
/// \tparam O
/// \tparam d
/// \param gridSize
/// \param domainSize
/// \param gridStretching
///
/// \return 
template<class S, class O, int d>
Teuchos::RCP<const CoordinatesGlobal<S,O,d> >
createCoordinatesGlobal(
		const Teuchos::RCP<const GridSizeGlobal<O> >& gridSize,
		const Teuchos::RCP<const DomainSize<S> >& domainSize,
		const Teuchos::Tuple<Teuchos::ParameterList,3>& gridStretching ) {

	return(
			Teuchos::rcp(
				new CoordinatesGlobal<S,O,d>(
					gridSize,
					domainSize,
					gridStretching ) ) );

}



/// \brief creates coarse coordinates
///
/// \tparam S
/// \tparam O
/// \tparam d
/// \param gridSize
/// \param coordinates
///
/// \return 
template<class S, class O, int d>
Teuchos::RCP<const CoordinatesGlobal<S,O,d> >
createCoordinatesGlobal(
		const Teuchos::RCP<const GridSizeGlobal<O> >& gridSize,
		const Teuchos::RCP<const CoordinatesGlobal<S,O,d> >& coordinates ) {

	return(
			Teuchos::rcp(
				new CoordinatesGlobal<S,O,d>(
					gridSize,
					coordinates )
				)
			);

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COORDINATESGLOBAL_HPP