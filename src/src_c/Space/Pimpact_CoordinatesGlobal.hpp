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
/// \tparam ScalarT scalar type
/// \tparam OrdinalT index type
/// \tparam dim computational dimension
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
template<class ScalarT, class OrdinalT, int dim>
class CoordinatesGlobal {

	template<class ST,class OT,int dT>
	friend Teuchos::RCP<const CoordinatesGlobal<ST,OT,dT> > createCoordinatesGlobal(
			const Teuchos::RCP<const GridSizeGlobal<OT> >& gridSize,
			const Teuchos::RCP<const DomainSize<ST> >& domainSize,
			const Teuchos::Tuple< Teuchos::RCP<Teuchos::ParameterList>, 3 >& gridStretching );

	template<class ST,class OT,int dT>
	friend Teuchos::RCP<const CoordinatesGlobal<ST,OT,dT> > createCoordinatesGlobal(
			const Teuchos::RCP<const GridSizeGlobal<OT> >& gridSize,
			const Teuchos::RCP<const CoordinatesGlobal<ST,OT,dT> >& coordinates );

protected:

	using TO = const Teuchos::Tuple< Teuchos::ArrayRCP<ScalarT>, dim >;

  TO xS_;
  TO xV_;

	Teuchos::Tuple< Teuchos::RCP<Teuchos::ParameterList>, 3 > stretchPara_;

  //TO dxS_; // dh: doesn't know where needed?
  //TO dxV_; // dh: locally computed in grid size local

	/// \name coordinate stretchings
	/// @{ 

	/// \brief equidistant grid
	///
	/// \f[\mathrm{ x[i] = \frac{iL}{M-1} - x0}\f]
	///
	/// \param[in] i
	/// \param[in] L
	/// \param[in] M
	/// \param[in] x0
	/// \param[out] x
	void coord_equi( const ScalarT& i, const ScalarT& L, const ScalarT& M, const ScalarT& x0, ScalarT& x/*, ScalarT& dx*/ ) {
		x  = i*L/( M-1. ) - x0;
		//dx =   L/( M-1. );
	}
  
	/// \brief coordinate stretching for parabulas
	///
	/// \f[ \mathrm{ x[i] = L\left( \frac{ i^2 }{ (M-1)^2 } + 2\alpha \frac{i}{M-1}  \right)\frac{1}{1+2\alpha)} - x0 } \f]
	/// 
	/// \param[in] i index
	/// \param[in] L length of domain
	/// \param[in] M number of global grid points
	/// \param[in] x0 origin
	/// \param[in] alpha parameter for parabola alpha=0 very parabolic alpha>>0 equidistant
	/// \param[out] x coordinate
	void coord_parab( const ScalarT& i, const ScalarT& L, const ScalarT& M, const ScalarT& x0, const ScalarT& alpha, ScalarT& x/*, ScalarT& dx*/ ) {
		x  = L*( std::pow(i,2)/std::pow(M-1.,2) + 2.*alpha*i/(M-1.) )/(1.+2.*alpha) - x0;
		//dx = L*(      2.*(i)  /std::pow(M-1.,2) + 2.*alpha  /(M-1.) )/(1.+2.*alpha);
	}



	/// \brief cos stretching
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
	/// \f[\mathrm{ wL = \frac\pi{2(i0L + iML-1)} }\f]
	/// \f[\mathrm{ xL = \frac{\cos(wL *i0L)}{wL} }\f]
	/// \f[\mathrm{ wU = \frac\pi{2(M + i0U - iML)} }\f]
	/// \f[\mathrm{ xU = \frac{\cos(wU * i0U)}{wU} }\f]
	/// \f[\mathrm{ x[i] = i*L/(M-1) - x0 }\f]
	/// 
	/// \note
	///	- i0L >= 0., i0U >= 0., iML >= 1 and iMU <= M is already tested.
	/// - so following is satisfied wL, wU >= 0..
	/// - identical to coord_tan except of std::cos functions.
	inline void coord_cos(
			const ScalarT& i,
			const ScalarT& L,
			const ScalarT& M,
			const ScalarT& x0,
			const ScalarT& iML,
			const ScalarT& iMU,
			const ScalarT& i0L,
			const ScalarT& i0U,
			ScalarT& x/*, ScalarT& dx*/ ) {

		ScalarT wL;
		ScalarT wU;
		ScalarT xL;
		ScalarT xU;

		ScalarT pi = 4.*std::atan( 1. );

		//--- parameters for grid stretching 
		if( iML<=1. ) {
			wL = 0.;
			xL = 0.;
		}
		else {
			wL = pi/( 2.*( iML + i0L - 1. ) );
			xL = std::cos( wL*i0L )/wL;
		}

		if( iMU>=M ) {
			wU = 0.;
			xU = 0.;
		}
		else {
			wU = pi/( 2.*( M - iMU + i0U ) );
			xU = std::cos( wU*i0U )/wU;
		}

		//--- coordinates in the physical space
		if( i<iML && wL!= 0. ) {
			if( (i + i0L - 1.) < 0. )
				// mirroring of the function
				x = - xL + std::cos( wL*( i + i0L - 1. ) )/wL;
			else
				x =   xL - std::cos( wL*( i + i0L - 1. ) )/wL;
		}
		else if( i>iMU && wU!=0. ) {
			if( (M - i + i0U) < 0. )
				// mirroring of the function
				x = ( 2. - std::cos( wU*( i - i0U - M ) ) )/wU;
			else 
				x =  std::cos( wU*( i - i0U - M ) )/wU;
			x = x + xL - iML + iMU;
		}
		else
			x = i + xL - iML;


		//--- Normalization
		x *= L/( xL + xU - iML + iMU );
		x -= x0;
		//std::cout << "i: " << i << "\n"
			//<< " wl: " << wL << "xl: " << xL << ", x_i: " << x << "\n";
		//std::cout << "\n"<<  x << "\n";

	}

	///  @} 



	/// \brief helper function getting number for switch statement
	/// from name 
	/// \param[in] name input name
	/// \return according int number
	int string2int( const std::string& name ) {
		std::string lcName = name;
		std::transform(lcName.begin(), lcName.end(), lcName.begin(), ::tolower);
		if( "none" == lcName ) return( 0 );
		else if( "parabola" == lcName ) return( 1 );
		else if( "parab" == lcName ) return( 1 );
		else if( "para" == lcName ) return( 1 );
		else if( "cos" == lcName ) return( 2 );
		else {
			const bool& Stertch_Type_not_known = true; 
			TEUCHOS_TEST_FOR_EXCEPT( Stertch_Type_not_known );
		}
		return( 0 );
	}



	/// \brief 
	///
	/// \param[in] gridSize
	/// \param[in] domainSize
	/// \param[in] stretchPara
	CoordinatesGlobal(
			const Teuchos::RCP<const GridSizeGlobal<OrdinalT> >& gridSize,
			const Teuchos::RCP<const DomainSize<ScalarT> >& domainSize,
			const Teuchos::Tuple< Teuchos::RCP<Teuchos::ParameterList>, 3 >& stretchPara ):
		stretchPara_(stretchPara) {

			for( int dir=0; dir<dim; ++dir ) {

				OrdinalT M = gridSize->get(dir);
				ScalarT Ms = M;

				xS_ [dir] = Teuchos::arcp<ScalarT>( M   );
				//dxS_[dir] = Teuchos::arcp<ScalarT>( M   );
				xV_ [dir] = Teuchos::arcp<ScalarT>( M+1 );
				//dxV_[dir] = Teuchos::arcp<ScalarT>( M+1 );

				if( dir<3 ) {

					ScalarT L  = domainSize->getSize(dir);
					ScalarT x0 = domainSize->getOrigin(dir);

					int stretchType = string2int( stretchPara_[dir]->get<std::string>( "Stretch Type", "none" ) );
					for( OrdinalT i=0; i<M; ++i ) {
						ScalarT is = i;

						switch( stretchType ) {
							case 0:
								coord_equi( is, L, Ms, x0, xS_[dir][i] );
								break;
							case 1:
								coord_parab( is, L, Ms, x0, stretchPara_[dir]->get<ScalarT>( "alpha", 0.5 ), xS_[dir][i] );
								break;
							case 2:
								coord_cos(
										is+1.,
										L,
										Ms,
										x0,
										stretchPara_[dir]->get<ScalarT>( "N metr L", 1. ),
										stretchPara_[dir]->get<ScalarT>( "N metr U", M  ),
										stretchPara_[dir]->get<ScalarT>( "x0 L", 0. ),
										stretchPara_[dir]->get<ScalarT>( "x0 U", 0. ),
										xS_[dir][i] );
								break;
							default:
								coord_equi( is, L, Ms, x0, xS_[dir][i] );
								break;
						}
						//std::cout << "after i: "<< i << " x: " << xS_[dir][i] << "\n\n";
					}
					for( OrdinalT i=0; i<=M; ++i ) {
						ScalarT is = static_cast<ScalarT>(i) - 0.5;
						switch( stretchType ) {
							case 0:
								coord_equi( is, L, Ms, x0, xV_[dir][i] );
								break;
							case 1:
								coord_parab( is, L, Ms, x0, stretchPara_[dir]->get<ScalarT>( "alpha", 0.5 ), xV_[dir][i] );
								break;
							case 2:
								coord_cos(
										is+1.,
										L,
										Ms,
										x0,
										stretchPara_[dir]->get<ScalarT>( "N metr L", 1. ),
										stretchPara_[dir]->get<ScalarT>( "N metr U", M  ),
										stretchPara_[dir]->get<ScalarT>( "x0 L", 0. ),
										stretchPara_[dir]->get<ScalarT>( "x0 U", 0.  ),
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

					ScalarT L = 4.*std::atan(1.);

					for( OrdinalT i=0; i<M; ++i ) {
						ScalarT is = i;
						xS_ [dir][i] = is*L/( Ms-1. );
						//dxS_[dir][i] =    L/( Ms-1. );
					}
					for( OrdinalT i=0; i<=M; ++i ) {
						ScalarT is = static_cast<ScalarT>(i) - 0.5;
						xV_ [dir][i] = is*L/( Ms-1. );
						//dxV_[dir][i] =  L/( Ms-1. );
					}
				}
			}
		};


	/// \brief constructor from fine grid
	///
	/// \param[in] gridSizeC
	/// \param[in] coordinatesF
	CoordinatesGlobal(
			const Teuchos::RCP<const GridSizeGlobal<OrdinalT> >& gridSizeC,
			const Teuchos::RCP<const CoordinatesGlobal<ScalarT,OrdinalT,dim> >& coordinatesF ) {

		for( int dir=0; dir<dim; ++dir ) {

			OrdinalT Mc = gridSizeC->get(dir);
			OrdinalT Mf = coordinatesF->xS_[dir].size();

			if( Mc==Mf ) {
				xS_ [dir] = coordinatesF->xS_ [dir];
				xV_ [dir] = coordinatesF->xV_ [dir];
			}
			else {

				xS_ [dir] = Teuchos::arcp<ScalarT>( Mc   );
				xV_ [dir] = Teuchos::arcp<ScalarT>( Mc+1 );

				OrdinalT d = 1;

				if( dir<3 ) {
					if( Mc>1 ) // shouldn't be neccessary when strategy makes its job correct. maybe throw exception
						d = ( Mf - 1 )/( Mc - 1);
				}
				else {
					if( Mc>0 ) // shouldn't be neccessary when strategy makes its job correct. maybe throw exception
						d = Mf / Mc;
				}

				for( OrdinalT j=0; j<Mc; ++j )
					xS_[dir][j] = coordinatesF->xS_[dir][j*d];

				for( OrdinalT j=1; j<Mc; ++j )
					xV_[dir][j] = coordinatesF->xS_[dir][j*d-1];

				xV_[dir][0 ] =   xS_[dir][0   ]-xV_[dir][1   ];
				xV_[dir][Mc] = 2*xS_[dir][Mc-1]-xV_[dir][Mc-1];

			}
		}

		//Teuchos::RCP<std::ostream> out = createOstream( "coord.txt", 0 );
		//print( *out );

	};


public:


	/// \name getter
	/// @{ 


	constexpr const ScalarT* getX( const int& dir, const int& ftype ) const {
		return(
				( EField::S==static_cast<EField>(ftype) || dir!=ftype )?
					xS_[dir].getRawPtr():
					xV_[dir].getRawPtr() 
				);
	}

	constexpr const Teuchos::Tuple< Teuchos::RCP<Teuchos::ParameterList> ,3>& getStretchParameter() const {
		return( stretchPara_ );
	}

	///  @} 
	
	void print( std::ostream& out=std::cout ) const {

		for( int i=0; i<dim; ++i ) {
			out << "Global coordinates of scalars in dir: " << i << "\n";
			out << "i\txS\n";
			OrdinalT j = 0;
			for( typename Teuchos::ArrayRCP<ScalarT>::iterator jp=xS_[i].begin(); jp<xS_[i].end(); ++jp )
				out << ++j << "\t" << *jp << "\n";
		}

		for( int i=0; i<dim; ++i ) {
			out << "Global coordinates of vectors in dir: " << i << "\n";
			out << "i\txV\n";
			OrdinalT j = 0;
			for( typename Teuchos::ArrayRCP<ScalarT>::iterator jp=xV_[i].begin(); jp<xV_[i].end(); ++jp )
				out << j++ << "\t" << *jp << "\n";
		}

	};


}; // end of class CoordinatesGlobal



/// \brief create Grid coordinates Global
/// \relates CoordinatesGlobal
///
/// \tparam ST
/// \tparam OT
/// \tparam d
/// \param gridSize
/// \param domainSize
/// \param gridStretching
///
/// \return 
template<class ST, class OT, int d>
Teuchos::RCP<const CoordinatesGlobal<ST,OT,d> >
createCoordinatesGlobal(
		const Teuchos::RCP<const GridSizeGlobal<OT> >& gridSize,
		const Teuchos::RCP<const DomainSize<ST> >& domainSize,
		const Teuchos::Tuple< Teuchos::RCP<Teuchos::ParameterList> ,3>& gridStretching ) {

	return(
			Teuchos::rcp(
				new CoordinatesGlobal<ST,OT,d>(
					gridSize,
					domainSize,
					gridStretching ) ) );

}



/// \brief creates coarse coordinates
///
/// \tparam ST
/// \tparam OT
/// \tparam d
/// \param gridSize
/// \param coordinates
///
/// \return 
template<class ST, class OT, int d>
Teuchos::RCP<const CoordinatesGlobal<ST,OT,d> >
createCoordinatesGlobal(
		const Teuchos::RCP<const GridSizeGlobal<OT> >& gridSize,
		const Teuchos::RCP<const CoordinatesGlobal<ST,OT,d> >& coordinates ) {

	return(
			Teuchos::rcp(
				new CoordinatesGlobal<ST,OT,d>(
					gridSize,
					coordinates )
				)
			);

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_COORDINATESGLOBAL_HPP
