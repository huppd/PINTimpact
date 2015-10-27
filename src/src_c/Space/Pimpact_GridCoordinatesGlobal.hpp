#pragma once
#ifndef PIMPACT_GRIDCOORDINATESGLOBAL_HPP
#define PIMPACT_GRIDCOORDINATESGLOBAL_HPP


#include <cmath>
#include <ostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_DomainSize.hpp"
#include "Pimpact_GridSizeGlobal.hpp"
#include "Pimpact_Types.hpp"




namespace Pimpact{



extern "C" {

/// \todo remove me
void PI_getGlobalCoordinates(
		const int& stretchType,
		const double& L,
		const int& M,
		const double& y_origin,
		double* const ys,
		double* const yv,
		double* const dys,
		double* const dyv );

}


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
/// \realtes GridCoordinatesGlobal
/// \ingroup SpaceObject
template<class Scalar, class Ordinal, int dim>
class GridCoordinatesGlobal {


	template<class ST,class OT,int dT>
	friend Teuchos::RCP<const GridCoordinatesGlobal<ST,OT,dT> > createGridCoordinatesGlobal(
			const Teuchos::RCP<const GridSizeGlobal<OT,dT> >& gridSize,
			const Teuchos::RCP<const DomainSize<ST> >& domainSize,
			const Teuchos::Tuple<EGridStretching,3>& gridStretching );

public:

	typedef const Teuchos::Tuple<Scalar*,dim> TO;

protected:

	Teuchos::RCP<const GridSizeGlobal<Ordinal,dim> > gridSize_;

  TO xS_;
  TO xV_;

  TO dxS_;
  TO dxV_;

	GridCoordinatesGlobal(
			const Teuchos::RCP<const GridSizeGlobal<Ordinal,dim> >& gridSize,
			const Teuchos::RCP<const DomainSize<Scalar> >& domainSize,
			const Teuchos::Tuple<EGridStretching,3>& gridStretching ):
		gridSize_( gridSize ) {

			for( int dir=0; dir<dim; ++dir ) {

				Ordinal M = gridSize_->get(dir);
				Scalar Ms = M;

				xS_[dir]  = new Scalar[ M   ];
				dxS_[dir] = new Scalar[ M   ];
				xV_[dir]  = new Scalar[ M+1 ];
				dxV_[dir] = new Scalar[ M+1 ];

				if( dir<3 ) {

					Scalar L  = domainSize->getSize(dir);
					Scalar x0 = domainSize->getOrigin(dir);

					for( Ordinal i=0; i<M; ++i ) {
						Scalar is = i;
						xS_ [dir][i] = is*L/( Ms-1. ) - x0;
						dxS_[dir][i] =    L/( Ms-1. );
					}
					for( Ordinal i=0; i<=M; ++i ) {
						Scalar is = i;
						xV_ [dir][i] = ( is - 0.5 )*L/( Ms-1) - x0;
						dxV_[dir][i] =              L/( Ms-1);
					}

				}
				else if( 3==dir ) {

					Scalar L = 4.*std::atan(1.);

					for( Ordinal i=0; i<M; ++i ) {
						Scalar is = i;
						xS_ [dir][i] = is*L/( Ms-1. );
						dxS_[dir][i] =    L/( Ms-1. );
					}
					for( Ordinal i=0; i<=M; ++i ) {
						Scalar is = i;
						xV_ [dir][i] = ( is - 0.5 )*L/( Ms-1. );
						dxV_[dir][i] =              L/( Ms-1. );
					}

				}
			}
		};

public:

	~GridCoordinatesGlobal() {
		for( int i=0; i<dim; ++i ) {
			delete[] xS_[i];
			delete[] xV_[i] ;
			delete[] dxS_[i];
			delete[] dxV_[i];
		}
	};

	/// \name getter
	/// @{ 

  /// \todo include getdx
	const Scalar* getX( ECoord dir, EField ftype ) const  {
		if( EField::S==ftype )
			return( xS_[static_cast<int>(dir)] );
		else if( static_cast<int>(dir)==static_cast<int>(ftype) )
			return( xV_[ static_cast<int>(dir) ] );
		else
			return( xS_[ static_cast<int>(dir) ] );
	}

	const Scalar* get( int dir, int ftype ) const  {
		return( getX( static_cast<ECoord>(dir), static_cast<EField>(ftype) ) );
	}
	const Scalar* get( ECoord dir, int ftype ) const  {
		return( getX( dir, static_cast<EField>(ftype) ) );
	}
	const Scalar* get( int dir, EField ftype ) const  {
		return( getX( static_cast<ECoord>(dir), ftype ) );
	}

	///  @} 
	
	void print( std::ostream& out=std::cout ) const {
		for( int i=0; i<dim; ++i ) {
			out << "Global coordinates of scalars in dir: " << i << "\n";
			out << "i\txS\n";
			for( int j=0; j<gridSize_->get(i); ++j )
				out << j+1 << "\t" << xS_[i][j] << "\n";
		}
		for( int i=0; i<dim; ++i ) {
			out << "Global coordinates of vectors in dir: " << i << "\n";
			out << "i\txV\n";
			for( int j=0; j<gridSize_->get(i)+1; ++j )
				out << j<< "\t" << xV_[i][j] << "\n";
		}
	};

}; // end of class GridCoordinatesGlobal



/// \brief create Grid coordinates Global
/// \relates GridCoordinatesGlobal
/// \todo input enum of streching
template<class S, class O, int d>
Teuchos::RCP<const GridCoordinatesGlobal<S,O,d> > createGridCoordinatesGlobal(
		const Teuchos::RCP<const GridSizeGlobal<O,d> >& gridSize,
		const Teuchos::RCP<const DomainSize<S> >& domainSize,
		const Teuchos::Tuple<EGridStretching,3>& gridStretching ) {

	return(
			Teuchos::rcp(
				new GridCoordinatesGlobal<S,O,d>(
					gridSize,
					domainSize,
					gridStretching ) ) );

}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDCOORDINATESGLOBAL_HPP
