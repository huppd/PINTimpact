#pragma once
#ifndef PIMPACT_GRIDCOORDINATESGLOBAL_HPP
#define PIMPACT_GRIDCOORDINATESGLOBAL_HPP



#include<cmath>
#include<ostream>


#include"Teuchos_RCP.hpp"
#include"Teuchos_Tuple.hpp"

#include"Pimpact_Types.hpp"




namespace Pimpact{


extern "C" {

void PI_getGlobalCoordinates(
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
/// \ingroup SpaceObject
template<class Scalar, class Ordinal, int dim>
class GridCoordinatesGlobal {


	template<class ST,class OT,int dT>
	friend Teuchos::RCP<const GridCoordinatesGlobal<ST,OT,dT> > createGridCoordinatesGlobal(
			const Teuchos::RCP<const GridSizeGlobal<OT,dT> >& gridSize,
			const Teuchos::RCP<const DomainSize<ST> >& domainSize );

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
			const Teuchos::RCP<const DomainSize<Scalar> >& domainSize ):
		gridSize_( gridSize ) {

			for( int i=0; i<dim; ++i ) {
				xS_[i]  = new Scalar[ gridSize_->get(i)   ];
				xV_[i]  = new Scalar[ gridSize_->get(i)+1 ];
				dxS_[i] = new Scalar[ gridSize_->get(i)   ];
				dxV_[i] = new Scalar[ gridSize_->get(i)+1 ];
				if( i<3 )
					PI_getGlobalCoordinates(
							domainSize->getSize(i),
							gridSize_->get(i),
							domainSize->getOrigin(i),
							xS_[i],
							xV_[i],
							dxS_[i],
							dxV_[i] );
				else if( 3==i )
					PI_getGlobalCoordinates(
							4.*std::atan(1.),
							gridSize_->get(i),
							0.,
							xS_[i],
							xV_[i],
							dxS_[i],
							dxV_[i] );
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
			out << "Global coordinates of scalars in dir: " << i << "\n";
			out << "i\txS\n";
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
		const Teuchos::RCP<const DomainSize<S> >& domainSize ) {

	return(
			Teuchos::rcp(
				new GridCoordinatesGlobal<S,O,d>( gridSize, domainSize ) ) );

}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_GRIDCOORDINATESGLOBAL_HPP
