!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !--- grid, index and index-range convention---
  !
  ! coordinates and indices on staggered grid (2D):
  !   - boundary conditions are set on pressure grid points
  !   - assignments to grid points outside of the given index ranges are ignored or overwritten
  !   - x1p, x1v, x2p, x2v, x3p, x3w: local  coordinates (inside each processor-block)
  !   - y1p, y1v, y2p, y2v, y3p, y3w: global coordinates (include all blocks)
  !
  !
  !    O: pressure, concentrations
  !    >: vel(:,:,:1)
  !    v: vel(:,:,:2)
  !   //: boundary
  !
  !               y1p(1)  y1p(2)   ...    y1p(M1)
  !           y1u(0)  y1u(1)  y1u(2)   ...    y1u(M1)
  !
  ! x2v(N22B)       v       v       v       v       y2v(M2)
  !               //|///////|///////|///////|//
  ! x2p(N2p)    >---O--->---O--->---O--->---O--->   y2p(M2)
  !  .            //|       |       |       |//      .
  !  .            //v       v       v       v//      .
  !  .            //|       |       |       |//      .
  ! x2p(S1p +2) >---O--->---O--->---O--->---O--->   y2p(3)
  !               //|       |       |       |//
  ! x2v(S22B+2)   //v       v       v       v//     y2v(2)
  !               //|       |       |       |//
  ! x2p(S1p +1) >---O--->---O--->---O--->---O--->   y2p(2)
  !               //|       |       |       |//
  ! x2v(S22B+1)   //v       v       v       v//     y2v(1)
  !               //|       |       |       |//
  ! x2p(S1p)    >---O--->---O--->---O--->---O--->   y2p(1)
  !               //|///////|///////|///////|//
  ! x2v(S22B)       v       v       v       v       y2v(0)
  !
  !              x1p(S1p) x1p(S1p+1) ...  x1p(N1p)
  !         x1u(S11B) x1u(S11B+1) x1u(S11B+2)...x1u(N11B)
  !
  !
  !
  ! pressure / concentrations:
  !
  !  start/end index
  !  |spatial direction
  !  ||
  !  S1p
  !  N1p
  !
  !
  ! velocity:
  !
  !  start/end index
  !  |spatial direction
  !  ||velocity component
  !  |||includes/excludes the outermost grid points on the boundary or outside the domain
  !  ||||                             (periodicity/symmetry is not treated as a boundary)
  !  S11B
  !  N11B
  !  S11
  !  N11
  !
  !
  ! Note that
  !   S12B = S13B = S1p
  !   N12B = N13B = N1p
  !   S21B = S23B = S2p
  !   N21B = N23B = N2p
  !   S31B = S32B = S3p
  !   N31B = N32B = N3p
  !
  
  
  !> \brief here can the user sets the coordinates according to his favorite stretching
  subroutine get_coords
  ! (basic subroutine)
  
  use mod_dims
  use mod_vars
  use usr_vars
  use usr_func
  
  implicit none
  
  integer                ::  i, j, k
  
  
  !--- specify global coordinates and grid spacings ---
  ! note: - local processor-block coordinates and grid spacings are automatically derived from global grid
  !       - dy3p = dy3w = 1. for 2D (may simplify programming)
  !       - ensure that for all i
  !            y1p(i) < y1p(i+1)
  !            y1u(i) < y1u(i+1)
  !            y1p(i) < y1u(i) < y1p(i+1)
  !               etc.
  !       - code is tested only for
  !            y1p(1 ) = 0.
  !            y1p(M1) = L1
  !               etc.
  !
  ! example:
  !
  do i = 1, M1
     call coord_tanh(L1,REAL(M1),0.,0.,REAL(i)    ,y1p(i),dy1p(i))
  end do
  do i = 0, M1
     call coord_tanh(L1,REAL(M1),0.,0.,REAL(i)+0.5,y1u(i),dy1u(i))
  end do
  y1p = y1p - y1_origin
  y1u = y1u - y1_origin
  
  do j = 1, M2
!     CALL coord_tanh(L2,REAL(M2),0.,1.,REAL(j)    ,y2p(j),dy2p(j))
     call coord_tanh(L2,REAL(M2),0.,0.,REAL(j)    ,y2p(j),dy2p(j))
  end do
  do j = 0, M2
!     CALL coord_tanh(L2,REAL(M2),0.,1.,REAL(j)+0.5,y2v(j),dy2v(j))
     call coord_tanh(L2,REAL(M2),0.,0.,REAL(j)+0.5,y2v(j),dy2v(j))
  end do
  y2p = y2p - y2_origin
  y2v = y2v - y2_origin
  
  
  if (M3 /= 2) then
  do k = 1, M3
     call coord_tanh(L3,REAL(M3),0.,0.,REAL(k)    ,y3p(k),dy3p(k))
  end do
  do k = 0, M3
     call coord_tanh(L3,REAL(M3),0.,0.,REAL(k)+0.5,y3w(k),dy3w(k))
  end do
  y3p = y3p - y3_origin
  y3w = y3w - y3_origin
  end if
  
  
  end subroutine get_coords
