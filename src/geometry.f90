!======================================================================!
!
       module geometry
       implicit none
!
       contains
!
!======================================================================!
!
       subroutine Rz(n,inpmat,angle,outmat)
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,n),intent(in)   ::  inpmat  !  Input matrix
       real(kind=8),dimension(3,n),intent(out)  ::  outmat  !  Output matrix
       real(kind=8),intent(in)                  ::  angle   !  Rotation angle
       integer,intent(in)                       ::  n       !  Dimension of the input matrix
!
! Declaration of the local variables
!
       real(kind=8),dimension(3,3)              ::  rota      !  Rotation matrix
!
! Rotating around z-axis
!
       rota(:,:) =  0.0d0
       rota(1,1) =  cos(angle)
       rota(2,2) =  cos(angle)
       rota(3,3) =  1.0d0
       rota(1,2) = -sin(angle)
       rota(2,1) =  sin(angle)
!
       outmat = matmul(rota,inpmat)
!
       return
       end subroutine Rz
!
!======================================================================!
!
       subroutine Ry(n,inpmat,angle,outmat)
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,n),intent(in)   ::  inpmat  !  Input matrix
       real(kind=8),dimension(3,n),intent(out)  ::  outmat  !  Output matrix
       real(kind=8),intent(in)                  ::  angle   !  Rotation angle
       integer,intent(in)                       ::  n       !  Dimension of the input matrix
!
! Declaration of the local variables
!
       real(kind=8),dimension(3,3)              ::  rota      !  Rotation matrix
!
! Rotating around y-axis
!
       rota(:,:) =  0.0d0
       rota(1,1) =  cos(angle)
       rota(2,2) =  1.0d0
       rota(3,3) =  cos(angle)
       rota(3,1) = -sin(angle)
       rota(1,3) =  sin(angle)
!
       outmat = matmul(rota,inpmat)
!
       return
       end subroutine Ry
!
!======================================================================!
!
       subroutine Rx(n,inpmat,angle,outmat)
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,n),intent(in)   ::  inpmat  !  Input matrix
       real(kind=8),dimension(3,n),intent(out)  ::  outmat  !  Output matrix
       real(kind=8),intent(in)                  ::  angle   !  Rotation angle
       integer,intent(in)                       ::  n       !  Dimension of the input matrix
!
! Declaration of the local variables
!
       real(kind=8),dimension(3,3)              ::  rota      !  Rotation matrix
!
! Rotating around x-axis
!
       rota(:,:) =  0.0d0
       rota(1,1) =  1.0d0
       rota(2,2) =  cos(angle)
       rota(3,3) =  cos(angle)
       rota(2,3) = -sin(angle)
       rota(3,2) =  sin(angle)
!
       outmat = matmul(rota,inpmat)
!
       return
       end subroutine Rx
!
!======================================================================!
!
       subroutine Syz(n,inpmat,outmat)
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,n),intent(in)   ::  inpmat  !  Input matrix
       real(kind=8),dimension(3,n),intent(out)  ::  outmat  !  Output matrix
       integer,intent(in)                       ::  n       !  Dimension of the input matrix
!
! Declaration of the local variables
!
       real(kind=8),dimension(3,3)              ::  rota      !  Rotation matrix
!
! Reflection in the x-axis
!
       rota(:,:) =  0.0d0
       rota(1,1) = -1.0d0
       rota(2,2) =  1.0d0
       rota(3,3) =  0.0d0
!
       outmat = matmul(rota,inpmat)
!
       return
       end subroutine Syz
!
!======================================================================!
!
       subroutine Sxz(n,inpmat,outmat)
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,n),intent(in)   ::  inpmat  !  Input matrix
       real(kind=8),dimension(3,n),intent(out)  ::  outmat  !  Output matrix
       integer,intent(in)                       ::  n       !  Dimension of the input matrix
!
! Declaration of the local variables
!
       real(kind=8),dimension(3,3)              ::  rota      !  Rotation matrix
!
! Reflection in the y-axis
!
       rota(:,:) =  0.0d0
       rota(1,1) =  1.0d0
       rota(2,2) = -1.0d0
       rota(3,3) =  0.0d0
!
       outmat = matmul(rota,inpmat)
!
       return
       end subroutine Sxz
!
!======================================================================!
!
       subroutine Sxy(n,inpmat,outmat)
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,n),intent(in)   ::  inpmat  !  Input matrix
       real(kind=8),dimension(3,n),intent(out)  ::  outmat  !  Output matrix
       integer,intent(in)                       ::  n       !  Dimension of the input matrix
!
! Declaration of the local variables
!
       real(kind=8),dimension(3,3)              ::  rota      !  Rotation matrix
!
! Reflection in the z-axis
!
       rota(:,:) =  0.0d0
       rota(1,1) =  1.0d0
       rota(2,2) =  1.0d0
       rota(3,3) = -1.0d0
!
       outmat = matmul(rota,inpmat)
!
       return
       end subroutine Sxy
!
!======================================================================!
!
       subroutine inv(n,inpmat,outmat)
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,n),intent(in)   ::  inpmat  !  Input matrix
       real(kind=8),dimension(3,n),intent(out)  ::  outmat  !  Output matrix
       integer,intent(in)                       ::  n       !  Dimension of the input matrix
!
! Declaration of the local variables
!
       integer                                  ::  i,j     !  Indexes
!
! Inversion of coordinates
!
       do i = 1, n
         do j = 1, 3
           outmat(j,i) = -inpmat(j,i)
         end do
       end do
!
       return
       end subroutine inv
!
!======================================================================!
!
       subroutine translate(m,n,inpmat,alpha,vec,outmat)
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(m,n),intent(in)   ::  inpmat  !  Input matrix
       real(kind=8),dimension(m,n),intent(out)  ::  outmat  !  Output matrix
       real(kind=8),dimension(m),intent(in)     ::  vec     !  Rotation angle
       real(kind=8),intent(in)                  ::  alpha   !  Direction
       integer,intent(in)                       ::  m       !  Real space dimension
       integer,intent(in)                       ::  n       !  Number of points
!
! Declaration of the local variables
!
       integer                                  ::  i,j     !  Indexes     
!
! Translating though space
!
       do i = 1, n
         do j = 1, m
           outmat(j,i) = inpmat(j,i)  + alpha*vec(j)
         end do
       end do
!
       return
       end subroutine translate
!
!======================================================================!
!
! DCOMVEC - Double precision Center Of Mass VECtor
!
! This function returns the Center of Mass of a molecule
!
       function dcomvec(nat,coord,mass) result(com)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,nat),intent(in)  ::  coord  !  Coordinates
       real(kind=8),dimension(nat),intent(in)    ::  mass   !  Masses
       real(kind=8), dimension(3)                ::  com    !  Center of Mass vector
       integer,intent(in)                        ::  nat    !  Number of atoms
!
! Declaration of the local variables
!
       real(kind=8)                              ::  totm   ! Total mass
       integer                                   ::  i,j    ! Indexes
!
! Calculating the center of mass
!
       totm = 0.0d0
       do i = 1, nat
         totm = totm + mass(i)
       end do
!
       com(:) = 0.0d0
       do i = 1, nat
         do j = 1, 3
           com(j) = com(j) + mass(i)*coord(j,i)
         end do
       end do
!
       com(:) = com(:) / totm
!
       return
       end function dcomvec
!
!======================================================================!
!
! SCOMVEC - Single precision Center Of Mass VECtor
!
! This function returns the Center of Mass of a molecule
!
       function scomvec(nat,coord,mass) result(com)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=4),dimension(3,nat),intent(in)  ::  coord  !  Coordinates
       real(kind=8),dimension(nat),intent(in)    ::  mass   !  Masses
       real(kind=8), dimension(3)                ::  com    !  Center of Mass vector
       integer,intent(in)                        ::  nat    !  Number of atoms
!
! Declaration of the local variables
!
       real(kind=8)                              ::  totm   ! Total mass
       integer                                   ::  i,j    ! Indexes
!
! Calculating the center of mass
!
       totm = 0.0d0
       do i = 1, nat
         totm = totm + mass(i)
       end do
!
       com(:) = 0.0d0
       do i = 1, nat
         do j = 1, 3
           com(j) = com(j) + mass(i)*coord(j,i)
         end do
       end do
!
       com(:) = com(:) / totm
!
       return
       end function scomvec
!
!======================================================================!
!
! DCENVEC - Double precision CENtroid VECtor
!
! This function returns the centroid of a set of points
!
       function dcenvec(n,m,point) result(cen)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(n,m),intent(in)  ::  point  !  Coordinates
       real(kind=8), dimension(n)              ::  cen    !  Center of Mass vector
       integer,intent(in)                      ::  n      !  Space dimension
       integer,intent(in)                      ::  m      !  Number of points
!
! Declaration of the local variables
!
       integer                                 ::  i,j    ! Indexes
!
! Calculating the centroid of the set of points
!
       cen(:) = 0.0d0
       do i = 1, m
         do j = 1, n
           cen(j) = cen(j) + point(j,i)
         end do
       end do
!
       cen(:) = cen(:)/m
!
       return
       end function dcenvec
!
!======================================================================!
!
! SCENVEC - Single precision CENtroid VECtor
!
! This function returns the centroid of a set of points
!
       function scenvec(n,m,point) result(cen)
!
        implicit none
!
! Declaration of the in/out variables
!
        real(kind=4),dimension(n,m),intent(in)  ::  point  !  Coordinates
        real(kind=4), dimension(n)              ::  cen    !  Center of Mass vector
        integer,intent(in)                      ::  n      !  Space dimension
        integer,intent(in)                      ::  m      !  Number of points
!
! Declaration of the local variables
!
       integer                                 ::  i,j    ! Indexes
!
! Calculating the centroid of the set of points
!
       cen(:) = 0.0d0
       do i = 1, m
         do j = 1, n
           cen(j) = cen(j) + point(j,i)
         end do
       end do
!
       cen(:) = cen(:)/m
!
       return
       end function scenvec
!
!======================================================================!
!
! DMINIMGVEC - Double precision MINimum IMaGe VECtor
!
! This function returns
!
       function dminimgvec(v1,v2,L) result(r)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3),intent(in)  ::  v1  !
       real(kind=8),dimension(3),intent(in)  ::  v2  !
       real(kind=8),dimension(3),intent(in)  ::  L   !
       real(kind=8),dimension(3)             ::  r   !
!
! Calculating the minimum image vector
!
       r(:)  = v2(:) - v1(:)
!
       r(:)  = r(:) - L(:)*anint(r/L)
!	     
       return
       end function dminimgvec
!
!======================================================================!
!
! SMINIMGVEC - Single precision MINimum IMaGe VECtor
!
! This function returns
!
       function sminimgvec(v1,v2,L) result(r)
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=4),dimension(3),intent(in)  ::  v1  !
       real(kind=4),dimension(3),intent(in)  ::  v2  !
       real(kind=4),dimension(3),intent(in)  ::  L   !
       real(kind=4),dimension(3)             ::  r   !
!
! Calculating the minimum image vector
!
       r(:) = v2(:) - v1(:)
!
       r(:) = r(:) - L(:)*anint(r/L)
!	     
       return
       end function sminimgvec
!
!======================================================================!
!
       end module geometry
!
!======================================================================!
