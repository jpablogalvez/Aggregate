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
       subroutine Sx(n,inpmat,outmat)
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
       end subroutine Sx     
!
!======================================================================!
!
       subroutine Sy(n,inpmat,outmat)
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
       end subroutine Sy
!
!======================================================================!
!
       subroutine Sz(n,inpmat,outmat)
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
       end subroutine Sz  
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
       subroutine translate(n,inpmat,alpha,vec,outmat)
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,n),intent(in)   ::  inpmat  !  Input matrix
       real(kind=8),dimension(3,n),intent(out)  ::  outmat  !  Output matrix
       real(kind=8),dimension(3),intent(in)     ::  vec     !  Rotation angle
       real(kind=8),intent(in)                  ::  alpha   !  Direction
       integer,intent(in)                       ::  n       !  Dimension of the input matrix
!
! Declaration of the local variables
!
       integer                                  ::  i,j     !  Indexes     
!
! Translating though space
!
       do i = 1, n
         do j = 1, 3
           outmat(j,i) = inpmat(j,i)  + alpha*vec(j)
         end do
       end do
!
       return
       end subroutine translate
!
!======================================================================!
!
       function cofm_vector(nat,coord,mass) result(cofm)
!
! This function returns the Center of Mass of a molecule
!
! Parameters:
!
!  nat
!  coord
!  mass
!
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=8),dimension(3,nat),intent(in)  ::  coord  !  Coordinates
       real(kind=8),dimension(nat),intent(in)    ::  mass   !  Masses
       integer,intent(in)                        ::  nat    !  Number of atoms
       real(kind=8), dimension(3)                ::  cofm   !  Center of Mass vector
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
       cofm(:) = 0.0d0
       do i = 1, nat
         do j = 1, 3
           cofm(j) = cofm(j) + mass(i)*coord(j,i)
         end do
       end do
!
       cofm(:) = cofm(:) / totm
!
       return
       end function cofm_vector
!
!======================================================================!
!
       function minimgvec(v1,v2,L) result(r)
       implicit none
!
! Declaration of the in/out variables
!
       real(kind=4),dimension(3),intent(in)  ::  v1  !
       real(kind=4),dimension(3),intent(in)  ::  v2  !
       real(kind=4),dimension(3),intent(in)  ::  L   !
       real(kind=4),dimension(3)             ::  r   !
!
! Declaration of the local variables
!
       integer                               ::  i   !
!
! Calculating the minimum image vector
!
       r  = v2 - v1
!
       r  = r - L*anint(r/L)
!	     
       return
       end function minimgvec
!
!======================================================================!
!
       end module geometry
!
!======================================================================!
